//*************************************************************************************************
//*  Filename:   cmoil_dib.cpp
//*
//*  Description: 
//*    Bimtmap operations specific for windows
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Nov. 1, 2000	Baohua Wang	Initial Release
//*  Nov. 4, 2003 Baohua Wang Use Magick++ lib
//*  Sep. 6, 2007 Thomas Blom  Update for Linux/OSX
//*************************************************************************************************
#ifdef __APPLE__
	#include <GLUT/glut.h>
#elif defined(WIN32)
	#include "GLwin32/glut.h"
#else
	#include "GL/glut.h"
#endif
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "cmoil_dib.h"
#include "cmoil_menu.h"

namespace CMOIL {

CmoilBitmap::CmoilBitmap ( ) 
{
   pixels = NULL;
   imageBytes = 0;
   imageWidth=0;
   imageHeight=0;
}

CmoilBitmap::~CmoilBitmap ( ) {
   if ( pixels != NULL)
      delete pixels;
}

//Get image as DIB, then convert to output filefomat and save
char *CmoilBitmap::PrintImg (const char*imgfilename)
{
  GetPixels( );

  // get RGBA format
  static char  geostr[20]="0x0"; 
  if ( pixels != NULL )
  {
     // switch rows since Magick::Geometry.yNegative doesn't work
     int rowsize= ((strcmp(pixelFormat, "RGB")==0)?3:4) * imageWidth;
     char *tmpcptr = new char[rowsize];

     if ( tmpcptr != NULL ) 
     {
       char *ipixels=pixels, *jpixels=pixels+imageBytes-rowsize;
       for ( int i=0; i<imageHeight/2; i++, ipixels+=rowsize, jpixels-=rowsize)
       {
          memcpy(tmpcptr, ipixels, rowsize);    // swap a row
          memcpy(ipixels, jpixels, rowsize);
          memcpy(jpixels, tmpcptr, rowsize);
       }
       delete[] tmpcptr;
     }
     sprintf(geostr, "%dx%d", imageWidth, imageHeight);

     // separate ImageMagick from cmoil code 
     int pixelSize;
     char mgkcmd[FILENAME_MAXLEN];

     sprintf(mgkcmd,"%scmoil_magick",ExeDir ? ExeDir : "" );

#ifdef _WIN32
     FILE *fp = POPEN(mgkcmd, "wb");
#else
     FILE *fp = POPEN(mgkcmd, "w");    // "wb" fails on linux
#endif

     if ( fp != NULL)
     {
       fwrite(imgfilename, FILENAME_MAXLEN, 1, fp);
       fwrite(&imageBytes, 4, 1, fp);
       fwrite(&imageWidth, 4, 1, fp);
       fwrite(&imageHeight, 4, 1, fp);
       pixelSize=strcmp(pixelFormat,"RGB")==0?3:4;
       fwrite(&pixelSize, 4, 1, fp);
       fwrite(pixels, 512, imageBytes/512+1, fp);
       PCLOSE(fp);
     }
  }
  return geostr;
}

// save image info to DIF structure
char *CmoilBitmap::GetPixels () 
{
   const int PIXEL_SIZE=4;			  // 3 for BGR in bytes, 4 for BRGA, bit_depth=8bits
   int lastBytes = imageBytes;

   imageWidth  =  glutGet(GLUT_WINDOW_WIDTH);   // length in pixels
   imageHeight =  glutGet(GLUT_WINDOW_HEIGHT)-MenuBar::MenuButtonHeight-1;    // take out menu bar
   	
   // 1. even length for converting to mpeg 
   // 2. there's a requirement that the number of bytes in each data row of a Windows DIB be a multiple of 4
   imageWidth  -=  (imageWidth%2);       
   imageHeight -=  (imageHeight%2);

   imageBytes = PIXEL_SIZE * imageWidth * imageHeight ; // total image size in bytes 
   if ( lastBytes < imageBytes )
   {
     if ( pixels != NULL)
       delete[] pixels;
   }
   pixels = new char[imageBytes];
   if ( pixels != NULL ) 
   {
     if (PIXEL_SIZE == 4 )
     {
       pixelFormat="RGBA";
       ::glReadPixels(0, 0, imageWidth, imageHeight, GL_RGBA,  GL_UNSIGNED_BYTE, pixels);
     }
     else // if (PIXEL_SIZE==3)
     {
       pixelFormat="RGB";
       ::glReadPixels(0, 0, imageWidth, imageHeight, GL_RGB, GL_UNSIGNED_BYTE, pixels);
     }
   }
   return pixels; 
}
}
