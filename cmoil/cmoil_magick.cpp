//*************************************************************************************************
//*  Filename:   cmoil_magick.cpp
//*
//*  Description: 
//*    Main loop of graphic display in cmoil program.
//*
//*  stdin stream:
//*    256 chars for filename,  followed by the image file data
//*
//*************************************************************************************************
#include <stdio.h>
#include <fcntl.h>
#include "Magick++.h"

#ifdef _WIN32
#include <io.h>
#endif

const unsigned int FILENAME_MAXLEN  =	256;		// max length for a filename

int main() 
{
  int imageBytes, imageWidth, imageHeight, pixelSize;
  char *pixels, *pixelFormat, imgfilename[FILENAME_MAXLEN];
  FILE *infp=stdin;

#ifdef _WIN32
  if ( _setmode(0, _O_BINARY) == -1 )
    return -1;
#endif

  fread(imgfilename, FILENAME_MAXLEN, 1, infp);
  fread(&imageBytes, 4, 1, infp);
  fread(&imageWidth, 4, 1, infp);
  fread(&imageHeight, 4, 1, infp);
  fread(&pixelSize, 4, 1, infp);

  if (pixelSize==3)
    pixelFormat="RGB";
  else 
    pixelFormat="RGBA";    // == 4

  if ( (pixels = (char*)malloc(imageBytes+512)) != NULL )
  {
    fread(pixels, 512, imageBytes/512+1,infp);
    Magick::InitializeMagick(NULL);
    Magick::Blob  blob(pixels, imageBytes);
    Magick::Image image(blob, Magick::Geometry(imageWidth,imageHeight), 8, pixelFormat) ;
    image.write(imgfilename);   // write to new file format
  }
  return 0;
}
 


