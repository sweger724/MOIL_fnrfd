//*************************************************************************************************
//*  Filename:   cmoil_movieprint.cpp
//*
//*  Description: 
//*    Save structure movie images for Windows
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  June 2001	Baohua Wang	Initial Release
//*  Spet 2007  Thomas Blom Update for Linux/Apple
//*************************************************************************************************

#ifdef __APPLE__
	#include "sys/time.h"	// for ru_utime
	#include "sys/wait.h"	// for id_t
	#include "i386/types.h"  // for u_int32_t
	#include "sys/stat.h"  // for dev_t
#endif

#include <stdio.h>
#include <string.h>
#include "cmoil_movieprint.h"
#include "gifmerge.h"
#include "cmoil_globals.h"
#include "cmoil_msgboard.h"

namespace CMOIL {
MoviePrint::MoviePrint(InputParameters &iParms)
{
   printing = 0;
   totalFrames = 0;
   framesPrinted = 0;
   delay = 0;
   rmGifDir=1;           // remove all intermediate gifs after merge
   init_path(iParms);
   geometry="0x0";
   fdist=1;
   oneFrameOnly=0;
}

// set up directory to save output
//
void MoviePrint::init_path(InputParameters &iParms)
{
   char *pChr;

   exedir = iParms.exedir;

   totalFrames = iParms.maxframes;
   gifloops = iParms.gifloops;
   strcpy (imageType, iParms.imagetype);

   delay = (int) (iParms.sleepseconds / 10 );
   framedist = iParms.framedistance;
   strcpy(imagePathBase, iParms.bmpName); 		// derive path from bmp filename
   if ( (pChr = strrchr(imagePathBase, '.') ) != NULL )   
     *pChr = '\0';
   strcat(imagePathBase, "_img");
   //_mkdir(imagePathBase); 
}

void MoviePrint::start_print(int PrintOneFrameOnly) 
{
  if ( rmGifDir )      // the img directory does not exist
    MKDIR(imagePathBase); 
  printing = 1;
  framesPrinted = 0;
  fdist=1;
  oneFrameOnly=PrintOneFrameOnly;
  if ( !strcmp(imageType, "gif") )  {
    char mrgGif[FILENAME_MAXLEN];
    sprintf(mrgGif, "%s.gif", imagePathBase);
    GifMerge::GIF_Merge_Start(mrgGif, delay, gifloops, 2);
  }
  print();
}

// movie name format:  "img%05d"
void MoviePrint::print(void) 
{
  if ( printing && (framesPrinted < totalFrames) && (framedist<=1 || framedist==fdist || oneFrameOnly != 0 ) ) 
  {
    CmoilBitmap  monoPicture;
    char *suffix;
   
    fdist=1;
    framesPrinted++;
    if ( strncmp(imageType, "mpg", 3)==0 ) {
      suffix="yuv";     // save to .yuv for mpg
    } else {
      suffix=imageType;
    }
    sprintf(imageTmpName, "%s%cimg%05d.%s", imagePathBase, DIRDELIM, framesPrinted, suffix); 
    geometry=monoPicture.PrintImg(imageTmpName); 	// Print picture as sequential file first
    
    if ( !strcmp(imageType, "gif") )  
    {
      if ( GifMerge::GIF_Merge(imageTmpName) != 0 )
        rmGifDir=0;
    }
    if ( framesPrinted >= totalFrames || oneFrameOnly != 0 ) 					
      end_print();
  } 
  else if (printing)
  {
    fdist++;
  }
}

#include "cmoil_msgboard.h"

//exedir     directory for cmoil as well as mpeg2enc parameter frame file
void MoviePrint::end_print() 
{
   int i;
   printing = 0; 		     	// call convert to merge gif files

   if ( !strcmp(imageType, "gif") ) 
   {   
     GifMerge::GIF_Merge_End();

      //clean up
     if ( rmGifDir != 0 )
     {
       for (i=1; i<= totalFrames; i++ )
       {
         sprintf(imageTmpName, "%s%cimg%05d.gif", imagePathBase, DIRDELIM, i);
         remove(imageTmpName);
       }
       RMDIR(imagePathBase);
       sprintf(imageTmpName, "%s.gif", imagePathBase); 
     }
     else
       rmGifDir=1;
   }
   else if ( !strcmp(imageType, "mpg") )
   {
     char cmd[256];
     char *parfile=GetMpeg2encPar();
    
     if ( parfile != NULL )
     {
       sprintf(cmd, "%smpeg2enc %s %s.mpg", exedir, parfile, imagePathBase);
       system(cmd);

       for (i=1; i<= totalFrames; i++ )
       {
         sprintf(imageTmpName, "%s%cimg%05d.yuv", imagePathBase, DIRDELIM, i);
         remove(imageTmpName);
       }
       RMDIR(imagePathBase);
       UNLINK(parfile);
       sprintf(imageTmpName, "%s.log", imagePathBase);
       remove(imageTmpName);              // remove statistics file

       sprintf(imageTmpName, "%s.mpg", imagePathBase);
     } 
   }
   else 
     strcpy(imageTmpName, imagePathBase);

   char msg[256];
   sprintf(msg, "Saved to: %s", imageTmpName);
   MsgBrd.set(msg);       // change later
}

// input the directory of the mpeg2enc paramenter frame file
char *MoviePrint::GetMpeg2encPar() {
  static char newpar[256]="";
  char   modelpar[256];
  FILE  *inpar=NULL, *outpar=NULL;
  int  w, h;

  sscanf(geometry, "%dx%d", &w, &h);
  sprintf(modelpar, "%s%cmpeg2enc.par", exedir, DIRDELIM); 
  sprintf(newpar,   "%s.par", imagePathBase);

  inpar=fopen(modelpar, "r");
  if (inpar == NULL)
    return NULL;
  outpar=fopen(newpar, "w");
  if (outpar== NULL)
  {
    fclose(inpar);
    return NULL;
  }

  // user modelpar[] as line buffer
  while ( (fgets(modelpar, 256, inpar)) != NULL ) {
    if ( !strncmp(modelpar, "$SOURCE$", 8) ) {
       fprintf(outpar, "%s%cimg%%05d\n", imagePathBase, DIRDELIM);
    } else if ( !strncmp(modelpar, "$STATISTICS$", 12) ) {
       fprintf(outpar, "%s.log\n", imagePathBase);
    } else if ( !strncmp(modelpar, "$NFRAMES$", 9) ) {
        fprintf(outpar, "%d\n", totalFrames);
    } else if ( !strncmp(modelpar, "$FIRSTFRAME$", 12) ) {
       fprintf(outpar, "%d\n", 1);            // may change to any starting frame later
    } else if ( !strncmp(modelpar, "$WIDTH$", 7) ) {
       fprintf(outpar, "%d\n", w);
    } else if ( !strncmp(modelpar, "$HEIGHT$", 8) ) {
       fprintf(outpar, "%d\n", h);
    } else {
      fprintf(outpar, "%s", modelpar);
    }
  }
  fclose(inpar);
  fclose(outpar);

  return newpar;
}

}
