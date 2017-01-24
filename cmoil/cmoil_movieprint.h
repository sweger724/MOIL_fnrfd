#ifndef _CMOIL_MOVIEPRINT_H
#define _CMOIL_MOVIEPRINT_H
//****************************************************************************************************
//*  Filename:   movieprint.hpp
//*
//*  Description: Class for printing Movie on Windows
//*
//*
//*  History:
//*  Date		Developer	Description
//*  ------------ ----------- ------------------------------------------------------------------------
//*  June 2001	Baohua Wang	Initial Development
//*
//****************************************************************************************************
#include "cmoil_const.h"
#include "cmoil_dib.h"
#include "cmoil_parm.h"

namespace CMOIL {

class MoviePrint : public State
{
private:
  int totalFrames;				// total number of frames to be printed
  int delay;					// frame delay in 1/100th second
  int gifloops;					// number of loops for gif animation
  char imageType[4];				// gif, bmp, jpg etc.
  char imagePathBase[FILENAME_MAXLEN];	// image file directory for group of images 
  char imageTmpName[FILENAME_MAXLEN];	// temp storage for file name
  char cmd[FILENAME_MAXLEN*3];  		// ImageMagick convert command
  int  rmGifDir;                          // remove all intermediate GIFs after merge
  int  framedist;                         // required frame distance
  int  fdist;					// current frame distance
  int  oneFrameOnly;                      // print currFrameOnly
  void init_path(InputParameters &iParms);	// used for constructor

  char *geometry;                        // pointer to a string of WIDTHxHEIGHT 
  char *exedir;                          // executable directory
  char *GetMpeg2encPar();                // compose parameter file for mpeg2enc.exe

public:
  MoviePrint(InputParameters &iParms);

  int printing;			// status : 0/1
  int framesPrinted;
  void start_print(int PrintOneFrameOnly=0);
  void print();
  void end_print();       
};

}

#endif


