#ifndef _CMOIL_DIB_H
#define _CMOIL_DIB_H
//****************************************************************************************************
//*  Filename:   cmoil_dib.hpp
//*
//*  Description: 
//*    CmoilBitmap class.  
//*
//*  History:
//*  Date		Developer	Description
//*  ------------ ----------- ------------------------------------------------------------------------
//*  June 2001	B. Wang	Initial Development
//*
//****************************************************************************************************
#include "cmoil_globals.h"
#include <stdio.h>

namespace CMOIL {

//for saving image to a bitmap
class CmoilBitmap : public State
{
private :
   char *pixels;	                  // pixels data array in RGB format
   char *pixelFormat;               // RGBA or RGB
   int imageWidth;
   int imageHeight;
   int imageBytes;			// number of pixels

   char *GetPixels ( );		// transfer image to DIB structure
public :
   CmoilBitmap ();

   ~CmoilBitmap ();
   char *PrintImg (const char *aFilename);	// Transform Image Pixels to new file format and save
                                                // return image geometry
} ;

}

#endif
