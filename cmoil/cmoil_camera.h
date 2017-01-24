#ifndef _CMOIL_CAMERA_H
#define _CMOIL_CAMERA_H

//*************************************************************************************************
//*  Filename:   cmoil_camera.h
//*
//*  Description: head file for Camera class
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. 02 2005	Baohua Wang	Reshape
//*
//*************************************************************************************************
#include "cmoil_crd.h"

namespace CMOIL {
class Camera
{
private:
   double ratio;
   double nearz;
   double farz;
   double wd2;
   double ndfl;

   double left;
   double right;
   double top;
   double bottom;
    
   double maxz;
   double minz;

   Coordinate r;

   Coordinate vp;             /* View position           */
   Coordinate vd;             /* View direction vector   */
   Coordinate vu;             /* View up direction       */
   Coordinate pr;             /* Point to rotate about   */
   double focallength;  	/* Focal Length along vd   */
   double radians;     	
   double eyesep;       	/* Eye separation  /2        */
   //int screenwidth;
   //int screenheight;
   
   void SetVp(double x=0.0, double y=0.0, double z=0.0);

public : 
   Camera();

   void SetSizes(int wWidth, int wHeight, double maxZ, double minZ);
   void LookAt(int fromLeftRigth, int stereoMode);
   void Clip(int front, int back);

   enum {LEFT=-1, MIDDLE=0, RIGHT=1};
   enum {ScreenSpace=0, ObjectSpace=1};
   int  eye;                  /* left or right */

   int  clipFront;     // total clip offset
   int  clipBack;      // total clip offset
   
   int  clipSpace;     // clip in object or window space, not used
   
};

#define Set_gluPerspective()  gluPerspective( 40.0,1.0,1.0,10000)

}

#endif
