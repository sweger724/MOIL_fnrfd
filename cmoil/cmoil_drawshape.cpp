//************************************************************************************
//* Filename:  cmoil_drawshape.cpp
//* 
//* Description:
//*     For diplaying structural information.
//*
//* Modification History:
//*  Date        	Developer	Description
//*  ------------	-----------	-------------------------------------------------------
//*  Nov 2 2005	B. Wang	Initial Development
//*    
//************************************************************************************
#include "cmoil.h"
#include "cmoil_globals.h"
#include "cmoil_drawshape.h"

namespace CMOIL {
 
const int    Stick::StickListStart=100000;
const double Stick::StickRadius=0.1;
const float Tube::TubeRadius=Stick::StickRadius*2;   // std tube radius
const float SplineShape::StdRibbonWidth=1.0; 

const DiskArray Tube::UnitCircle={   {TubeRadius,0.,0.},
                                           {TubeRadius*cos(PI10th),TubeRadius*sin(PI10th),0},                                           
                                           {TubeRadius*cos(PI10th*2),TubeRadius*sin(PI10th*2),0},
                                           {TubeRadius*cos(PI10th*3),TubeRadius*sin(PI10th*3),0},                                           
                                           {TubeRadius*cos(PI10th*4),TubeRadius*sin(PI10th*4),0},
                                           {TubeRadius*cos(PI10th*5),TubeRadius*sin(PI10th*5),0},
                                           {TubeRadius*cos(PI10th*6),TubeRadius*sin(PI10th*6),0},
                                           {TubeRadius*cos(PI10th*7),TubeRadius*sin(PI10th*7),0},
                                           {TubeRadius*cos(PI10th*8),TubeRadius*sin(PI10th*8),0},
                                           {TubeRadius*cos(PI10th*9),TubeRadius*sin(PI10th*9),0},
                                           {-TubeRadius,0.,0},
                                           {-TubeRadius*cos(PI10th),-TubeRadius*sin(PI10th),0},
                                           {-TubeRadius*cos(PI10th*2),-TubeRadius*sin(PI10th*2),0},
                                           {-TubeRadius*cos(PI10th*3),-TubeRadius*sin(PI10th*3),0},
                                           {-TubeRadius*cos(PI10th*4),-TubeRadius*sin(PI10th*4),0},
                                           {-TubeRadius*cos(PI10th*5),-TubeRadius*sin(PI10th*5),0},
                                           {-TubeRadius*cos(PI10th*6),-TubeRadius*sin(PI10th*6),0},
                                           {-TubeRadius*cos(PI10th*7),-TubeRadius*sin(PI10th*7),0},
                                           {-TubeRadius*cos(PI10th*8),-TubeRadius*sin(PI10th*8),0},
                                           {-TubeRadius*cos(PI10th*9),-TubeRadius*sin(PI10th*9),0} };   //for speed

//==================================================================================//

void Tube::AddToQuadricList(GLUquadricObj *qobj)
{
  if ( !glIsList(TubeEndListNum) )              // add sphere for ending tube
  {
    glNewList(TubeEndListNum, GL_COMPILE); 
    gluSphere(qobj, /*Radius*/ TubeRadius, /*Slices*/ 12,  /*Stacks*/ 12 ); 
    glEndList();
  }
}

void Tube::DrawEndBall(Coordinate &center) 
{
  glTranslatef(center.x, center.y, center.z);	   // Move to the proper location
  glCallList(TubeEndListNum); 			   // Render sphere display list
  glTranslatef(-center.x, -center.y, -center.z); // Move back
}

// rotate unit cirle from Z-axis to "direction" firstat then translate it to "center"
void Tube::RotTransDisk(Coordinate center, Coordinate direction, DiskArray &newDisk)
{
   glMatrixMode(GL_MODELVIEW);
   glPushMatrix();
   glLoadIdentity();

   double angle= acos((ZAxis^direction.norm()))*PI2DEGREE;  // rotate unit cirle to the position
   Coordinate normal=(ZAxis*direction).norm(); 
   glRotated(angle, normal.x, normal.y, normal.z); 
   GLfloat m[16];
   glGetFloatv(GL_MODELVIEW_MATRIX, m );   // get rotation matrix

   for (int i=0; i<TubeDiskSlice; i++)                // get coordinates
   {
      newDisk[i][0]=UnitCircle[i][0]*m[0]+UnitCircle[i][1]*m[4] +center.x;       // m[] is column major
      newDisk[i][1]=UnitCircle[i][0]*m[1]+UnitCircle[i][1]*m[5] +center.y; 
      newDisk[i][2]=UnitCircle[i][0]*m[2]+UnitCircle[i][1]*m[6] +center.z; 
   }
   glPopMatrix();
}

//preCntr, preDisk will be updated in the method
void Tube::DrawTubeBtw2Centers(Coordinate center, Coordinate &preCntr, DiskArray &preDisk, bool isFirst)
{
    DiskArray disk;

    Coordinate direction=preCntr-center;
    if (isFirst)      // get first disk as well
       Tube::RotTransDisk(preCntr,direction,preDisk);   //first disk takes current direction
    Tube::RotTransDisk(center,direction,disk);
    
    glBegin(GL_QUAD_STRIP);
    int minii=0;   // added for the problem with simple draw algorithm
    {
      float mind=100, d;
      for (int ii=0; ii<TubeDiskSlice; ii+=2)
      {
          d = (preDisk[0][0]-disk[ii][0])*(preDisk[0][0]-disk[ii][0])
            + (preDisk[0][1]-disk[ii][1])*(preDisk[0][1]-disk[ii][1])
            + (preDisk[0][2]-disk[ii][2])*(preDisk[0][2]-disk[ii][2]);

        if ( d<mind)
        {
          mind=d;
          minii=ii;
        } 
      }
    }

    int j;
    Coordinate v;
    for(int i=0; i<TubeDiskSlice; i++)
    {
      v.set(preDisk[i][0]-preCntr.x, preDisk[i][1]-preCntr.y, preDisk[i][2]-preCntr.z);
      v.norm();
      glNormal3d(v.x, v.y, v.z);
      glVertex3fv(preDisk[i]);
      j=(minii==0)?i:(i+minii)%TubeDiskSlice;
      v.set(disk[j][0]-center.x, disk[j][1]-center.y, disk[j][2]-center.z);
      v.norm();
      glNormal3d(v.x, v.y, v.z);
      glVertex3fv(disk[j]);
    }
    v.set(preDisk[0][0]-preCntr.x, preDisk[0][1]-preCntr.y, preDisk[0][2]-preCntr.z);
    v.norm();
    glNormal3d(v.x, v.y, v.z);
    glVertex3fv(preDisk[0]);
    v.set(disk[0][0]-center.x, disk[0][1]-center.y, disk[0][2]-center.z);
    v.norm();
    glNormal3d(v.x,v.y,v.z);
    glVertex3fv(disk[0]);
    glEnd();

    memcpy(preDisk, disk, sizeof(preDisk)); 
    preCntr = center;    
}


//
void Tube::Draw(Atom *myatom)
{
  bss[0].NextPoint(myatom);
  for(int i=0; i<bss[0].numOut; i++)
  {
    setColor(bss[0].out_atom[i]->rn, colorbase); 
    Draw(bss[0].a_out[i], bss[0].out_atom[i]->secStruct); 
  }
}

void Tube::Draw(Coordinate center, SecStructType secstruct) 
{    
  if ( count > 0 ) // count==1, no direction info yet
  { 
     DrawTubeBtw2Centers(center,preCntr,preDisk, (count==1)); 
  } 
  else   
  {
    DrawEndBall(center);   // count==0: draw starting ball
    preCntr=center;
  }
  count++;
}

//==================================================================================//

void Curve::Draw(Atom *myatom)
{
  bss[0].NextPoint(myatom);
  for(int i=0; i<bss[0].numOut; i++)
  {
    setColor(bss[0].out_atom[i]->rn, colorbase); 
    glVertex3d(bss[0].a_out[i].x, bss[0].a_out[i].y, bss[0].a_out[i].z); 
  }
}

//==================================================================================//
void Ribbon::Draw(Coordinate newCA, Coordinate newO, int rn, SecStructType secstruct)
{
    Atom *at;

    nvp = quad(newCA, newO, 0);
    if ( nvp != NULL )                 // output mid-point
    {
        currCaIndex = (currCaIndex+1)%4;
        at = input[0]+ currCaIndex;
        *at=nvp[2]; 
        at->rn = rn;                   // assign current rn to it
        bss[0].NextPoint(at);             // bspline for P1s

        at = input[1]+ currCaIndex;
        *at=nvp[3]; 
        at->rn = rn;                   // assign current rn to it
        bss[1].NextPoint(at);             // bspline for P2s
        DrawDo(); 
    }
}

void Ribbon::DrawDo()
{  
   int numOut=bss[0].numOut;

   if ( numOut < 1 )
      return;
   if ( numOut == 1 )
   {
     for ( int ii=0; ii<numBss; ii++)
       lastDraw[ii]=bss[ii].a_out[0];        // reset lastDraw1 upon first point of the chain
   }
 
   for (int k=0; k<numOut; k++ )
   {
     Coordinate nmv=normal(lastDraw[0],bss[1].a_out[k],bss[0].a_out[k]);

     if ( (lastNorm1^nmv) < 0 )
       nmv = (-1)*nmv;
     glNormal3d(nmv.x, nmv.y, nmv.z);
     lastNorm1=nmv;

     setColor(bss[0].out_atom[k]->rn, colorbase); 
     glVertex3d(bss[0].a_out[k].x, bss[0].a_out[k].y ,bss[0].a_out[k].z);
     glVertex3d(bss[1].a_out[k].x, bss[1].a_out[k].y ,bss[1].a_out[k].z); 
   }
   for ( int jj=0; jj<numBss; jj++)
   {
     lastDraw[jj]=bss[jj].a_out[numOut-1];
     bss[jj].numOut=0;     // finished the drawing from output, reset output array
   }
}

void Ribbon3D::DrawDo()
{  
   if (bsplineOutCount >= 1)
   {
     tRGBA clr = *colorbase; 
     Draw4Quadstrips(0,bsplineOutCount-1,&clr);  
     DrawQuadstripEnd(0, &clr);   // seal the ends
     DrawQuadstripEnd(bsplineOutCount-1, &clr);

     for (int j=0;j<numBss;j++)
       bss[j].reset();              // finished the drawing from output, reset output array
   }
}

//==================================================================================//
void SecStruct::ChainBegin()
{
  SplineShape::ChainBegin();
  oldStructLow=STRUCT_UNKNOWN;
}

void SecStruct::DrawBegin(SecStructType secstruct) 
{
   bsplineOutStruct=secstruct|STRUCT_INIT;
   Draw(); 
}

void SecStruct::DrawEnd()
{
    if ( bss[0].numOut > 1 )                                                                
      Draw();
    glEnd(); 
}

SecStruct::SecStruct(AtomArray *aaptr, int windowSize, int colorfilter) : SplineBelt2D(aaptr,windowSize,colorfilter)
{
   oldStructLow=STRUCT_UNKNOWN; 
   bsplineOutStruct=STRUCT_UNKNOWN;
   prevStruct=STRUCT_UNKNOWN;
   structFactor=0.0;      
   secStruct=aaptr->atomList[0].secStruct;
}

// use rn=-1 for newCA|newO not found
void SecStruct::Draw(Coordinate newCA, Coordinate newO, int newRn, SecStructType newStruct)
{
    Atom *at;

    if (newRn > -1 )
    {
      nvp = quad(newCA, newO, 0);        // unit ribbon spline quad points, between preCA,O <-> newCA

      if ( nvp != NULL )                       // output mid-point
      {
        currCaIndex = (currCaIndex+1)%4;

        at = input[0]+ currCaIndex;
        *at=nvp[2]; 
        at->rn = newRn;              // assign current rn to it
        bss[0].NextPoint(at);           // bspline for P1s

        at = input[1]+ currCaIndex;
        *at=nvp[3]; 
        at->rn = newRn;              // assign current rn to it
        bss[1].NextPoint(at);           // bspline for P2s

        Draw();
      }
    }
    bsplineOutStruct=prevStruct;    // advance struct Types
    prevStruct=newStruct;
    secStruct=newStruct;
}

void  SecStruct::Draw()
{
   Coordinate    S1, S2, deltaS, nmv;
   int    k, mid, nhelix=0;
   double direction=1.;
   float deltapoints;
   SecStructType secStructLow = sLowType(bsplineOutStruct);
 
   if (bsplineOutStruct&STRUCT_INIT)  
   {
     StructEndBegin(oldStructLow, secStructLow); 
     oldStructLow=secStructLow;
   }
   else if (secStructLow!=oldStructLow) 
   {
     StructEndBegin(oldStructLow, secStructLow); 
     oldStructLow=secStructLow;
   }

   // start to draw, and save points to lastDrawX
   if ( (secStructLow==STRUCT_LASTOFSHEET) && bss[0].numOut > 1 )   
   {                                                  // start an arrow, modify width
     deltaS = (structFactor-0.45)*(lastDraw[0]-lastDraw[1]);
     S1 = lastDraw[0] + deltaS;
     S2 = lastDraw[1] - deltaS;
     glVertex3d(S1.x, S1.y ,S1.z );
     glVertex3d(S2.x, S2.y ,S2.z );
     lastDraw[0]=S1;
     lastDraw[1]=S2;

     deltapoints = 1.1/(bss[0].numOut-1.);              // arrow edge delta
   }

   if (bss[0].numOut == 1 )
   {
     for ( int ii=0; ii<numBss; ii++)
       lastDraw[ii]=bss[ii].a_out[0];    // reset lastDraw1 upon first point of the chain
   }
   else if ( sIsHelix(secStructLow) /* && (secStructLow!=oldStructLow) */    )
   {  
      if (secStructLow!=oldStructLow)
        nhelix=0;
      else
        nhelix++;
      if ( nhelix < 2 )           // check helix's normal
      {    
        mid = bss[0].numOut/2;
        if ( ( normal(lastDraw[0],bss[1].a_out[mid],bss[0].a_out[mid])^
             (bss[0].a_out[mid]-avg(lastDraw[0],bss[0].a_out[bss[0].numOut-1]))   )  < 0 )
          direction=-1;
        else
          direction=1;
      }
   }

   for (k=0; k<bss[0].numOut; k++ )
   {
     nmv= direction * normal(lastDraw[0],bss[1].a_out[k],bss[0].a_out[k]);
     glNormal3d(nmv.x, nmv.y, nmv.z);
     
     if ( colorsetting != NULL)
     {
       setColor(bss[0].out_atom[k]->rn, colorbase);   
     }
     
     S1=bss[0].a_out[k];
     S2=bss[1].a_out[k];
     if ( (secStructLow==STRUCT_LASTOFSHEET) )
     {
       structFactor -= deltapoints;           // get arrow edge
       lastDraw[1]=S2;
     }
     deltaS = structFactor * ( S1 - S2 );
     S1 = S1 + deltaS;
     glVertex3d(S1.x, S1.y ,S1.z );
     S2 = S2 - deltaS;
     glVertex3d(S2.x, S2.y ,S2.z );           //overhead: if ( secStruct != STRUCT_UNKNOWN)
   }
   if ( bss[0].numOut > 0 )
   {
     lastDraw[0]=S1;
     lastDraw[1]=S2; 
     for ( int ii=0; ii<numBss; ii++) 
        bss[ii].numOut=0;              // finished output
   }
}  // Draw()

//
void SecStruct::StructEndBegin( SecStructType oldStruct, SecStructType newStruct )
{
  // set Color
  //
  tRGBA clr(255,0,0);       // red

  switch (newStruct)
  {
    case STRUCT_HELIX:  structFactor=0.0; 
                        clr.set(colorbase->r,colorbase->g,colorbase->b); 
                        break;
    case STRUCT_SHEET:  structFactor=0.3; 
                        clr.set(255,0,0);
                        break;
    case STRUCT_LASTOFSHEET:                   // draw same color arrow
                        structFactor=1.0;
                        break;
    default:            structFactor=-0.50;   //-0.45; 
                        clr.set(255,255,0);
   }
   if (newStruct != STRUCT_LASTOFSHEET ) 
      setColor(aaPtr->atomList->rn,&clr);    

    // switch between GL_QUAD_STRIP(helix|sheet)  and GL_LINE_STRIP (loop)
    //
    if ( sIsUnknown(newStruct)|| sIsUnknown(oldStruct))
    {
      glEnd();   //pairing 
      if ( sIsUnknown(newStruct)) 
      {
        glLineWidth(5.0);
        glBegin(GL_LINE_STRIP);
        glNormal3d(0,0,0);
        if (bss[0].numOut > 1 )
        {                                 // repeat last point as the beginning of new struct type
          lastDraw[0]=avg(lastDraw[0], lastDraw[1]);
          glVertex3d(lastDraw[0].x, lastDraw[0].y ,lastDraw[0].z );
        } 
      }
      else
      {
        glBegin(GL_QUAD_STRIP);
        if (bss[0].numOut > 1 )
        {
          glVertex3d(lastDraw[0].x, lastDraw[0].y ,lastDraw[0].z );
          glVertex3d(lastDraw[1].x, lastDraw[1].y ,lastDraw[1].z );
        }
      }
   }
}

//==================================================================================//
const float SplineBelt3D::thickness=Tube::TubeRadius*2;

// draw a face upon bsplineOut[0-3], edge1 & edge2 are indexes of bsplineOut[]
//
void SplineBelt3D::DrawQuadstrip(int startIndex, int stopIndex, int edge1, int edge2, bool SetDrawDelta, tRGBA *dflColor)
{
   Coordinate nmv, *out1, *out2, *out3, *out4;

   out1=bsplineOut[edge1]+startIndex;
   out2=bsplineOut[edge2]+startIndex;
   if (SetDrawDelta)        // merged SetOutputArrays() for speed
   {
     out3=bsplineOut[edge1+2]+startIndex;  // 0->2
     out4=bsplineOut[edge2+2]+startIndex;   // 1->2
   }
   glBegin(GL_QUAD_STRIP);
   for (int i=startIndex; i<=stopIndex; i++, out1++, out2++)
   {
     setColor(resNum[i], dflColor); 
     if ( i != stopIndex )
       nmv= normal(*out1 , *out2, *(out2+1) ); 
     glNormal3d(nmv.x, nmv.y, nmv.z);
     if (SetDrawDelta)                   // get the sheet which is parallel, delta away, to the edge1-edge2 face
     {
        *out3 = *out1 + (thickness*nmv);    // also set the coordinate
        glVertex3d(out3->x, out3->y, out3->z); 
        *out4 = *out2 + (thickness*nmv);
        glVertex3d(out4->x, out4->y, out4->z); 
        out3++;
        out4++;
     }
     else
     {
       glVertex3d(out1->x, out1->y, out1->z); 
       glVertex3d(out2->x, out2->y, out2->z); 
     }
   }
   glEnd();
}

// seal the ribbon Ends
void SplineBelt3D::DrawQuadstripEnd(int outIndex, tRGBA *dflColor)
{
   tRGBA clr = 0.7 * (*dflColor);
   setColor(resNum[outIndex], &clr); 

   glBegin(GL_QUADS);
   Coordinate *out1=bsplineOut[0]+outIndex; 
   Coordinate *out2=bsplineOut[1]+outIndex;
   Coordinate *out3=bsplineOut[2]+outIndex;
   Coordinate *out4=bsplineOut[3]+outIndex; 
   Coordinate   nmv=normal(*out1, *out2, *out3);
   glNormal3d(nmv.x, nmv.y, nmv.z);
   glVertex3d(out1->x, out1->y, out1->z);
   glVertex3d(out2->x, out2->y, out2->z);
   glVertex3d(out4->x, out4->y, out4->z);
   glVertex3d(out3->x, out3->y, out3->z);
   glEnd();
}

//dflColor changed to darker color after call this method
void SplineBelt3D::Draw4Quadstrips(int startIndex, int stopIndex, tRGBA *dflColor)
{
   glEnable(GL_LIGHTING);

   // draw two main faces
   DrawQuadstrip(startIndex,stopIndex,0,1,false,dflColor);  // bsplineOut[0,1], also set bsplineOut[2,3], must called first 
   DrawQuadstrip(startIndex,stopIndex,0,1,true, dflColor);

   //draw side faces with darker color
   tRGBA clr = 0.7 * (*dflColor);
   DrawQuadstrip(startIndex,stopIndex,0,2,false,&clr);
   DrawQuadstrip(startIndex,stopIndex,1,3,false,&clr);
}


// get points of two bsplines, and save them into bsplineOut[0,1]
// DoSpline();
void SplineBelt3D::Draw(Coordinate newCA, Coordinate newO, int rn, SecStructType newStruct)
{
    Atom *at;

    nvp = quad(newCA, newO, 0);
    if ( nvp != NULL )                 // output mid-point
    {
        currCaIndex = (currCaIndex+1)%4;

        at = input[0]+ currCaIndex;
        *at=nvp[2];                    // upper point
        at->rn = rn;                   // assign current rn to it
        at->secStruct=newStruct;
        bss[0].NextPoint(at);             // bspline for P1s
        int numOut=bss[0].numOut;

        memcpy(bsplineOut[0]+bsplineOutCount, bss[0].a_out, sizeof(Coordinate)*numOut); 
   
        for ( int j=0; j<numOut; j++) 
	  {
          resNum[bsplineOutCount+j] = bss[0].out_atom[j]->rn; 
          structTypes[bsplineOutCount+j] = bss[0].out_atom[j]->secStruct;
	  }

        at = input[1]+ currCaIndex;    // lower point
        *at=nvp[3]; 
        at->secStruct=newStruct;
        at->rn = rn;                   // assign current rn to it
        bss[1].NextPoint(at);             // bspline for P2s
        memcpy(bsplineOut[1]+bsplineOutCount, bss[1].a_out, sizeof(Coordinate)*numOut); 

        bsplineOutCount+=numOut;
    }

    bsplineOutStruct=prevStruct;    // advance struct Types
    prevStruct=newStruct;
    secStruct=newStruct;
}

SplineBelt3D::SplineBelt3D(AtomArray *aaptr, int windowSize, int colorfilter):SplineShape(aaptr,windowSize,colorfilter,2)
{
   for (int i=0; i<NumEdges; i++)
     bsplineOut[i]=new Coordinate[5*aaptr->numRes];
   resNum = new int[5*aaptr->numRes];
   structTypes = new char[5*aaptr->numRes];
   oldStructLow=STRUCT_UNKNOWN;
   bsplineOutCount=0; 
}

SplineBelt3D::~SplineBelt3D() 
{
  for (int i=0; i<NumEdges; i++)
    if ( bsplineOut[i] != NULL)
      delete[] (bsplineOut[i]);
  if (resNum != NULL)
      delete[] resNum;
  if (structTypes !=NULL)
      delete[] structTypes;
}

//==================================================================================//
void SecStruct3D::DrawHelix(int startIndex, int stopIndex)
{
   tRGBA clr= *colorbase;      // default to ribbon color
   Draw4Quadstrips(startIndex,stopIndex,&clr);  
   DrawQuadstripEnd(startIndex, &clr);   // seal the ends
   DrawQuadstripEnd(stopIndex, &clr);
}

void SecStruct3D::DrawSheet(int startIndex, int stopIndex)
{
   // modify output for sheet width
   int i;
   Coordinate last0=bsplineOut[0][stopIndex];  // last point are modified during drawing
   Coordinate last1=bsplineOut[1][stopIndex];  

   for (i=startIndex; i<=stopIndex-6; i++)
   {
      Coordinate delta =  0.05*(bsplineOut[0][i]-bsplineOut[1][i]);  // 0.6*2 of ribbon
      bsplineOut[0][i] = bsplineOut[0][i]+delta;
      bsplineOut[1][i] = bsplineOut[1][i]-delta;
   } 
   for (; i<=stopIndex; i++)   // arrow
   {
      Coordinate delta = ((0.45*(stopIndex-i)-1.)*0.55)*(bsplineOut[0][i]-bsplineOut[1][i]) ;
      bsplineOut[0][i] = bsplineOut[0][i]+delta;
      bsplineOut[1][i] = bsplineOut[1][i]-delta;  
   }

   tRGBA clr(255,0,0);                  // default to red
   Draw4Quadstrips(startIndex,stopIndex, &clr);
   DrawQuadstripEnd(startIndex, &clr);  // seal the beginning end

   bsplineOut[0][stopIndex]=last0;      // recover last point for computing next nmv 
   bsplineOut[1][stopIndex]=last1;
}

// refer to SplineBelt3D::DrawQuadstrip() for the thickness
void SecStruct3D::GetLoopCenters(int startIndex, int stopIndex )   // for tube centers
{
   Coordinate *out1=bsplineOut[0]+startIndex;
   Coordinate *out2=bsplineOut[1]+startIndex;
   Coordinate *out3=bsplineOut[2]+startIndex;
   for (int i=startIndex; i<=stopIndex; i++, out1++, out2++, out3++)
   {
     if ( i != stopIndex )
     {
       Coordinate nmv= normal(*out1 , *out2, *(out2+1) ); 
       *out3 = 0.5*(*out1 + *out2 + (thickness*nmv) );    // set back to bsplineOut[0]
     }
   }
}

// draw a tube upon bsplineOut[0,1], different than drawing tube from feeding points in 1D
void SecStruct3D::DrawLoop(int startIndex, int stopIndex)
{
  if (stopIndex >= bsplineOutCount) 
     stopIndex=bsplineOutCount-1;

  GetLoopCenters(startIndex, stopIndex);   // output in  bsplineOut[2]

  Coordinate preCntr=bsplineOut[2][startIndex];
  DiskArray  preDisk;    // initial in DrawTubeBtw2Centers()

  tRGBA clr(255,255,0);    // default to yellow
  if ( stopIndex>=startIndex)
  {
     setColor(resNum[startIndex], &clr); 
     Tube::DrawEndBall(preCntr);    //draw starting ball
  }

  for(int i=startIndex+1; i<stopIndex; i++)
  {
     setColor(resNum[i], &clr); 
     Tube::DrawTubeBtw2Centers( bsplineOut[2][i], preCntr, preDisk, (i==startIndex+1));  
  }

  if ( stopIndex>startIndex)
  {
    setColor(resNum[startIndex], &clr); 
    Tube::DrawEndBall(preCntr);    //draw starting ball
  }
}

void SecStruct3D::DrawDo() 
{
   AtomArray *aptr=bss[0].aaptr;
   int i, ii=0, jj=0; 
   while ( ii<bsplineOutCount)
   {
     bsplineOutStruct= (structTypes[ii]&0x3);
     for (i=ii; i<bsplineOutCount; i++)
     {
       if (bsplineOutStruct != (structTypes[i]&0x3))
         break;
     } 
     jj=(i==bsplineOutCount)?i-1:i;
     if ( sIsHelix(bsplineOutStruct) )
        DrawHelix(ii,jj);
     else if (sIsSheet(bsplineOutStruct) )  
        DrawSheet(ii,jj);
     else 
        DrawLoop(ii,jj+1);
     ii=i;
   }
   for (int j=0;j<numBss;j++)
     bss[j].reset(); 
}

void SecStruct3D::DrawBegin(SecStructType secstruct) 
{
   bsplineOutStruct=secstruct|STRUCT_INIT;
   Draw(); 
}

void SecStruct3D::DrawEnd()
{
   if ( bss[0].numOut > 1 )                                                                
     Draw();
   DrawDo();
}

SecStruct3D::SecStruct3D(AtomArray *aaptr, int windowSize, int colorfilter) : SplineBelt3D(aaptr,windowSize,colorfilter)
{
   oldStructLow=STRUCT_UNKNOWN; 
   bsplineOutStruct=STRUCT_UNKNOWN;
   prevStruct=STRUCT_UNKNOWN;
   structFactor=0.0;      
   secStruct=aaptr->atomList[0].secStruct;
}

//==================================================================================//
void Stick::Draw(Atom *base, Atom *top)   
{
   glEnable(GL_LIGHTING);
   Coordinate direction(top->x-base->x, top->y-base->y, top->z-base->z);
   double r = direction.len();
   int l = listNum(r);
   if ( !glIsList(l) )
      AddToQuadricList(base, top);
     
   double angle= acos((ZAxis^direction.norm()))*PI2DEGREE;  // rotate unit cirle to the position
   Coordinate normal=(ZAxis*direction).norm(); 

   glPushMatrix(); 
   glTranslated(base->x, base->y, base->z);
   glRotated(angle, normal.x, normal.y, normal.z); 
   glCallList(l);
   glPopMatrix();
}

void Stick::AddToQuadricList(Atom *ptr1,  Atom *ptr2)
{ 
  Coordinate crd(ptr1->x-ptr2->x, ptr1->y-ptr2->y, ptr1->z-ptr2->z);
  double r= crd.len();       // atom distance distance
  int    l=listNum(r); 
  if ( !glIsList(l) )
  {
    glNewList(l, GL_COMPILE);  // Create tube display list.
    gluCylinder(Main::qobj, StickRadius, StickRadius, r*0.5, 8, 1);  // half of r
    glEndList();
  }
}

//==================================================================================//
// if colorsetting=NULL, then use the defaultColor
void SplineShape::setColor(int rn, tRGBA *defaultColor)
{
   tRGBA *c=(colorsetting==NULL)?defaultColor:(colorsetting+rn);
   if (c != NULL ) 
      glColor3ub( Color2R(colorFilter,c),Color2G(colorFilter,c),Color2B(colorFilter,c) ); 
}

//  Get next group of vectors for drawing quads.
//  Every run, store CA, O coordinates as static variables, so do to 
//  preCA & preO (init as {0,0,0}s). 
//  Output:  {PeptidePlaneNormal, middle point,  deltaMiddle-, deltaMidell+} 
// 
//  A = ( CA - preCA ) ;
//  B = ( preO - preCA );
//  C = A x B                  :  define peptide plan
//  D = C x A                  :  in peptide plan and vertical to delta(CA)
//  P = ( preCA+CA)/2          :  middle point
//  P1 = P - D/|D|             :  up point
//  P2 = P + D/|D|             :  down point
//  P3 = P1+xx                 :  up front
//  P4 = P2+yy                 :  down front
//  return :   { C/|C|,  P,  P1, P2, P3, P4 }
//
//  reset=0:   not reset, continue bspline  (using)
//  reset=1:   take current input as first point
//  reset=-1:  ignore CA & O, reset all static variables  (using)
//  reset=-2:  do not take current input, just recaculated upon new ribbonWidth
//
Coordinate* SplineShape::quad(Coordinate CA, Coordinate O,  int reset)
{
  Coordinate A, B, C, D;
 
  if (reset!=0 && reset!=-2)
    numCA=0;

  if (reset==-1 )
    return NULL;
     
  if (reset != -2 )
  {
    numCA++;
 
    if ( numCA>1 )
    {   // output
      A = CA -  preCA;    	     // ( CA - preCA )
      B = preO - preCA ;	          // ( preO - preCA )
      C = A * B ;                  // C=AxB,  normal of the peptide plane

      v[0]=C.norm();

      D =  C * A ;                // D=CxA
      if ( ( preD ^ D )< 0 )      // peptide planes should flip less than 90 degree
        D= (-1) * D;
      D.norm();
      D = ribbonWidth * D;        //D=RibbonWidth * D/|D|
      preD=D;

      v[1] = avg(CA, preCA);      // P
      // v[1] = v[1] + ( ribbonRadiusFactor * v[0] );              // avg(preO, v[1]);
      v[1] = v[1] + v[0];         // avg(preO, v[1]);
      v[2] = v[1] - D;
      v[3] = v[1] + D;
    }
    else
      preD=CrdZero;

    // first CA, store CA, O info only
    preCA = CA;
    preO = O; 
  }
  else if ( numCA > 1) 
  {
    D=preD;
    D.norm();
    D = ribbonWidth*D;
    v[2] = v[1]-D;
    v[3] = v[1]+D;
  }

  if ( numCA>1)
  {
    return v;
  }
  else
    return NULL;
}  // quad()

// reset spline
void SplineShape::ChainBegin()   
{
   currCaIndex=-1;
   quad(CrdZero,CrdZero,-1); 
   for ( int i=0; i<numBss; i++)
   {
     bss[i].reset();
     //lastDraw[i].set(aaptr->atomList[0]);  
   }
}

void SplineShape::ChainEnd()
{
   for ( int i=0; i<numBss; i++ )
      bss[i].NextPoint(NULL);
}

SplineShape::SplineShape(AtomArray *aaptr, int windowsize, int colorfilter, int numbss) : numBss(numbss) 
{ 
   aaPtr = aaptr;

   currCaIndex=-1;
   quad(CrdZero,CrdZero,-1 );
   for ( int i=0; i<numBss; i++)
   {
     bss[i].init(aaptr);
     lastDraw[i].set(aaptr->atomList[0]);  
   }
   lastNorm1.set(0.,0.,1.);
   if (numBss>1)
     ribbonWidth=StdRibbonWidth*aaptr->InParm.ribbonWidthFactor;
   else
     ribbonWidth=0.0;

   colorFilter=colorfilter; 
   windowSize=windowsize;
   glLineWidth((float)(2.0*windowSize));
   colorbase = &(aaptr->InParm.ribbonColorBase);               // default ribbon color
   colorsetting = aaptr->ribbonColor;
   setColor(aaptr->atomList->rn, colorbase); 
           
   glEnable(GL_LIGHTING);
   glPolygonMode(GL_FRONT_AND_BACK,GL_FILL );
}

//==================================================================================//

void AtomArray::backbone(int windowSize, int colorFilter)
{
  if ( DrawBackboneAsCurve )
    backboneCurve(windowSize, colorFilter);
  else
    backboneTube (windowSize, colorFilter);
}

void AtomArray::backboneCurve(int windowSize, int colorFilter)
{
  Curve curve(this, windowSize, colorFilter);
  backbone(&curve, windowSize, colorFilter); 
}

void AtomArray::backboneTube(int windowSize, int colorFilter)
{
  Tube  tube(this, windowSize, colorFilter);
  backbone(&tube, windowSize, colorFilter); 
}

void AtomArray::backbone(SplineShape *Loop, int windowSize, int colorFilter)
{
  Atom *listi=atomList, *previ=atomList;
  
  Loop->ChainBegin(); 
  for (int I=0; I<numAtom; I++, listi++ ) 
  {
    if ( previ->rn != listi->rn && sIsChainEnd(previ->secStruct) ) // ending a chain
    {
      Loop->ChainEnd();
      Loop->ChainBegin();
    } 
    if( ! strcmp(listi->atomName, "CA") ) 
    {
      Loop->Draw(listi);   
    }

    if ( previ->rn != listi->rn )
       previ=listi;
  } 
  Loop->ChainEnd(); 
}  // backbone()   /* end of interpolation function */

//==================================================================================//

//
//  isRibbon=true: display ribbon only
//          =false: display secondary structure
//  pre:  AtomArray::setSecStruct() has been called, thus listi->secStruct is set
//          |-drawing--|
//          secStruct  oldRes               currRes    listi
//          |          |<-oldRes.secStruct->|          |
// -X--.----.--X--.----.--X--.--------------.--X--.----.--.--.----.
//  CA C    N  CA C    N  CA C              N  CA C    N  CA C
//  X1      X2<--Fit-->X3                   X4      <----- bspline
//
//
void AtomArray::quadStripStruct(SplineShape *splineShape, int windowSize, int colorFilter, bool isRibbon)
{
  Atom *listi=atomList,  *currRes=atomList;
  Coordinate currCA, currO;
  bool foundCA=false, foundO=false;

  splineShape->ChainBegin();
  splineShape->DrawBegin(listi->secStruct);
  for (int i=0; i<numAtom; i++, listi++)
  {
    if (listi->skip!=0 )
        continue;
  
    if ( currRes->rn != listi->rn)         // start a new residual, draw the ribbon/2ndStruct
    {
      if ( foundCA && foundO)
         splineShape->Draw(currCA, currO, currRes->rn, currRes->secStruct); // currCA and currO should be read in
      else
         splineShape->Draw(CrdZero, CrdZero, -1, currRes->secStruct);  // pass in secStruct
      if ( sIsChainEnd( currRes->secStruct))  
      {
        splineShape->ChainEnd();                    // draw next=last quad
        splineShape->DrawEnd();
        splineShape->ChainBegin(); // (listi->secStruct);
        splineShape->DrawBegin(listi->secStruct);
      }
      currRes=listi; 
      foundCA=false;
      foundO=false;
    } // if (currRes->rn)

    if( !strcmp(listi->atomName, "CA"))           // consider CAE, CAx
    {
      currCA = *listi;
      foundCA=true;
    } 
    else if ( !strcmp(listi->atomName, "O" ))        
    {
      currO = *listi;
      foundO=true;
    }
  } // for(i)
  if ( foundCA && foundO)
    splineShape->Draw(currCA, currO, currRes->rn, currRes->secStruct); // currCA and currO should be read in
  else
    splineShape->Draw(CrdZero, CrdZero, -1, currRes->secStruct);  // pass in secStruct
  splineShape->ChainEnd(); 
  splineShape->DrawEnd();
}  // quadStripStruct()

void AtomArray::ribbon(int windowSize, int colorFilter)
{
  Main::lightposition(); 
  if ( DrawBackboneAsCurve )
    ribbon2D(windowSize, colorFilter);
  else
    ribbon3D(windowSize, colorFilter);
}

void AtomArray::ribbon2D(int windowSize, int colorFilter)
{
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1.0);

  Ribbon ribbon(this,windowSize,colorFilter);
  quadStripStruct(&ribbon, windowSize, colorFilter, true); 

  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,0.0);

}  // ribbon()

void AtomArray::ribbon3D(int windowSize, int colorFilter)
{
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1.0);
  Ribbon3D ribbon(this,windowSize,colorFilter);
  quadStripStruct(&ribbon, windowSize, colorFilter, true); 
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,0.0);

}  // ribbon()

void AtomArray::struct2nd(int windowSize, int colorFilter)
{
  Main::lightposition(); 
  if ( DrawBackboneAsCurve )
    struct2nd2D(windowSize, colorFilter);
  else
    struct2nd3D(windowSize, colorFilter);

}  // struct2nd()

void AtomArray::struct2nd2D (int windowSize, int colorFilter)
{
  SecStruct secstruct(this,windowSize,colorFilter);
  quadStripStruct(&secstruct, windowSize, colorFilter, false); 
}

void AtomArray::struct2nd3D (int windowSize, int colorFilter)
{
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,1.0);
  SecStruct3D secstruct(this,windowSize,colorFilter);
  quadStripStruct(&secstruct, windowSize, colorFilter, false); 
  glLightModelf(GL_LIGHT_MODEL_TWO_SIDE,0.0);

}

}
