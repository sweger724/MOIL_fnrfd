//************************************************************************************
//* Filename:  cmoil_bspline.cpp
//* 
//* Description:
//*    This file includes Bspline class which generates bspline fit to input points.
//*    It also includes methods which are in class AtomArray but are related to use
//*    bspline to generate smooth backbone and ribbon display.
//*
//*    The class Bspline is derived from Richard Gillilan's "bspline.c" code for XTM
//*    (ftp://bridgedec.chess.cornell.edu). Now it does not need all known points in 
//*    input array. It takes dynamic inputs and outputs N interpolated values for each 
//*    segment defined by current input point. 
//*
//*    The ribbon algorithm is from "Algorithm for ribbon models of proteins",  by
//*    M. Carson & C. Bugg,  J. Mol. Graphics Vol 4(2) June 1986 pp 121-122.
//* 
//* Note:
//*    Maximum number of interpolated points is set to 10 for cmoil.
//*    Do not use OpenGL evaluators for speed and flexibility.
//*
//* Modification History:
//*  Date        	Developer	Description
//*  ------------	-----------	-------------------------------------------------------
//*  May 30 2003	B. Wang	Initial Development
//*    
//************************************************************************************
#include "cmoil.h"
#include "cmoil_globals.h"
#include "cmoil_bspline.h"

namespace CMOIL {
 
Bspline::Bspline(AtomArray *inArray, int N)
{
  if (inArray == NULL)
    init(inArray, 0);
  else
    init(inArray, N);    
}

void Bspline::init(AtomArray *inArray, int inN)
{
  aaptr=inArray;
  N=inN;                   // number of interpolated points per segment
  numIn=0;
  numOut=0;
  currIn=-1;

  if (N>0 && N<=MaxNumSegments)      // max(N)=10
  {
   
    S[0][0] = 6.0/(N*N*N); S[0][1] = S[0][2] = S[0][3] = 0.0;
    S[1][0] = S[0][0]; S[1][1] = 2.0/(N*N); S[1][2] = S[1][3] = 0.0;
    S[2][0] = 1.0/(N*N*N); S[2][1] = 1.0/(N*N); S[2][2] = 1.0/N; S[2][3] = 0.0;
    S[3][0] = S[3][1] = S[3][2] = 0.0;  S[3][3] = 1.0;

    B[0][0] = -1; B[0][1] =  3; B[0][2] = -3; B[0][3] = 1;
    B[1][0] =  3; B[1][1] = -6; B[1][2] = 3; B[1][3] = 0;
    B[2][0] = -3; B[2][1] =  0; B[2][2] = 3; B[2][3] = 0;
    B[3][0] =  1; B[3][1] =  4; B[3][2] = 1; B[3][3] = 0;

    for (int i=0; i<NP; ++i) {  
      for (int j=0; j<NP; ++j) { 
        B[i][j] /= 6.0; 
      }
    } 
  } 
}

// if aptr==NULL, it's an end of a chain

void Bspline::NextPoint(Atom *aptr)
{
  if (aptr!=NULL)
  {
     numIn++;
     currIn = (currIn+1)%NP;
     a_in[currIn]=aptr;
     StepFit(false);
  }
  else                        
  {
    if (numIn > 2 )   // push last point out
    {
       psudoAtom.x=2*a_in[currIn]->x-a_in[(currIn+3)%4]->x;
       psudoAtom.y=2*a_in[currIn]->y-a_in[(currIn+3)%4]->y;
       psudoAtom.z=2*a_in[currIn]->z-a_in[(currIn+3)%4]->z;
       psudoAtom.rn=currIn;
       psudoAtom.secStruct=a_in[currIn]->secStruct;
       numIn++;
       currIn = (currIn+1)%4;
       a_in[currIn]=&psudoAtom;
    }
    StepFit(true);
    numIn=0;
  }
}

// This subroutine performs cubic-spline interpolation
// on the input array (3D points) and fills the output
// array with the interpolated values.  
//
// lastPoint=true:  the input is the end of a chain
//
void Bspline::StepFit(bool lastPoint) 
{
  float x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4;
  int tt=currIn;
  int myreturn=0;
  Atom *atom1=NULL, *atom2=NULL, *at;

  if (lastPoint==false || numIn>2)
  {
    if (numIn==3) 
    {             // first output point
        at = a_in[tt];

        x4 = at->x; y4 = at->y; z4 = at->z;
        tt = (tt+3)%NP;         // tt-1
        at = a_in[tt];

        x3 = at->x; y3 = at->y; z3 = at->z;
        atom2 = at;
        tt = (tt+3)%NP;         // tt-1
        at = a_in[tt];

        x2 = at->x; y2 = at->y; z2 = at->z;
        atom1 = at;
        x1 = 2*x2-x3;       y1 = 2*y2-y3;         z1 = 2*z2-z3;
    }
    else if ( numIn > 3 )
    {
        at = a_in[tt];

        x4 = at->x; y4 = at->y; z4 = at->z;
        tt = (tt+3)%NP;         // tt-1
        at = a_in[tt];

        x3 = at->x; y3 = at->y; z3 = at->z;
        atom2 = at;
        tt = (tt+3)%NP;         // tt-1
        at = a_in[tt];
        x2 = at->x; y2 = at->y; z2 = at->z;

        atom1 = at;
        tt = (tt+3)%NP;         // tt-1
        at = a_in[tt];

        x1 = at->x; y1 = at->y; z1 = at->z;
    } 
    else if (numIn==1)
    {                           // output first point
      a_out[0] = *(a_in[currIn]);   // currIn=tt
      out_atom[0] = a_in[currIn];
      numOut=1;        
      myreturn=1;
    }
    else
    {
      numOut=0;
      myreturn=1;
    }
  }
  else  // 1 or 2 points, no spline fit, numIn never be 3 here
  {
    numOut=numIn;
    tt= (tt+5-numOut)%NP;
    for ( int t=0; t<numOut; t++, tt++)
    { 
      tt %= NP;
      a_out[t] = *(a_in[tt]);  
      out_atom[t] = a_in[tt];
    }
    myreturn=1;           // flag to return
    numIn=0;    //reset
    currIn=-1;
  }

  if ( myreturn == 0 ) 
  {
    int i,j,k;
    numOut=N;

    G[0][0] = x1; G[0][1] = y1; G[0][2] = z1; G[0][3] = 1.0;
    G[1][0] = x2; G[1][1] = y2; G[1][2] = z2; G[1][3] = 1.0;
    G[2][0] = x3; G[2][1] = y3; G[2][2] = z3; G[2][3] = 1.0;
    G[3][0] = x4; G[3][1] = y4; G[3][2] = z4; G[3][3] = 1.0;

    /* matrix multiply to construct M */

    for (i = 0; i < NP; ++i) {
      for (j = 0; j < NP; ++j) {
        BG[i][j] = 0.0;
        for (k = 0; k < NP; ++k) {
          BG[i][j] += B[i][k]*G[k][j];
        }
      }
    }

    for (i = 0; i < NP; ++i) {
      for (j = 0; j < NP; ++j) {
        M[i][j] = 0.0;
        for (k = 0; k < NP; ++k) {
          M[i][j] += S[i][k]*BG[k][j];
        }
      }
    }

    /* forward difference loop given in Carson's paper */

    for (k=0; k<N; k++) {  /* inner loop over points in interval */
      for (i=NP-1; i>0; --i) {
        for (j=0; j<NP; ++j) {
          M[i][j] += M[i-1][j]; 
	  }
      }
      float m=M[NP-1][NP-1]; 
      if (m< 1.0e-5) {
        fprintf(stderr,"M = %f\n",m); 
        numOut=k;
        break;
      }
      a_out[k].set(M[3][0]/m,M[3][1]/m,M[3][2]/m);
      if ( k <= N/2 )
        out_atom[k]=atom1;
      else
        out_atom[k]=atom2;
    }
  }
}  /* end of interpolation function SepFit() */

}


