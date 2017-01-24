#ifndef _CMOIL_RMS_H
#define _CMOIL_RMS_H
//*************************************************************************************************
//*  Filename:   cmoil_rms.hpp
//*
//*  Description: head file for cmoil_rms.cpp
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Jul. 25 2000	Baohua Wang	Extract from nloopp code
//*
//*************************************************************************************************
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include "cmoil_atom.h"
#include "cmoil_globals.h"
 
namespace CMOIL {
const int MAX_SEQ=3000;           /* Max number of residues in a protein */
const int YES=1;
const int NO=0;

#define NO_MISSING_COORD(prot,i) \
 ((prot->coord[set]->x[i]<COORD_INF) && (prot->coord[set]->y[i]<COORD_INF) && (prot->coord[set]->z[i]<COORD_INF))

const double RMS_ZERO=0.0001;
const int    ALIGN_MAXLEN=200;

/* -------------------------------------------------------------
   COORDINATES
   ===========
   -------------------------------------------------------------
*/
typedef struct
{
  double  x[MAX_SEQ];
  double  y[MAX_SEQ];
  double  z[MAX_SEQ];
} COORDINATES;


// add CE for auto-alignment
class CeAlign  : public State
{
public:
  int   numAlignedMonos;
  int   numGaps;
  int   numAlignedAtoms;
  float rmsd;          // calling overlap only 

  CeAlign();
  void  matches(AtomArray &aa1, AtomArray &aa2, COORDINATES &match1, COORDINATES &match2, char alignType='B');
  int   alignAA2Match(AtomArray &aa, COORDINATES &match, int seq_1or2 );
  int   overlap(AtomArray &aa1, AtomArray &aa2);

private:
  char  alignType;                                 // 'B': backbone; 'C': Calpha
  int   align[ALIGN_MAXLEN];                       // maximum alignment segments=50

  int   align_2pdbs(char *pdb1, char chain1, char *pdb2, char chain2, char *alignfile );
  void  get_tmppdb(AtomArray &aa, char *tmppdb, char *chain);
  void  outfilename(AtomArray &aa1, AtomArray &aa2, char *outfile);

};

class OverlapAlign
{
public :
  OverlapAlign();
  int     overlap(AtomArray &aa1, AtomArray &aa2);
  char    ErrMsg[100];

  char    structs[2][10];     // store strunct:chain info
  int     numAlignedAtoms;
  float   rmsd;  
  CeAlign autoAlign;
  char    ovtype;           // align by: 'A'_seleted atoms, 'B'_backbone, 'C'_CA
  
private:
  void   OverlapMoveStruct(AtomArray &aa);
  int    cpyAA2Match(AtomArray &aa, COORDINATES &match, bool isOrg=true);
  void   setErrMsg();

  void   eigsrt(double d[], double **v, int n);
  int    jacobi(double **a, int n, double d[], double **v, int *nrot);
  double *vector(long nl, long nh);
  void   free_vector(double *v, long nl, long nh);
  double **convert_matrix(double *a, long nrl,long nrh, long ncl, long nch);
  void   free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
  int    diag_kabsch(double kab[3][3], double *eigval);
  int    rms_dist(COORDINATES *match1, COORDINATES *match2, int n);

  double  OverlapMatrix[5][3];    //  translate1[3], translate2[3], rotate[3][3] in order

};

}

#endif


