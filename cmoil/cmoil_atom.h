#ifndef _CMOIL_ATOM_H
#define _CMOIL_ATOM_H
//*************************************************************************************************
//*  Filename:   cmoil_atom.h
//*
//*  Description: head file for AtomArray class
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Oct. 09 2000	Baohua Wang	Extract from cmoil.cpp and add const definitions
//*  Aug. 02 2005	Baohua Wang	reshape
//*
//*************************************************************************************************
#include <stdio.h>
#include "cmoil_const.h"
#include "cmoil_parm.h"
#include "cmoil_crd.h"
#include "cmoil_movieprint.h"
#include "cmoil_globals.h"

namespace CMOIL { 

class Atom
{
public:
  Atom(void);
  double x;       /* coordinates for display */
  double y;
  double z;
  double orgX;     /* coordinates read from file */
  double orgY;
  double orgZ;
  double radius;
  int rn;			// residue number
  int rnPDB;		// residue number to display or the real PDB resNum, for supporting PDB file
  int chainIndex;		// chain index of the chain array
  NameString atomName;  // discard atom if atomName=""
  int neighbors[ATOM_MAX_NEIGHBORS];
  short int numNeighbors;
  tRGBA color;
  float charge;
  bool skip;
  bool skipsurf;        // for generating surface
  bool skipsurfdisp;    // for displaying surface
  bool align;		// selected atoms for overlap
  int  dispmode;
  SecStructType secStruct;    // type of 2nd structure

  const Atom &operator=(const Coordinate &a);
  friend double angleT(Atom &a0, Atom &a1, Atom &a2, Atom &a3) ;
};


// for display correct serial number from PDB file, map pdb num to unique internal index
//
class Chain
{
public:
  char   id;			// 22th char of a PDB Atom record
  int    model;
  int    resPdbStart;		// the start residue number in a PDB file
  int    resMapStart;		// the start residue index in the resName array
  int    atomPdbStart;        // the start atom numver in a PDB file
  int    atomMapStart;        // the start atom index in Atom array
};

class AtomArray : public State
{
  friend bool linePlaneCrossPoint(float *plane, Coordinate &pt1, Coordinate &pt2, float *crossPt);

private:
  // format of binary surface vertex file
  typedef struct _SurfHeader
  {
    int numSurfTris;		      // number of triangles in surface
    short int outerProced;		// outer surface processed
    short int cavityNormalReversed;   // reverse normal vector for cavity atoms
  } SurfHeader;         		// binary vertex file format: maxNumSurfTris + [header+SurfTriangle*numSurfTris]*M

  typedef struct _SurfTriangle      // line structure of *surface triangle
  {
    int   atomNum;
    float nm[3][3];                // normals
    float vt[3][3];                // vertexes
    int   color;                  
    int   nextAtom;                // neighbor atom number
    int   edge;
    int   out;                     // set to 1 if it's on outer surface, not cavity    
  } SurfTriangle;

  typedef struct _SurfEdgeIndex      // line structure of *surfIndex
  {
    int        startTri;             // first tri index in  *surface
    short int  triCnt;               // number of tris in this edge
    bool       out;
  } SurfEdgeIndex;

  enum { SURF_TRI_LEN=sizeof(SurfTriangle),     
         SURF_HDR_LEN=sizeof(SurfHeader) } ;

  enum { MaxNumSurfEdgePerAtom=25 };  // =MAX_DENSITY in surf.exe

  static double BOND_MAXLEN; //1.9 
  static double BOND_MINLEN; //0.85

  double versionNumber;		// version number of connectivity file

  FILE *crdfp;			// ReadPth() , Dcd, Crd
  FILE *vertexfp;             // surface vertex 
  int  *dcdNFrzAtoms;         // for read in dcd file
  int   nNoFrzAtoms;

  int   maxNumSurfTris;       // for loading path's single surf at a time
                              // if it's 0, the "surf" not run yet (flag)
  float Coeff[4];             // surf cutoff plane coeff Ax+By+Cz=D

  float maxAtomCharge;
  float minAtomCharge;

  int startprinting;          // flag for Movie printing
  SurfEdgeIndex *surfIndex;
  int           numSurfEdges; // number of element in surfIndex

  MoviePrint  *pMovie;

  void ReReadConn();
  void ReadConn();				// allocate resName & list
  void ReadBinXYZ(FILE *fp , bool isSingleRecord=true );
  void ReadPth();
  void ReadDcd();
  void ReadCrd();
  void ReadPdb();
  void ReadXYZ();
  void XYZMaxScale();
  void ReadCoordinate();		

  void CalcNeighbors();			// get neibhbors upon distance and radius
  void CalcDistNeighbors();         // get neighbors upon distance only

  tRGBA  inPickColor;
  bool inPick(int atomIndex, char *pickList, bool getColor=false);

  int  isRibbonProtein(char *pname);
  void setStick();
  void setRibbon();
  void setSpace();

  void initAllocNames(int maxNumMonos, int maxNumAtoms);

  void FindmaxZ();
  void FindmaxAngle();

  void setColor();
  void setSkip();
  void setSphere();
  void setSphereSkip();

  void setSkipSurf();
  void setDefaultDispMode();

  void setAlign();

  void setSecStruct();              // set flags for secondary structure ( secStruct) and chainEnd
  void setSecStructStep1();         //
 
  void getSurface(int isurface);    // write surface i's vertexes into file surfName
  void readSurface();               // read surface vertex into memory
  void markCavity();             // for show Cavity only
  int  setOneOuterSurfAtom();       // ''
  void buildSurfOutConnectivity(int maxzAtom, int maxzEdge);     // ''
  void reorderEdgePair(int j, int k, SurfEdgeIndex  *edgeJ, SurfEdgeIndex *edgeK);
  void setSurfAtomEdgeSelfGroup(int *maxzEdge);
  bool surfEdgeConnected(SurfEdgeIndex *edge1, SurfEdgeIndex *edge2);
  void surfAtomEdgeConnect(int outerEdge, int atomIndex);
  void fixSurfConnectivity();

  bool sameatom(char *a, char *b);	// first atom is the same as second atom
  char* NextOp(char *astring); 	// next operator of input pick string

  char *pdbNameTrim(char *astr, int maxlen);  // remove blanks from begin-end of the string and
     // return the new string.  If maxlen >= 10, then return the input string without changes

  int triangleToDraw(SurfTriangle **cptr, float *cutoffPlane);

  static void DisplayMsgWrtCrd(char *filename);

  void backbone(SplineShape *Loop, int windowSize, int colorFilter);
  void backboneCurve(int windowSize, int colorFilter);
  void backboneTube(int windowSize, int colorFilter);
  void ribbon2D(int windowSize, int colorFilter);
  void ribbon3D(int windowSize, int colorFilter);
  void struct2nd2D(int windowSize,int colorFilter);  
  void struct2nd3D(int windowSize,int colorFilter);  

  void quadStripStruct(SplineShape *splineShape, int windowSize, int colorFilter, bool isRibbon);  // draw ribbon or 2nd struct

public:
   AtomArray();
  ~AtomArray();

  InputParameters InParm;
  int   dispmodeComb;         // dispmode combination for all atoms
  
  int 		numRes;
  int 		numAtom;
  char 		moleculeName[MAXLEN_NAMES];   // name in the beginning of wcon file
  NameString 	*resName;				// [MONO_MAXNUM];
  Atom 		*atomList;				// [ATOM_MAXNUM];
  bool            orgCrdModified;               // =true only if orgX,Y,Zs in atomList been modified

  tRGBA           *ribbonColor;                 // index of rn
     
  SurfHeader      surfHeader;                   //
  SurfTriangle	*surface;                     // 
  int             nSurfPlanePoints;             // select points for the surf cutoff plane
  float           totalSubArea;                 // total area on selected sub surfaces
  int             numSubAreaAtoms;              // number of atoms on the selected area
  void            setSurfCutoffPlane(int isAxe=0);   // get last 3 picked atoms for cutoff plane
  void            flipSurfCutoffPlane();        // reverse normal vector

  int             currentmol;			      // current reading in molecular
  int             currentsurf;                  // current reading surface vertexes

  int             numChains;                    // for PDB reading   
  Chain           *chains;  
  int             displayChainid; 
  // int             displayModelid;      

  float           OverlapRmsd;

  void            getSurfaces();               // set surface for all structures
  void            getSubArea();
  void  		getSurfAtoms();

  void SubtractCoMass();
  void SubtractCoMass(double cX, double cY, double cZ);
  void ResetCrd();

  double 	maxZ;
  double 	minZ;
  double 	angle;

  double CoMX, CoMY, CoMZ;	// for SubtractCoMass()
  int    pickSphereCenter;    // atom index of pick sphere center;

  //double 	min[3];	/* minX, minY, minZ: for XYZ display */
  //double 	max[3];     /* maxX, maxY, maxZ: for XYZ display */
  
  void 	Print();

  double Distance(Atom *atom1, Atom *atom2);  // calculate the distance between two atoms

  /*DisplayMode*/ int DispMode();    // allow bit combination
  int BondMode();

  void ReadNextStruct();
  void DisplayNextPdbModel();	
  void DisplayNextPdbChain();	

  int pick_atoms(int x, int y);
  bool checkNeigh(int no, char *a, int *sn);
  int setNeigh();

  int GetInputParameters (char *inputfile, int crdIndex=0, char *cmoil=NULL);	// return crd index
  void ProcessData();
  void ReProcessData();   // default to reprocess without overlap matrix

  void mprint_start(int PrintOneFrameOnly=0);     // movie print
  void mprint_frame();

  // manage chains
  void ChainAdd(char chainid, int resStart, int resMap, int atomStart, int atomMap, int model=-1);
  int  resNum (Atom *anatom);      // get display resNumber
  int  atomIndex(int pdbAtomNum, int *mapArray);  // convert pdb atom number to internal index

  void Connect(int atom1, int atom2, int checkFirst=0);   // connect two atoms indicated by the atom index 
  double AtomRadius(char *atomName);   // get estimated radius from the input atom name 

  //Coordinate *quad(Coordinate CA, Coordinate O, int reset, double ribbonWidth=1.0);   // get ribbon middle point
  void backbone (int windowSize=2,int colorFilter=ColorFilterNone);      // draw backbone
  void struct2nd(int windowSize=2,int colorFilter=ColorFilterNone);   // draw 2nd structure,  helix, sheet, etc  or a ribbon
  void ribbon   (int windowSize=2,int colorFilter=ColorFilterNone);   // draw ribbon

  void surf(int colorFilter=ColorFilterNone);           // draw Connolly Surface
  void WriteOrgCrd(char *extraTxt="");                  // write current structure to a file in CRD forma
  void OutputCrdInPDB(FILE *crdf, int chainIndexID=0, char *extraTxt="");

  void SphereNet(int colorFilter=ColorFilterNone);

  void DisplayStructure(int &imgIndex, int colorFilter, int lineWidth=1);

};


}

#endif
