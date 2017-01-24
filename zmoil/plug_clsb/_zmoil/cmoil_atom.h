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
#include "cmoil_globals.h"

#include "zvec.h"
	// for DMat4, DVec3 in ApplyTransform()
#include "zhashtable.h"
#include "ztl.h"
#include "ztime.h"
#include "zhash3d.h"

extern void trace( char *msg, ... );
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

	// optional velocity of atom (read from dvd file)
  double vx;
  double vy;
  double vz;
  double velocityMag() {
	  return sqrt( vx * vx + vy * vy + vz * vz );
  }

  int rn;			// residue number
  int rnPDB;		// residue number to display or the real PDB resNum, for supporting PDB file
  int chainIndex;		// chain index of the chain array
  NameString atomName;  // discard atom if atomName=""
  int neighbors[ATOM_MAX_NEIGHBORS];
  short int numNeighbors;
  int isNeighbor( int atom ) {
	  for( int i=0; i<numNeighbors; i++ ) {
		  if( atom == neighbors[i] ) {
			  return 1;
		  }
	  }
	  return 0;
  }

  int hydrogenBonds[ATOM_MAX_HBONDS];
  int numHBonds;
	// july 2011: hydrogen bonding

  tRGBA color;
  float charge;
  double weight;
	// weight added may 2011; is the 'weighting' or 'b-factor' column of CRD file

  bool skip;
  bool skipsurf;        // for generating surface
  bool skipsurfdisp;    // for displaying surface
  bool align;		// selected atoms for overlap
  int  dispmode;
  SecStructType secStruct;    // type of 2nd structure

  const Atom &operator=(const Coordinate &a);
  friend double angleT(Atom &a0, Atom &a1, Atom &a2, Atom &a3) ;
};


class AtomArray;
class Model {
public:
	int    modelName;			// the value from the pdb MODEL entry - not our index in models[]!
	int    resPdbStart;			// the start residue number in a PDB file
	int    resMapStart;			// the start residue index in the resName array
	int    atomPdbStart;        // the start atom numver in a PDB file
	int    atomMapStart;        // the start atom index in Atom array
	int    display;				// should this model be displayed?
	Model() { 
		display=1;
	}
};

class Chain {
public:
	char   id;				  // 22th char of a PDB Atom record
	int    modelIdx;			  // index of the model in models[]
	int    resPdbStart;		  // the start residue number in a PDB file
	int    resMapStart;		  // the start residue index in the resName array
	int    atomPdbStart;        // the start atom numver in a PDB file
	int    atomMapStart;        // the start atom index in Atom array
	int    display;			  // should this chain be displayed?
	Chain() { 
		display = 1;
	}
};

class AtomArray : public State {
  friend bool linePlaneCrossPoint(float *plane, Coordinate &pt1, Coordinate &pt2, float *crossPt);

public:
  // format of binary surface vertex file
  typedef struct _SurfHeader
  {
    int numSurfTris;		      // number of triangles in surface
    short int outerProced;		// outer surface processed
    short int cavityNormalReversed;   // reverse normal vector for cavity atoms
  } SurfHeader;         		// binary vertex file format: maxNumSurfTris + [header+SurfTriangle*numSurfTris]*M

  struct SurfTriangle;
  struct CavityInfo {
	  double surfaceArea;
	  ZHashTable atoms;
		// a hash aids clustering the atoms that belong to the cavity
	  ZTLVec< SurfTriangle * > tris;
	  CavityInfo() {
		  surfaceArea = 0.0;
	  }
  };
  ZTLPVec< CavityInfo > cavities;

  // Spatial Hashing: at present we don't require information about who is in the same
  // box, we just want to know, for a given coordinate, what it's address is so that we
  // can use this address (3-tuple) to look up other information that pertains to this box.
  // So use of the ZHash3d is perhaps overkill, but does not hurt us.
  // TODO: we might make use of the spatial hash for the cavity-functionality above?
  //
  // NOTE: I have made all of this static such that one boxinfo file applies to all
  // loaded structures.  We'll make it per-structure later if it makes more sense.
  static ZHash3D spatialHash;
  static void initSpatialHash( int xBins, int yBins, int zBins, FVec3 &min, FVec3 &max );
  static double *boxInfo;
  static double boxInfoMin;
  static double boxInfoMax;
	// boxInfo will be alloc'd for gridsize^3 entries, and will be used to lookup
	// the data for the box based on the index-address of a 3-space coordinate.
  static void getColorFromBoxInfo( float x, float y, float z, tRGBA *c );
  static int readBoxInfo( char *filename );
	// read dimensions & values from a file generated by MOIL

	

  struct SurfTriangle      // line structure of *surface triangle
  {
    int   atomNum;
    float nm[3][3];                // normals
    float vt[3][3];                // vertexes
    int   color;                  
    int   nextAtom;                // neighbor atom number
    int   edge;
    int   out;                     // set to 1 if it's on outer surface, not cavity   
	CavityInfo* cavity;					// used to group discrete cavities
	SurfTriangle() {
		cavity = 0;
	}
	void dump() { trace( "atomNum %d, nextAtom %d, edge %d, cavity %X\n", atomNum, nextAtom, edge, cavity ); }
  };

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

  FILE *crdfp;				  // ReadPth() , Dcd, Crd
  FILE *velfp;				  // for option dvd files
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

  void ReReadConn();
  void ReadConn( int createBonds );				// allocate resName & list
  bool ReadBinXYZ(FILE *fp , bool isSingleRecord=true, bool isVelocity=false );
  void ReadPth( int index = -1 );
  void getVelocityStats();
  void buildSpatialHash( float gridsize );
  void ReadDcd( int index = -1 );
  void ReadCrd();
  void ReadPdb();
  void ReadXYZ();
  void XYZMaxScale();
  void ReadCoordinate( int index=-1 );		

  void CalcNeighbors();			// get neibhbors upon distance and radius
  void CalcDistNeighbors();     // get neighbors upon distance only

  void calcHydrogenBonds();		// Look for hydrongen bonding

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
public:
		// public for debug code
  void setSkip();
  void setSphere();
  void setSphereSkip();
   int WriteTransformedPdb();


private:

  void setSkipSurf();
  void setDefaultDispMode();

  void setAlign();

  void setSecStruct();              // set flags for secondary structure ( secStruct) and chainEnd
  void setSecStructStep1();         //
 
  int  getSurface(int isurface);    // write surface i's vertexes into file surfName
  void readSurface();               // read surface vertex into memory

  void readRaviSurfaceNormals();
  void readRaviAtoms( char *ptr );

  public:
  void printCavity();
  private:
	int clustersOverlap( int c1, int c2 );
	void collapseClusters();
	CavityInfo * addToCavityCluster( int atom1, int atom2 );
	void clusterCavities();
	int surfTrianglesConnected( SurfTriangle *t1, SurfTriangle *t2 );
	void partitionClusterByGeometry( CavityInfo *ci );

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

  void backbone(SplineShape *Loop, int windowSize, int colorFilter, int type );
  void backboneCurve(int windowSize, int colorFilter, int type );
  void backboneTube(int windowSize, int colorFilter, int type );
  void ribbon2D(int windowSize, int colorFilter);
  void ribbon3D(int windowSize, int colorFilter);
  void struct2nd2D(int windowSize,int colorFilter);  
  void struct2nd3D(int windowSize,int colorFilter);  

  void quadStripStruct(SplineShape *splineShape, int windowSize, int colorFilter, bool isRibbon);  // draw ribbon or 2nd struct

public:
   AtomArray();
  ~AtomArray();

  InputParameters InParm;
  int   dispmodeComb;							// dispmode combination for all atoms

  int alpha;									// a global alpha to be applied to all draw operations
												// for use when displaying multiple molecules (tfb)
  bool bDraw;									// used to disable drawing of this structure
  bool bBlink;									// model will alternate draw/nodraw
  bool bCycle;									// model will cycle through alternative MODEL def'ns
  ZTime lastBlinkOn;							// use to control blink rate
  ZTime lastStructAnimate;						// use to control animation rate

  ZHashTable properties;
		// rather than continue to add binary members to AtomArray
		// I'm going to start storing properties here which will make
		// loading and saving easier should we chose to do that.

  enum CoMUpdateFrequency { CoM_Never, CoM_First, CoM_Always };
  CoMUpdateFrequency updateCoM;					// when should CoM be auto-updated?
  
  int 		numRes;
  int 		numAtom;
  char 		moleculeName[256];			// name in the beginning of wcon file
  NameString 	*resName;						// [MONO_MAXNUM];
  Atom 		*atomList;							// [ATOM_MAXNUM];
  int		firstVisibleAtom;
  int		lastVisibleAtom;
	// optimization for avoiding iteration over all atoms when whole 
	// models or chains are turned off.

  bool            orgCrdModified;               // =true only if orgX,Y,Zs in atomList been modified

  tRGBA           *ribbonColor;                 // index of rn
     
  SurfHeader      surfHeader;                   //
  SurfTriangle	 *surface;						// 
  int             nSurfPlanePoints;             // select points for the surf cutoff plane
  float           totalSubArea;                 // total area on selected sub surfaces
  int             numSubAreaAtoms;              // number of atoms on the selected area
  void            setSurfCutoffPlane(int isAxe=0);   // get last 3 picked atoms for cutoff plane
  void            flipSurfCutoffPlane();        // reverse normal vector

  int			   raviAtomCount;	
  double		  *raviAtoms;
  double          *raviPoints;
  double	      *raviNormals;
  double		  *raviNormalsCoarse;
  int			  raviNormalCount;
					// this stuff above is experimental surface normal information
					// used by Ravi for his docking algorithms.  These normals are
					// loaded from a separate file, and have nothing to do with
					// the normals of the surface generated by this program.  See
					// the function readRaviSurfaceNormals() to learn how the text
					// file is parsed to find the vectors to be drawn.

  int             currentmol;					// current reading in molecular
  long			 *structOffsets;				// holds offsets to structs in pthfile for random access (Tfb)
  long			 *surfOffsets;					// holds offsets to surfaces in .vtx file for random acces (tfb)

  int           numChains;                    // for PDB reading   
  Chain         *chains;  
  int			numModels;
  Model			*models;

  float           OverlapRmsd;

  void            getSurfaces();               // set surface for all structures
  void            getSubArea();
  void  		  getSurfAtoms();

  void UpdateBounds();

  double    maxDisplacementXYZ;
  DVec3		boundsMin;
  DVec3		boundsMax;
	// vecs define a bounding box that contains all atom coords
  double 	angle;


  DVec3 translate;
	// tracks the cumulative translation of atomic coordinates from their original
	// values as specified in the crd file.

  DVec3 centroid;
	// previous referred to as "center of mass"; the arithmetic mean of the original
	// coordinates read from the pdb/crd file; does not include 'skipped' atoms - used
	// to center the molecule onscreen, etc.

  void resetCrd();
	// reset the atomic coordinates to those as specified in the crd file.
	// clears translate member.  
	
  void computeCentroid();
	// compute mean coordinate of all displayed atoms

  void applyTranslation( DVec3 &t, bool negate=false );
	// apply the translation to atomic coordinates; translate member is adjusted
	// by this amount.  Typically used to translate by -centroid to center molecule.
	// If you desire the translation to be absolute, you may need to call ResetCrd
	// first to ensure any previous translations are reset.  optional negate for 
	// common practice of wanting to subtract a vector (e.g. centroid)

   void applyTransform( DMat4 & mat );
	// destructively apply the transform to the x,y,z coords of each atom.  
	
  int    pickSphereCenter;    // atom index of pick sphere center;

  AtomArray *alignToModel;
	// when this is nonNull, it indicates we should align ourselves to this structure.  
  DMat4 alignMatrix;
	// the matrix provided by tmalign for aligning our structure to the target
  char alignSequence[1024];
  bool setAlignmentTarget( AtomArray *aaptr, int trivial );
	// align our atoms to this target; this is currently done with external program tmalign,
	// which gives us a rotation matrix; this is stored in alignMatrix above, and we 
	// transform our atoms by this for display (in addition to the COM translation of
	// the target.
  void setPickedResidue( int resIndex );
	// used to highlight a residue that has been picked from UI

  // User transform per model: to allow models to be moved independently of one another,
  // we allow a per-model viewing transform.  This is *NOT* destructively applied to
  // atoms like the alignment transforms above are; this may present a problem in distance
  // metrics since the distance functions (and probably other things) make use of the atomic x,y,z coords;

	FQuat viewQuat;
		// the rotation in the coordinate system of this model
	FVec3 viewTranslate;
		// the translation in the world coordinate system

	void clearAlignments();
		// clear the alignToModel as well as the viewQuat/viewTranslate

  // -----------------------------------------------------------------------------------
  void 	Print();

  double Distance(Atom *atom1, Atom *atom2);  // calculate the distance between two atoms

  /*DisplayMode*/ int DispMode();    // allow bit combination
  int BondMode();

  void ReadNextStruct();
  void ReadStruct( int index );
//  void DisplayNextPdbModel();	
//  void DisplayNextPdbChain();	
  void setModelChainDisplay();
  void DisplayPdbModel( int model, bool show );
  void DisplayNextPdbModel();
  void DisplayPdbChain( int chain, bool show );

  bool checkNeigh(int no, char *a, int *sn);
  int setNeigh();

  int GetInputParameters (char *inputfile, int crdIndex=0, char *cmoil=NULL);	// return crd index
  void ProcessData();
  void ReProcessData();   // default to reprocess without overlap matrix

  void mprint_start(int PrintOneFrameOnly=0);     // movie print
  void mprint_frame();

  // manage chains
  void ChainAdd(char chainid, int resStart, int resMap, int atomStart, int atomMap, int model );
  void ModelAdd(int resStart, int resMap, int atomStart, int atomMap, int modelName );

  int  resNum (Atom *anatom);
	// get display resNumber
  int  atomIndex(int pdbAtomNum, int *mapArray);
	// convert pdb atom number to internal index
  char * getResName( int atomIndex ); 
	// get residue name from atom index

  int Connect(int atom1, int atom2, int checkFirst=0);   // connect two atoms indicated by the atom index 
  void connectWater();
	// checks water residues and generates bonds if not already bonded

  double AtomRadius(char *atomName);   // get estimated radius from the input atom name 

  //Coordinate *quad(Coordinate CA, Coordinate O, int reset, double ribbonWidth=1.0);   // get ribbon middle point
  void backbone (int windowSize=2,int colorFilter=ColorFilterNone, int type=BACKBONE );      // draw backbone
  void struct2nd(int windowSize=2,int colorFilter=ColorFilterNone );   // draw 2nd structure,  helix, sheet, etc  or a ribbon
  void ribbon   (int windowSize=2,int colorFilter=ColorFilterNone );   // draw ribbon
  void peptidePlaneVectors();	
	// debug: draw the vectors used to compute the normals to the peptide plane
  void velocityVectors();
  void velocityColors();

  void surf(int colorFilter=ColorFilterNone, tRGBA *defaultColor=0 );           // draw Connolly Surface
  void WriteOrgCrd(char *extraTxt="");                  // write current structure to a file in CRD forma
  void OutputCrdInPDB(FILE *crdf, int chainIndexID=0, char *extraTxt="");

  void SphereNet(int colorFilter=ColorFilterNone);

  void DisplayStructure(int &imgIndex, int colorFilter, int lineWidth=1);
  void DisplayPickedResidue( tRGBA &c, int colorFilter );

  void raviNormalDraw( int colorFilter=ColorFilterNone, bool bCoarseNormals=false );
  void raviAtomDraw( int colorFilter=ColorFilterNone );

};


}

#endif
