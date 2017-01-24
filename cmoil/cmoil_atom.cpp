//*************************************************************************************************
//*  Filename:   cmoil_atom.cpp
//*
//*  Description: 
//*    Methods for Atom and AtomArray classes
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. ?? 2000	??		Initial Development
//*  Jul. 18 2001	Baohua Wang Seperatd from cmoil.cpp. Reshape to support mutiple structure display 
//*              	           	and overlap.
//*
//*************************************************************************************************
#include <ctype.h>
#include "cmoil.h"
#include "cmoil_atom.h"
#include "cmoil_dib.h"
#include "cmoil_aminoacids.h"
#include "cmoil_pick.h"
#include "cmoil_rms.h"

//#include <iostream>


namespace CMOIL {

extern bool linePlaneCrossPoint(float *plane, Coordinate &pt1, Coordinate &pt2, float *crossPt);

double AtomArray::BOND_MAXLEN = 1.90;  //1.645 for some H-H distance=1.648  while 1.9 for S-C bond 1.85
double AtomArray::BOND_MINLEN = 0.85; 

/**************************************************************************************************/
//here are the declarations for all the necessary global variables...

Atom::Atom()
{
  numNeighbors=0;
  skip=0;				// false
  skipsurf=0;			// false
  skipsurfdisp=0;
  dispmode=SKIP;
  align = false;

  chainIndex = -1 ;           // 0;	            // chain index of multiple chain structure
  secStruct=STRUCT_UNKNOWN;
  charge=0.0;
}

const Atom &Atom::operator=(const Coordinate &a)
{
  x=a.x;
  y=a.y;
  z=a.z;
  return *this;
}

AtomArray::AtomArray(void)
{
  currentmol = 0;
  currentsurf=0;
  pMovie = NULL;
  CoMX = CoMY = CoMZ = 0.0;
  pickSphereCenter=-1;             // no picksphere
  orgCrdModified=false;
  crdfp=vertexfp=NULL;
  resName  = NULL;
  atomList = NULL;
  ribbonColor=NULL;
  surface  = NULL;
  surfIndex = NULL;
  totalSubArea=0.;
  numSubAreaAtoms=0;
  surfHeader.numSurfTris = 0;
  surfHeader.outerProced=0;
  surfHeader.cavityNormalReversed=1;
  maxNumSurfTris=0;
  nSurfPlanePoints=0;
  strcpy(moleculeName, "BULK");    // for some PDB without header name

  numChains = 0;
  chains = NULL;
  displayChainid = -1 ; // 0;     // index of chains, display all chains
  // displayModelid = 0;          // display all models
  
  dispmodeComb=0;                 // display mode bit combination

  dcdNFrzAtoms=NULL;
  nNoFrzAtoms=0;

  setSurfCutoffPlane(0);
}

AtomArray::~AtomArray(void)
{
  if (pMovie != NULL )
    delete pMovie;
  if (resName != NULL)
    delete[] resName;
  if (atomList != NULL)
    delete[] atomList;
  if (ribbonColor != NULL)
    delete[] ribbonColor;
  if ( dcdNFrzAtoms != NULL )
    delete[] dcdNFrzAtoms;
  if ( surface != NULL )
    delete[] surface ;
  if (  surfIndex != NULL)
    delete[]  surfIndex; 
  if (crdfp != NULL)
    fclose(crdfp);
  if (vertexfp != NULL)
    fclose(vertexfp);
}

void AtomArray::ChainAdd(char chainid, int resStart, int resMap, int atomStart, int atomMap, int modelid)
{
  Chain *newchains=NULL;

  if ( numChains %10 == 0 )
  {
    newchains = new Chain[numChains+10];		// allocate 10 chain at once
    if ( newchains == NULL)
    {
      fprintf(stderr, "CMOIL ERROR: uable to allocate memory for chians!\n");
      exit(-1);
    }

    if (numChains > 0 )
    {
      memcpy(newchains, chains, sizeof (Chain)*numChains);
      delete[] chains;
    }
    chains = newchains;
  }

  newchains = chains+numChains;    // move to the new element

  newchains->id = chainid;
  newchains->model = modelid;          //  default  is  0 
  newchains->resPdbStart=resStart;
  newchains->resMapStart=resMap;
  newchains->atomPdbStart=atomStart;
  newchains->atomMapStart=atomMap;
  numChains++;
}

// residure number to display
int AtomArray::resNum(Atom *myatom)
{
  if (this->InParm.filetag != fPDB )         // not reading from a pdb file
    return (myatom->rn+1);
  else if (myatom != NULL)
    return (myatom->rnPDB);		// display real PDB residue number
  else
    return 0;
}

void AtomArray::mprint_start(int oneFrameOnly)
{
  if ( InParm.structpstart>0 && currentmol>InParm.structpstart)
  {
    // rewind to the first picture
    currentmol=InParm.structno-1;
    startprinting=0;
  }
  if (pMovie == NULL) 
    pMovie = new MoviePrint(InParm);
  pMovie->start_print(oneFrameOnly);
}

void AtomArray::mprint_frame( )
{
  if  ( InParm.structpstart>0 && currentmol==InParm.structpstart )
  {
     startprinting=1;
  }

  if (pMovie != NULL && 
      (InParm.structpstart<=0 || startprinting) )
  {
    pMovie->print();
  }
}

const unsigned short int FORTRAN_INTSIZE = 4;		// FORTRAN integer size

void AtomArray::ReadPth()
{
  float dummy(0);
 
  // fprintf(stdout, "\n");
	// why?  (tfb removed oct4_2007)

  if (currentmol >= InParm.structno) 
    currentmol = 0;
  
  if (currentmol == 0) 
  {
    if ( crdfp != NULL)
      rewind(crdfp);
    else if ((crdfp = fopen(InParm.pthName, "rb")) == NULL) 
    {
      printf("Error opening file\n");
      exit(-1);
    } 
  } 
  currentmol++;

  fread(&dummy, sizeof(float),1,crdfp);
  fread(&dummy, FORTRAN_INTSIZE, 1,crdfp);

  ReadBinXYZ(crdfp, true);
}	// ReadPth()

void AtomArray::ReadBinXYZ(FILE *fp, bool isPth) 
{
   int  counter;
   int  crd_size;
   int  arraySize;
   char dummy[FORTRAN_INTSIZE+FORTRAN_INTSIZE];
   double element;

   if ( fp == NULL || feof(fp) )
       return;

   fread(&arraySize, FORTRAN_INTSIZE, 1,fp);
   if (!isPth )
   {
     crd_size = sizeof(float);
     if (InParm.swapmode == 1 )
     	 InParm.SwapBite(&arraySize, FORTRAN_INTSIZE);
     arraySize /= crd_size; 
   }
   else
   {
      crd_size=sizeof(double);
      arraySize=numAtom;
   }

   for (counter=0; counter<arraySize ; counter++)
   {
      fread(&element, crd_size ,1,fp);
      if (InParm.swapmode)    
         InParm.SwapBite(&element, crd_size );

      if (isPth )
         atomList[counter].x= element ;
      else 
      {
         if ( arraySize==numAtom)
           atomList[counter].x = *((float*)&element);
         else if (dcdNFrzAtoms != NULL )
           atomList[dcdNFrzAtoms[counter]-1].x=*((float*)&element);
      } 
   }

   if ( !isPth)
      fread(dummy, FORTRAN_INTSIZE+FORTRAN_INTSIZE, 1,fp);
        
   for (counter=0; counter<arraySize; counter++)
   {
      fread(&element, crd_size ,1,fp);
      if (InParm.swapmode)    
         InParm.SwapBite(&element, crd_size );

      if (isPth )
         atomList[counter].y= element ;
      else 
      {
         if ( arraySize==numAtom)
           atomList[counter].y = *((float*)&element);
         else if (dcdNFrzAtoms != NULL )
           atomList[dcdNFrzAtoms[counter]-1].y=*((float*)&element);
      } 
   }

   if ( !isPth)
      fread(dummy, FORTRAN_INTSIZE+FORTRAN_INTSIZE, 1,fp);

   for (counter=0; counter<arraySize; counter++)
   {
      fread(&element, crd_size ,1,fp);
      if (InParm.swapmode)    
         InParm.SwapBite(&element, crd_size );

      if (isPth )
         atomList[counter].z= element ;
      else 
      {
         if (arraySize==numAtom)
           atomList[counter].z = *((float*)&element);
         else if (dcdNFrzAtoms != NULL )
           atomList[dcdNFrzAtoms[counter]-1].z=*((float*)&element);
      } 
   }

   /* Read a final dummy */
   fread(dummy, FORTRAN_INTSIZE,1,fp);

}	// ReadBinXYZ()

void AtomArray::ReadDcd()
{
   unsigned char dummyStr[256];
   int    nFixedAtoms=0;
   int    nptx = 0;
   int    record_size=0;
   int    i;
   Atom   *listi;
     
   // fprintf(stdout, "\n");
			// why? (tfb oct4_2007)

   if (currentmol >= InParm.structno)
     currentmol=0;

   if (currentmol == 0) 
   {
      if (crdfp != NULL)
        rewind(crdfp);
      else if ((crdfp = fopen(InParm.dcdName, "rb")) == NULL) {
         fprintf(stderr, "Error opening file\n");
         exit(-1);
      } 
      fread(&record_size, FORTRAN_INTSIZE, 1, crdfp) ;            // read record size(array size)
      if (InParm.swapmode == 1 )
     	   InParm.SwapBite(&record_size, FORTRAN_INTSIZE);

      fread(dummyStr, record_size+FORTRAN_INTSIZE, 1, crdfp) ;    // head, vari[20s]

      nFixedAtoms = *((int*)(dummyStr + FORTRAN_INTSIZE*9 ));	// set nFixedAtom(vari[9th])
	if (InParm.swapmode == 1 )
     	   InParm.SwapBite(&nFixedAtoms, FORTRAN_INTSIZE);
  
      fread(&record_size, FORTRAN_INTSIZE, 1, crdfp) ; 
      if (InParm.swapmode == 1 )
     	   InParm.SwapBite(&record_size, FORTRAN_INTSIZE);

      fread(dummyStr, record_size+FORTRAN_INTSIZE, 1, crdfp); 		// i, title

      //printf("5. read in array of %d\n", record_size+FORTRAN_INTSIZE);

      fread(&record_size, FORTRAN_INTSIZE, 1, crdfp) ;   		// record_size = FORTRAN_INTSIZE
      if (InParm.swapmode == 1 )
     	   InParm.SwapBite(&record_size, FORTRAN_INTSIZE);

      if (record_size != FORTRAN_INTSIZE )
      {
         printf("Incorrect DCD file format\n");
         exit(-1);
      }

      fread(&nptx, FORTRAN_INTSIZE, 1, crdfp); 	      		// nptx(Int) = 04->93
      if (InParm.swapmode == 1 )
         InParm.SwapBite(&nptx, FORTRAN_INTSIZE);

      if (nptx != numAtom )
      {
         fprintf(stderr,  "Incorrect number of particles in the DCD file: %d\n", nptx);
         exit(-1);
      }

      fread(&record_size, FORTRAN_INTSIZE, 1, crdfp) ; 

      nNoFrzAtoms = numAtom - nFixedAtoms;

      // printf("7, now : numAtom=%d, nFixedAtoms=%d, nNofrz=%d\n", numAtom, nFixedAtoms,nNofrz );

  		if( dcdNFrzAtoms != NULL ) {
			delete[]dcdNFrzAtoms;
			dcdNFrzAtoms = 0;
		}
		if( nNoFrzAtoms>0 && nNoFrzAtoms != numAtom ) {	// read in list of nNofrz atom indeces
			dcdNFrzAtoms = new int[nNoFrzAtoms];
			fread( &record_size, FORTRAN_INTSIZE, 1, crdfp );
			if( InParm.swapmode == 1 )
				InParm.SwapBite( &record_size, FORTRAN_INTSIZE );

			for( i = 0; i < nNoFrzAtoms; i++ ) {
				fread( &record_size, FORTRAN_INTSIZE, 1, crdfp );
				if( InParm.swapmode == 1 )
					InParm.SwapBite( &record_size, FORTRAN_INTSIZE );
				dcdNFrzAtoms[i] = record_size;	// store atom index
			}
			fread( &record_size, FORTRAN_INTSIZE, 1, crdfp );
			if( InParm.swapmode == 1 )
				InParm.SwapBite( &record_size, FORTRAN_INTSIZE );

		}
   } 
   else if ( nNoFrzAtoms < numAtom)
   {      // keep orginal crd
     for (i=0, listi=atomList; i<numAtom; i++, listi++)
     {
       listi->x=listi->orgX;
       listi->y=listi->orgY;
       listi->z=listi->orgZ;
     }
   }

   ReadBinXYZ(crdfp, false);	

   // copy to orgXYZ
   if ( currentmol > 0 && nNoFrzAtoms<numAtom )
   {
     for (i=0; i<nNoFrzAtoms; i++)
     {
       listi=atomList+dcdNFrzAtoms[i];
       listi->orgX=listi->x;
       listi->orgY=listi->y;
       listi->orgZ=listi->z;
     } 
   }

   currentmol++;
}	// ReadDcd()

void AtomArray::ReadCrd()
{
  int  counter;
  char aline[LINE_MAXLEN];
  Atom *myatom;

  // crd only has one structure, thus InParm.structno=currentmol-1
  currentmol=1;

  if ( crdfp != NULL)
      rewind(crdfp);
  else if ((crdfp = fopen(InParm.crdName, "r")) == NULL) 
  {
    fprintf (stderr, "Unable to open file file : %s\n", InParm.crdName );
    exit(-1);
  }

  for (counter=0; counter<numAtom && fgets(aline,LINE_MAXLEN,crdfp)!= NULL ; )
  {
     myatom = atomList+counter; 
     if ( aline[0] != '*' &&  strlen(aline) >= 50 ) 
     {	// coordinates starts at 21st ends at 51st
	 sscanf(&aline[20], "%10lf%10lf%10lf", &(myatom->x),&(myatom->y),&(myatom->z));
       counter++;
     }
  }
}  //ReadCrd()


// read coordinates without read in wcon // 
// const unsigned int ATOM_MAXNUM  	=   100000;

char *AtomArray::pdbNameTrim(char *astr, int maxlen)
{
  static char myName[11];
  char *cptr;

  if ( astr == NULL)
    return NULL;

  if ( maxlen >= 10 )
    return astr;
 
  strncpy(myName, astr, maxlen);
  cptr= myName + maxlen;
  *cptr = '\0';			/* NULL terminated */
  
  cptr--;
  for ( ; cptr != myName && (*cptr == ' ' || iscntrl(*cptr)) ; cptr--)
    *cptr = '\0';
  
  for (cptr = myName; (*cptr==' '|| iscntrl(*cptr)); cptr++ )
    ;
  return cptr; 
}

#include "math.h"
double AtomArray::Distance(Atom *atom1, Atom *atom2)
{
  double tmp, sum=0.0;
 
  tmp =atom1->x-atom2->x;
  sum += tmp*tmp;

  tmp =atom1->y-atom2->y;
  sum += tmp*tmp;

  tmp =atom1->z-atom2->z;
  sum += tmp*tmp;

  return sqrt(sum);
}

// get estimated radius for input atom name, radius = PSGM/2.0 in all.prop
double AtomArray::AtomRadius(char *atomname)
{
   if ( strstr(atomname, "CL") != NULL  )
     return 1.50;
   else if ( strstr(atomname, "FE") != NULL )
     return 0.58;
   else if ( *atomname == 'N')
     return 1.62;
   else if ( *atomname == 'C' )
     return 1.90;
   else if ( *atomname == 'O' )
     return 1.48;
   else if ( strchr(atomname, 'P') != NULL )
     return 1.87;
   else if ( strchr(atomname, 'S') != NULL )
     return 1.77;
   else if ( *atomname == 'H' || ( atomname[1] == 'H' && !isalpha(*atomname) ) )
     return 0.15;
   else
     return 1.0;           // for any unknowns
}

// default: maxNumMonos=MONO_MAXNUM, maxNumAtoms=ATOM_MAXNUM);
void AtomArray::initAllocNames(int maxNumMonos, int maxNumAtoms)
{
  if ( resName == NULL )
  {
    resName = new NameString[maxNumMonos]; 
    if ( resName == NULL )
    {
      fprintf(stderr,"Uable to allocate NameString[%d]\n", maxNumMonos);
      exit(-1);
    }
    memset(resName, 0, sizeof(NameString)*maxNumMonos);
    numRes = maxNumMonos;                     
  } 
   
  if (atomList == NULL)
  {
    atomList  = new Atom[maxNumAtoms];   // default to ATOM_MAXNUM
    if ( atomList == NULL)
    {
      fprintf(stderr, "Uable to allocate Atom[%d]\n", maxNumAtoms);;
      exit(-1);
    }
    memset(atomList, 0, sizeof(Atom)*maxNumAtoms);
    numAtom = maxNumAtoms;
  }

  if ( InParm.ribboncolor  != NULL  )     // add to readPdb also
  {
    ribbonColor = new tRGBA[numRes];
    if ( ribbonColor == NULL)
    {
      fprintf(stderr, "Unable to allocate memory!\n" );
      exit(-1);
    }
  }
}

void AtomArray::ReadPdb()
{
  // since in PDB file only one structure exists, the file will read only once
  currentmol = 1;
  if ( crdfp != NULL)
    fclose(crdfp);
  
  if ((crdfp = fopen(InParm.pdbName, "r")) == NULL) 
  {
    fprintf(stderr, "Unable to open file to read : %s", InParm.pdbName);
    exit(-1);
  }
  initAllocNames(MONO_MAXNUM, ATOM_MAXNUM);
  if (chains != NULL)      // clean up the chains structure
  {
     delete[] chains;
     chains=NULL;
     numChains=0;
  }

  char        aline[LINE_MAXLEN]="";
  // read in PDB fields upon PDB format definition
  char *atom_name = aline + 12;    // len = 4 
  char *res_name  = aline + 17;    // len = 4
  char *alt_loc   = aline + 16;    // len = 1
  char *atom_num  = aline + 6;     // len = 5
  char *res_num   = aline + 22;    // len = 4
  char *insert    = aline + 26;    // len = 1
  char *xyz       = aline + 30;    // 3 floats
  char *charge    = aline + 78;    // len = 2

  int  resIndex=-1, atmIndex=-1;   // all atoms indexed together, residue indexed by chain
  int  modelid = 0;			// real model id start at 1 
  char lastInsert=' ';             
  char chainid='\0';
  int  pdbAtomNum, fisrtAtomIndexOfRes;
  int  pdbResNum,  lastPdbResNum=-999; 

  NameString *myres=NULL;
  Atom       *myatom=NULL; 
  int anPDB[ATOM_MAXNUM];      // only need for build connectivity

  while (fgets(aline,LINE_MAXLEN,crdfp)!= NULL )
  { 
    if (strncmp(aline, "ATOM", 4)==0 || strncmp(aline, "HETATM", 6)==0) 
    {
      // new atoms
      atmIndex++;
      myatom = atomList + atmIndex;
      pdbAtomNum = atoi(pdbNameTrim(atom_num, 5));    // number in PDB file
      pdbResNum = atoi(pdbNameTrim(res_num, 4));      // number in PDB file
      anPDB[atmIndex]=pdbAtomNum;
       
      // new residures
      if ( lastPdbResNum != pdbResNum ||
           lastPdbResNum==pdbResNum &&  *insert != ' ' && *insert != lastInsert )
      {
        resIndex++;
        myres  = resName + resIndex;    
        strcpy((char*)myres, pdbNameTrim(res_name, 4));	 // new residure 
        fisrtAtomIndexOfRes=atmIndex;
        lastInsert = *insert;
        lastPdbResNum = pdbResNum;
      }

      // handle alternative atoms
      if ( *alt_loc !=' ')        // alternative
      {  
         int i;
         for (i=fisrtAtomIndexOfRes; i<atmIndex ; i++ )
         {
           if ( strncmp(atomList[i].atomName, pdbNameTrim(atom_name,4), 4)==0)
             break;            // found atom already exists
         }
         if ( i<atmIndex) {
           atmIndex--;          // discard the current atom (alternative)
           continue;                   
         }
      } 
    
      // new chains
      if ( chainid != aline[21] )	      // insert at least one chain
      {
         chainid = aline[21];
         ChainAdd(chainid, pdbResNum, resIndex, pdbAtomNum, atmIndex, modelid);
      }

      if ( atmIndex >= ATOM_MAXNUM-1 || resIndex >= MONO_MAXNUM-1 )
          break;                                           // exceeds maximum

      // now clean resIndex, atomIndex, assign Atom structure
      // 
      strcpy(myatom->atomName, pdbNameTrim(atom_name, 4));    // new atom always
      myatom->rn = resIndex;
      myatom->rnPDB= pdbResNum;
      sscanf(xyz, "%lf%lf%lf", &(myatom->x), &(myatom->y), &(myatom->z) );
      myatom->chainIndex = numChains-1;
      if ( isdigit(*charge) ) {
        myatom->charge = *charge-'0';
        if (*(charge+1) == '-' )
          myatom->charge = - myatom->charge ;
      }
      // estimate the radius: average PSGM/2.0 in all.prop
      myatom->radius = AtomRadius(myatom->atomName);
      //printf("%d(%d) %d(%d) %s\n", myatom->rn, myatom->rnPDB, atmIndex, pdbAtomNum, myatom->atomName);
    }
    else if ( !strncmp(aline, "MODEL", 5) ) 
    {
        modelid = atoi(pdbNameTrim(aline+10, 4));
    } 
    else if ( !strncmp(aline, "TER", 3)  ||  !strncmp(aline, "ENDMDL", 6)) 
    {
        lastPdbResNum = -999; 		  // to start a new model
        lastInsert=' ';
        chainid='\0';
    }
    else if ( !strncmp(aline, "HEADER", 6) )
    {
        strcpy(moleculeName, pdbNameTrim(aline+62, 4));
    }
    else if ( !strncmp(aline, "CONECT", 6) )
    {
        break;
    }
  }  // break on CONNECT
  numAtom=atmIndex+1;                  

  CalcDistNeighbors();

  // continue process CONECT section in PDB file
  int  nargs, nei1, nei2, nei3, nei4;
  do 
  {
     if ( !strncmp(aline, "CONECT", 6) )
     {
       aline[31]='\0';		// truncate hydrogen bond
       nargs = sscanf ( aline + 6, "%d%d%d%d%d", &pdbAtomNum, &nei1, &nei2, &nei3, &nei4 );

       switch ( nargs )		// convert to atom index
       {
         case 5: nei4 = atomIndex(nei4,anPDB);
         case 4: nei3 = atomIndex(nei3,anPDB);
         case 3: nei2 = atomIndex(nei2,anPDB);
         case 2: nei1 = atomIndex(nei1,anPDB);
                 pdbAtomNum = atomIndex(pdbAtomNum,anPDB);
       }

       if (pdbAtomNum>=0)
         switch ( nargs )		// connect neighbors
         {
           case 5: if (nei4>=0) Connect(pdbAtomNum, nei4, 1);	// check for existing bond first
           case 4: if (nei3>=0) Connect(pdbAtomNum, nei3, 1);
           case 3: if (nei2>=0) Connect(pdbAtomNum, nei2, 1);
           case 2: if (nei1>=0) Connect(pdbAtomNum, nei1, 1);
         } 
     }
     else
       break;
  }  while (fgets(aline,LINE_MAXLEN,crdfp)!= NULL );
  fclose(crdfp);   // unlock the file
  crdfp=NULL;
}  // ReadPdb()

void AtomArray::DisplayNextPdbModel()    // not for overlapping
{
  int maxModelId=-1, maxChainid=-1, j;
  Atom *atptr=atomList;

  if (InParm.filetag != fPDB) 
     return;

  if (chains != NULL )    // all chain displayed
  {
      maxChainid = numChains-1;              // first chain in last Model
      maxModelId =(chains+maxChainid)->model;   
      for (j=maxChainid; j>0 && (maxModelId==(chains+j-1)->model); j--);
      maxChainid=j;                 
  }
            
  if ( maxChainid<0 || displayChainid>=maxChainid )             
  {                                          // display all models after individual display
     for (j=0; j<numAtom ; j++, atptr++ )
       atptr->skip = 0;
     if ( maxModelId <= 1 )
       MsgBrd.set("Single Model");
     else
       MsgBrd.set("All Models");
     displayChainid=-1;               // mark the current chain
  }
  else		 // displayChainid can be <0       // display single model
  {
    int  displayModelid=-1, modelid=-1, newModelId=-1;

    if (displayChainid>=0)
       displayModelid=(chains+displayChainid)->model;
    for (j=0; j<numAtom; j++ , atptr++)
    {
       modelid = (chains+atptr->chainIndex)->model; 
       if (newModelId<0 && modelid>displayModelid )
       {
         newModelId = modelid;                         // next model
         displayChainid=atptr->chainIndex;      // mark the current chain
       }
       if (newModelId==modelid)
         atptr->skip = 0;
       else
         atptr->skip = 1;
    } 
    char textbuff[20]= "Model -";
    if ( newModelId >0 ) 
      sprintf(textbuff+6, "%d", newModelId);
    MsgBrd.set(textbuff);
  }
}

void AtomArray::DisplayNextPdbChain()    // not for overlapping, pdb only
{
   if (InParm.filetag != fPDB) 
     return;

   // display single chain in PDB file
   if ( ++displayChainid==numChains )
      displayChainid=0;    // back to first chain

   Atom *atptr=atomList;
   for (int j=0; j<numAtom ; j++, atptr++ )
   {
      if (chains==NULL || atptr->chainIndex==displayChainid )
        atptr->skip = 0;
      else
        atptr->skip = 1;
   } 

   char textbuff[200] = "chain -";
   if ( chains!=NULL && chains->model <= 0 )
     sprintf(textbuff, "chain %c", (chains+displayChainid)->id );
   else
     sprintf(textbuff, "model %d chain %c", 
            (chains+displayChainid)->model,(chains+displayChainid)->id);
   MsgBrd.set(textbuff);
}

// compute PDB neighbors without using radius info.
void AtomArray::CalcDistNeighbors( )
{
  // 1. connect within a monomer
  int currentrn=-1, firsti=0, lasti=-1, i;
  int iPrevC=lasti ;
  Atom *myatom, *jatom, *katom;
  double dd=0;

  for (i=0; i<numAtom; i++ )
  {
     myatom = atomList+i;
     myatom->numNeighbors=0;

     if (currentrn != myatom->rn )    // another residue    
     {
       firsti = i;
       currentrn = myatom->rn;       // don't compare for the first atom of a residue
     }
     else
     {
       for (int j=firsti; j<i; j++ )  // compare [firsti, i-1] atoms with the i atom
       {
         jatom = atomList+j;
         dd = Distance(myatom, jatom);    // atom distance
         if (dd <= BOND_MAXLEN && dd>=BOND_MINLEN )
           Connect(i, j);
       }
     }

     // 2. connect monomers
     if (strcmp("C", myatom->atomName) == 0 )
       iPrevC = i;				// remember C for binding to next monomer
     else if ( strcmp("N", myatom->atomName) == 0 )
     {
       if ( iPrevC >=0 )
       { 
         dd = Distance(myatom, atomList+iPrevC ); 
         if ( dd <= BOND_MAXLEN && dd >= BOND_MINLEN )    /* <=1.4 */
           Connect(i, iPrevC ); 
       }
     }
  }

  // delete H loops e.g. (H,H) distance can be 1.648 < BOND_MAXLEN
  //
  for (i = 0; i<numAtom ; i++ )
  {
    myatom = atomList+i;
    if (strlen(myatom->atomName) == 0 )
      myatom->skip = 1;						// skip non-existing atoms

    else if (  ( myatom->atomName[0] == 'N' || myatom->atomName[0] == 'C') &&
                 myatom->numNeighbors > 1 )			// GLY has two Hs on Ca	
    {
      // detect and eliminate possible H-H ring for C
      // (donot use "H" name)
      for (int j=0; j< myatom->numNeighbors; j++)		// for all neighbors of C
      {
        jatom = atomList + myatom->neighbors[j];		       // a neighbor of C has 2 or more neighbors
        if (jatom->numNeighbors == 1 )
          continue;
           		     
        for ( int ii=j+1; ii < myatom->numNeighbors; ii++ )     // for other C's neighbors
        {
          for ( int jj=0; jj<jatom->numNeighbors; jj++ )   	 // for C's neighbor's neighbors
          {
            if ( jatom->neighbors[jj]== myatom->neighbors[ii] ) // H loop detected
            {
              jatom->numNeighbors--;				// delete the bond from jatom(H)
              if ( jj < jatom->numNeighbors )		// rearrange the atom array
                jatom->neighbors[jj] = jatom->neighbors[jatom->numNeighbors];
					 					// delete the bond from another C's neighbor
              katom = atomList + myatom->neighbors[ii];	// the neighbor(should be H) with the loopped bond
              for ( int l=0; l<katom->numNeighbors; l++)
              {
                if ( katom->neighbors[l] == myatom->neighbors[j] )
                {
                  katom->numNeighbors--;
                  if ( l < katom->numNeighbors)
                    katom->neighbors[l] = katom->neighbors[katom->numNeighbors];
                    break;
                 }
               }
               jj--;							// undo jj++ for new neighbors
             }
           }
         }
       } // for (j)
     }
  }

  /* check neighoring info */
  for (i = 0; i<numAtom ; i++ )
  {
    if ( (atomList+i)->numNeighbors > ATOM_MAX_NEIGHBORS )
      fprintf(stderr, "PDB format error at atom#%d, %.4s #neighbors=%d \n", 
         i, (atomList+i)->atomName, (atomList+i)->numNeighbors );
  }
}

// connect two atoms, default to nocheck, i.e. checkfirst=0
void AtomArray::Connect(int atom1, int atom2, int checkfirst )        
{
  if ( checkfirst )
  {
    int i;

    for (i=0; i < (atomList+atom1)->numNeighbors && checkfirst; i++ )
       if ( (atomList+atom1)->neighbors[i] == atom2 )
         checkfirst=0;			// used as a flag for found

    if (checkfirst) 
      for (i=0; i < (atomList+atom2)->numNeighbors && checkfirst; i++ )
        if ( (atomList+atom2)->neighbors[i] == atom1 )
          checkfirst = 0;

    if (checkfirst == 0 )			// already connected
    {
       return;
    }

    // no max numNeighbors checking now
  }

  (atomList+atom1)->neighbors[(atomList+atom1)->numNeighbors] = atom2;
  ((atomList+atom1)->numNeighbors)++;
  (atomList+atom2)->neighbors[(atomList+atom2)->numNeighbors] = atom1;
  ((atomList+atom2)->numNeighbors)++;
}

int AtomArray::atomIndex(int pdbAtomNum, int *anPDB)   // binary search to find the internal atomIndex 
{
  int minIndex=0;
  int maxIndex=(numAtom>0)?numAtom-1:0;
  int midIndex=minIndex;
  while ( minIndex != maxIndex && maxIndex<numAtom)
  {
    midIndex=(minIndex+maxIndex)/2;
    if ( pdbAtomNum<anPDB[midIndex] ) 
      maxIndex=midIndex-1;
    else if (pdbAtomNum>anPDB[midIndex] ) 
      minIndex=midIndex+1;
    else
      break;
  }
  if (anPDB[midIndex]==pdbAtomNum)
    return midIndex ; 
  else
    return -1;
} 

// read coordinates without read in wcon, then connected them one after another
// 
//const unsigned int ATOM_MAXNUM  	=   100000;
void AtomArray::ReadXYZ(void) 
{
  int  counter=0;
  char aline[LINE_MAXLEN];
  Atom *myatom;
  int  fields;

  // reset coordinates when reach currentmol=0 
  // since in Crd file only one structure exists, the file will read only once
  if ( crdfp != NULL)
    return;

  currentmol = 1;         // only one structure
  if ((crdfp = fopen(InParm.crdName, "r")) == NULL) 
  {
    fprintf(stderr, "Unable to open file file : %s\n"); 
    exit(-1);
  }

  initAllocNames(1, ATOM_MAXNUM);      // put all crd into 1 residue
  strcpy((char*)resName, "NONE");      // any name

  // read in coordinates
  counter=0;
  while (fgets(aline,LINE_MAXLEN,crdfp)!= NULL && counter < numAtom )
  {
    myatom = atomList + counter;
    if ( aline[0] != '*' ) 
    {	
	fields =  sscanf(aline, "%lf%lf%lf", &(myatom->x), &(myatom->y), &(myatom->z) );
      if ( fields > 1 )
      {
        if ( fields == 2 )
          myatom->z = 0.0;
        myatom->radius = 0.05; 
        if ( fields >= 2 ) 
        {
          myatom->rn = 0;
          strcpy(myatom->atomName, "ALL-");	// an psudo atom name for build neighbor info
          if (InParm.dispmode == SPACEBALL ) 
            myatom->dispmode=SPACEBALL;
          else
            myatom->dispmode=STICK;
          counter++;
        }
      }
    }
  }  
  numAtom = counter;                  

  // connect each other as neighbor
  // setBackbone("ALL-");
  atomList[0].numNeighbors=0;
  for (int m=1; m<numAtom; m++)
  {
    atomList[m].numNeighbors =0;
    Connect(m, m-1);
  }
}  //ReadXYZ()


// rescale the coordinates for XYZ format file 
void AtomArray::XYZMaxScale() 
{
  double min[3], max[3];       // for (x, y, z)
  int i;
  Atom *myatom;

  myatom=atomList;
  min[0]=max[0]=myatom->x;
  min[1]=max[1]=myatom->y;
  min[2]=max[2]=myatom->z;
  for (i=0; i<numAtom ; i++, myatom++ )
  {
     if ( myatom->x < min[0] ) min[0] = myatom->x;
     if ( myatom->x > max[0] ) max[0] = myatom->x;

     if ( myatom->y < min[1] ) min[1] = myatom->y;
     if ( myatom->y > max[1] ) max[1] = myatom->y;

     if ( myatom->z < min[2] ) min[2] = myatom->z;
     if ( myatom->z > max[2] ) max[2] = myatom->z;
  }

  for (i=0; i< 3; i++)
  {
     if (max[i] != min[i] )   
       max[i] -= min[i] ;		/* delta d */
     else 
       max[i]=1.0;
  }

  for (i=0, myatom=atomList; i<numAtom ; i++, myatom++ )
  {
     myatom->x = 10.*(myatom->x-min[0])/max[0];
     myatom->y = 10.*(myatom->y-min[1])/max[1];
     myatom->z = 10.*(myatom->z-min[2])/max[2];
  }
}

void AtomArray::ReadCoordinate( void ) 
{
  switch(InParm.filetag) 
  {
  case fPTH : 
    if ( InParm.structend <= 0 || currentmol < InParm.structend )
      ReadPth();    			// read path file 
    break;
  case fDCD: 
    if ( InParm.structend <= 0 || currentmol < InParm.structend )
      ReadDcd();
    break;
  case fCRD:
    ReadCrd();				// read crd file
    break;
  case fPDB: 
    ReadPdb();				// read pdb file
    break;
  case fXYZ: 
    ReadXYZ();				// read xyz file
    break;
  default:
    fprintf(stderr, "Illegal coordinate file type\n");
    exit(-1);
  }

  // assign orgX,Y,Z for display purpose
  Atom *myatom=atomList;
  for ( int i=0;  i<numAtom; i++, myatom++ )
  {
    myatom->orgX=myatom->x;
    myatom->orgY=myatom->y;
    myatom->orgZ=myatom->z;
  }
  orgCrdModified=false;

  if ( InParm.filetag == fXYZ ) 
    XYZMaxScale();			// rescale for XYZ file format
  else 
    setSecStruct();                // set secondary structure flag

  SubtractCoMass();
}

void AtomArray::Print()
{

#if 0
  cout<<"Version Number: "<<versionNumber<<endl;
  cout<<"Molecule Name: "<<moleculeName<<endl;
  cout<<"Number of Residues: "<<numRes<<endl;
  int i, j, counter;
  for (i=0; i<numRes; i++)
    {
      cout<<"Residue Name: "<<resName[i]<<endl;
    }
  for (j=0; j<numAtom; j++)
    {
      cout<<"Atom #"<<j<<endl;
      cout<<"Radius: "<<atomList[j].radius<<endl;
      cout<<"Residue Number: "<<atomList[j].rn<<endl;
      cout<<"Atom Name: "<<atomList[j].atomName<<endl;
      cout<<"x: "<<atomList[j].x<<endl;
      cout<<"y: "<<atomList[j].y<<endl;
      cout<<"z: "<<atomList[j].z<<endl;
      cout<<"Neighbors: ";
      for (counter=0; counter<atomList[j].numNeighbors; counter++)
	{
	  cout<<atomList[j].neighbors[counter]<<"  ";
	}
      cout<<endl;
    }
#endif

}

void AtomArray::ReadConn()
{
  int counter;
  int i;
  char c[LINE_MAXLEN];
  FILE *fp;
	
  fp = fopen(InParm.connName, "r");
  if ( fp == NULL )
    return;

  fgets(c, LINE_MAXLEN, fp);
  sscanf(c, "%lg",&versionNumber);  // 9.0
  fgets(c, LINE_MAXLEN, fp);        // ~ CONNECTIVITY....
  fgets(c, LINE_MAXLEN, fp);        // ~ totmon....
  fgets(c, LINE_MAXLEN, fp);        //   data

  if ( sscanf(c, "%d%d", &numRes, &numAtom) <0  || numRes <=0  ||  numAtom <=0  )
  {
    fprintf(stderr, "Incorrect connectivity file format!\n");
    exit(-1);
  }
  initAllocNames(numRes,numAtom); 

  fgets(c, LINE_MAXLEN, fp);         // <molecular name>
  sscanf(c, "%s", &moleculeName);
  
  fgets(c, LINE_MAXLEN, fp);         // ~ Pointers...
  fgets(c, LINE_MAXLEN, fp);         //    data
  fgets(c, LINE_MAXLEN, fp);         // ~ Monomer names
  for(counter=0; counter<numRes; counter++)
    fscanf(fp, "%s", &resName[counter]);    // multi-line data

  fgets(c, LINE_MAXLEN, fp);         // skip newline of currrent line 
  fgets(c, LINE_MAXLEN, fp);         // ~ Pointers.....
  
  //int dummy(0);
  //double dummy2(0.0);
  double charge, tempRad1(0.0), tempRad2(0.0);
  double rat(0.0);

  for(counter=0; counter<numRes; counter++)
      fscanf(fp, "%*d");             // &dummy);     //     data

  fgets(c, LINE_MAXLEN, fp);        // skip newline of current line
  fgets(c, LINE_MAXLEN, fp);        // Properties.....
  fgets(c, LINE_MAXLEN, fp);        // ~ pt mono....
  for(i=0; i<numAtom; i++)
  {
    fgets(c, LINE_MAXLEN, fp);     // for supporting different wcon format    
    sscanf(c, "%*d%d%*d%*d%s%*lg%lg%lg%lg", &atomList[i].rn, &atomList[i].atomName, 
	                                                   &charge, &tempRad1, &tempRad2 );
    atomList[i].rn--;        // change order to index
    if ( fabs(tempRad1) < 1e-16 )
    {
      atomList[i].radius=0.5;      // for those tempRad1=0, set the radius to 0.5
    }
    else 
    {
      rat=tempRad2/tempRad1;              
      atomList[i].radius=pow(rat,(1.0/6));
    }
    atomList[i].charge = (float) charge;
  }

  // read in neighbor info 
  int first(0), second(0);
  fgets(c,LINE_MAXLEN,fp);
  sscanf(c,"%d %d",&first,&second);
  while (!strstr(c,"Angle") && !feof(fp))
  {
    if (sscanf(c,"%d %d",&first,&second) == 2 && first>0 && second>0)
    {
      first--;
      second--;
      Connect(first, second);
    }
    fgets(c,LINE_MAXLEN,fp);
  } 
  fclose(fp);
}

void AtomArray::ReReadConn()
{
  int i;
  char c[LINE_MAXLEN];
  FILE *fp;
  int first(0), second(0);
  int inBonds=0;

  // init	
  for(i=0; i<numAtom; i++)
  {
    atomList[i].numNeighbors = 0;
    atomList[i].skip=0;
    if (DispMode()== SPACEBALL ) 
      atomList[i].dispmode=SPACEBALL;
    else
      atomList[i].dispmode=STICK;
  }

  fp = fopen(InParm.connName, "r");
  if ( fp == NULL )
    return;

  fgets(c,LINE_MAXLEN,fp);
  while (!strstr(c,"Angle") && !feof(fp))
  {
    switch (inBonds )
    {
	case 0 : 
	  if (strstr(c, "Bonds") != NULL )
		inBonds = 1; 
	  break;
	case 1:
	  inBonds++;
	  break;
	case 2:
        if ( sscanf(c,"%d %d",&first,&second) == 2 && first>0 && second>0)
	  {
          first--;
          second--;
          Connect(first, second);
        }
        break;
    }

    fgets(c,LINE_MAXLEN,fp);
  } 
  fclose(fp);
}

void AtomArray::CalcNeighbors()
{
  int i, j, k(0);
  double ix, iy, iz, jx, jy, jz;
  for (i=0; i<numAtom; i++)
    {
      for (k=j=0; j<numAtom; j++)
	{
	  ix=atomList[i].x; iy=atomList[i].y; iz=atomList[i].z;
	  jx=atomList[j].x; jy=atomList[j].y; jz=atomList[j].z;
	  if ((((ix-jx)*(ix-jx)+(iy-jy)*(iy-jy)+(iz-jz)*(iz-jz)) <
	       ((atomList[i].radius+atomList[j].radius)*(atomList[i].radius+atomList[j].radius))/4.0)&&(i!=j))
	    {
	      atomList[i].neighbors[k]=j;
	      k++;
	    }
	  atomList[i].numNeighbors=k;
	}
    }
}

// check the cutoff distance InParm.pickradius to set skip
void AtomArray::setSphereSkip()
{
  if ( InParm.picksphere ==  NULL)
    return;

  Atom *atom, *atomB, *firstAtom;
  int m, mb, lastRn=-1, inSphere=0;
  double dmin = InParm.pickradius;     // cutoff distance
  for (m=0,atom=atomList; m<numAtom; m++, atom++)
  {
    if (lastRn != atom->rn)
    {
      // check last residur for any atom in the selected sphere
      if (lastRn>=0 && inSphere==0)
      {
        while ( firstAtom != atom )
        {
          if (firstAtom->skip == 0 )
            firstAtom->skip=1;
          firstAtom++;
        }
      }
      lastRn=atom->rn;
      inSphere=0;      // any atom is selected in a residue, then the residue is selected
      firstAtom=atom;
    }
    else if (inSphere==1)
      continue;
         
    if ( (atom->dispmode&INSPHERE)==INSPHERE)  // not in selection
      inSphere=1;
    else
    {
       for (mb=0, atomB=atomList; mb<numAtom; mb++, atomB++)
       {
           if ( (atomB->dispmode&INSPHERE)==INSPHERE )  // loop over selections
           {
             if ( Distance(atom, atomB) <= dmin )
             {
               inSphere=1;  // the atom is not cutoff
               break;
             }
           }
       }
    }
  }
  if (lastRn>=0 && inSphere==0)       // check last residue
  {
    while ( firstAtom != atom )
    {
      if (firstAtom->skip == 0 )
        firstAtom->skip=1;
      firstAtom++;
    }
  }
}

void AtomArray::setSphere()
{
  if ( InParm.picksphere ==  NULL)
    return;

  Atom *atom, *atomB;
  int  m;

  double x=0,y=0,z=0, dd=0, dmin=-1;
  int n=0, i=0;

  for (m=0, atom=atomList; m<numAtom; m++, atom++)
  {
    if ( inPick(m, InParm.picksphere) )
    {
      x +=  atom->x;
      y +=  atom->y;
      z +=  atom->z;
      atom->dispmode |= INSPHERE;
      n++;
    } else if ( (atom->dispmode&INSPHERE)==INSPHERE ) 
      atom->dispmode ^= INSPHERE;    // reset bit
  }

  if ( n > 0 )   // get selection center
  {
    x /= n;         // average the sum, == center
    y /= n;
    z /= n;
    for (m=0,i=0,atom=atomList; m<numAtom; m++, atom++)
    {
      if ( (atom->dispmode&INSPHERE)==INSPHERE )
      {
          atomB = atomList+i;
          dd=(atom->x-atomB->x)*(atom->x-atomB->x)+(atom->y-atomB->y)*(atom->y-atomB->y)+(atom->z-atomB->z)*(atom->z-atomB->z);
          if ( dmin<0 || dd<dmin )
          {
             dmin = dd;
             i=m;
          }
      }
    }
    pickSphereCenter=i; // get center index
    //atom = atomList+i;         // get center coordinate
    //SubtractCoMass(atom->x, atom->y, atom->z);    // move to sphere center
  }
  setSphereSkip();

} // setSphere()


void AtomArray::setSkip()//////SKIPS THE SPECIFIED RES'S
{
  if (InParm.skip != NULL )
  {
    for (int m=0; m<numAtom; m++)
    {
      if ( inPick(m, InParm.skip) )
      {
        atomList[m].skip=1;
      }
    }
  }
}

void AtomArray::setSkipSurf()//////SKIPS THE SPECIFIED res on computing surface
{
  if (InParm.skipsurf != NULL )
  {
    for (int m=0; m<numAtom; m++)
    {
      if ( inPick(m, InParm.skipsurf) )
      {
        atomList[m].skipsurf=1;
      }
    }
  }
}

void AtomArray::setStick()//SETS CERTAIN ATOMS TO STICK MODE
{
  if (InParm.stick != NULL )
  {
    for (int m=0; m<numAtom; m++)
    {
      if ( atomList[m].skip==0 && inPick(m, InParm.stick) )
      {
        atomList[m].dispmode |= STICK;
        if ( (dispmodeComb&STICK) == 0 ) 
          dispmodeComb |= STICK;
      }
    }
  }
}

void AtomArray::setSpace()	//SETS CERTAIN ATOMS TO SPACE MODE
{
  if (InParm.space != NULL )
  {
    for (int m=0; m<numAtom; m++)
    {
      if (atomList[m].skip==0 &&  inPick(m, InParm.space) )
      {
        atomList[m].dispmode |= SPACEBALL;
        if ( (dispmodeComb&SPACEBALL) == 0 ) 
          dispmodeComb |= SPACEBALL;
      }
    }
  }
}

// set all protein to be ribbon mode if global displayMode is ribbon
void AtomArray::setRibbon()  //SETS ONLY THE 20 PROTEIN RESIDUES TO BE IN RIBBON MODE
{
  if (InParm.ribbon != NULL )
  {
    for (int m=0; m<numAtom; m++)
    {
      if (atomList[m].skip==0 && /* isRibbonProtein(resName[atomList[m].rn]) && */  inPick(m, InParm.space)!=0 )
      {
        atomList[m].dispmode |= RIBBON;
        if ( (dispmodeComb&RIBBON) == 0 ) 
          dispmodeComb |= RIBBON;
      }
    }
  }
}

int AtomArray::isRibbonProtein(char *pname)
{
  return AminoAcidsObj.IsAA(pname);
}

// set displayMode upon global dispMode inParm
void AtomArray::setDefaultDispMode()
{
  int i;
  Atom *listi;

  dispmodeComb=InParm.dispmode;

  for (i=0,listi=atomList; i<numAtom; i++, listi++)
  {
    if ( listi->skip==0)
	listi->dispmode=dispmodeComb;
  }
}

// Set Secondary Structure type(.secStruct) upon torsional angle for display purpose :
//   STRUCT_HELIX, STRUCT_SHEET, STRUCT_UNKNOWN
// Also set the flag STRUCT_CHAINEND.   
//
// ATOM(i).secStruct indicates the connection type between ATOM(i) and ATOM(i+1).
//
// phi=[-180, 0] && psi=[-170,0]: helix,  three continuoues H-bonds [i, i+4] dist(N,O)<3.5A
// phi=[-180, 0] && psi=[0, 190]: sheet
//                        others: unknow
//         
//  pass 2 is the current pass to scan atom list.  pass 1 stores previous pass 2
//
const int    MinNumHelixRes=5;     // min number of (phi-psi) pair residues in a helix =4
const int    MinNumSheetRes=4;     // min number of residues in a sheet =3
void AtomArray::setSecStruct()
{
  //return;

  int           i=0,      numResRead=0;
  Atom          *listi=atomList, *firstAtom=atomList, *currRes=atomList, *oldRes=atomList;
  SecStructType currType, newType;

  if ( atomList == NULL )  
    return;

  setSecStructStep1();

  currType=atomList->secStruct;
  newType=currType; 
  for ( ; i<=numAtom; i++, listi++)
  {
     // check for a new residue, no read in atom info yet
    if (i>0 && currRes->rn == listi->rn && i<numAtom)
      continue;

    numResRead++;                  // #residues read in for current structure 
    if ( i < numAtom)
      newType=listi->secStruct;      // next secStructType, current residue is currRes
    else 
      newType=STRUCT_UNKNOWN;

    // get a new residue
    if ( sIsUnknown(newType) && (newType&0x0F) )
    {
      newType &= 0xF0;
      if ( i<numAtom)
        listi->secStruct &= 0xF0;
    }
    else if (numResRead==2 )
      firstAtom->secStruct |= (newType&0x0F);

    if (sLowType(newType)!=sLowType(currType) || 
        sIsChainEnd(currRes->secStruct) || i==numAtom)  
                // Type changed, set secStruct now
    {
      if ( sIsHelix(currType) && numResRead<MinNumHelixRes ||
           sIsSheet(currType) && numResRead<MinNumSheetRes   )  // helix must have more than 5 residues
        currType = STRUCT_UNKNOWN;   // STRUCT_UNKNOWN;                    
      else
        currType = sLowType(currType);
      
      if ( currType != sLowType(newType) || sIsChainEnd(currRes->secStruct) )
      {
        while (firstAtom!=listi)  
        {
          if ( sIsSheet(currType) && 
               ( (sIsChainEnd(currRes->secStruct) && firstAtom==oldRes) ||  firstAtom==currRes ) )
            firstAtom->secStruct |= STRUCT_END;
          else
            firstAtom->secStruct = ( firstAtom->secStruct&0xF0|currType);
          firstAtom++;
        }  
        currType=newType;
        numResRead=0;
      }
    } 
    oldRes=currRes;
    currRes=listi;
  } // for(i)  

  for (i=0,listi=atomList; i<numAtom; i++, listi++)
  {
    if ( !sFoundO(listi->secStruct) &&  !sIsUnknown(listi->secStruct) )
      listi->secStruct &= 0xF0;              //unable to determine the peptide plane
  }      

  // printf("set secondary structure\n");

#if 0
  for (i=0,listi=atomList; i<numAtom; i++, listi++)
    if (i==0 || (listi-1)->rn != listi->rn) 
      printf("rn=%d\tsecStruct=%d, chainEnd=%c\n",  
                   listi->rn, 
                   listi->secStruct&0xF, 
                   ((listi->secStruct&STRUCT_CHAINEND)!=0)?'1':' ');        
#endif

}  // setSecStruct()

// get all surfaces for multi-structures before setting the 
                                   
void AtomArray::getSurfaces()
{
  static bool msg=false;
  int cmol = currentmol;

  setSkipSurf();                          // check if an atom should be ignore in surf.exe

  if (NUM_AA_TYPE > 125 )
  {
    fprintf(stderr,  "Error: Protein color type occupy more than one byte!\n");
    return;
  }
  if ( InParm.surfName[0] == '\0' )      // only read in once
  {
    strncpy(InParm.surfName, InParm.crdName, FILENAME_MAXLEN);
    if ( strlen(InParm.surfName) >= FILENAME_MAXLEN-5 )
      InParm.surfName[FILENAME_MAXLEN-5] = '\0';
    strcat(InParm.surfName, ".vtx");
    if (    (InParm.useExistingVtxFile==false) 
            ||  ((vertexfp=fopen(InParm.surfName, "rb"))==NULL) )   // file not exist
    {
#ifdef _WIN32
      vertexfp=fopen(InParm.surfName, "wb");     // merged output vertex files of surf.exe in binary
#else
      vertexfp=fopen(InParm.surfName, "w");     // merged output vertex files of surf.exe in binary
#endif
      if ( vertexfp!= NULL )
      {
        fwrite(&maxNumSurfTris, sizeof(int), 1, vertexfp);      // leave first 4 bytes for maxNumSurfTris
        currentmol=0;
        for (int i=0; i<=InParm.structno; i++)
        {
          ReadCoordinate();  // only path and dcd support movie
          getSurface(i);     // run surf.exe for each structure, assign maxNumSurfTris
        }
        rewind(vertexfp);               // move to beginning of file
        fwrite(&maxNumSurfTris, sizeof(int), 1, vertexfp);   // write maxNumSurfTris

        fclose(vertexfp);

        // cleanup intermediate files
        UNLINK("surf_in_tmp");  
        UNLINK("surf_in_tmp.log");
        UNLINK("surf_in_tmp.tri");
      }
    } 
    else
    { 
      // read existing vertex file to get maxNumSurfTris
      maxNumSurfTris=0;
      fread(&maxNumSurfTris, sizeof(int), 1, vertexfp);        // read first integer to get maxNumSurfTris 
      fclose(vertexfp);

      if ( maxNumSurfTris <= 0 )
        printf("Corrupted surface vertex file: %s\n", InParm.surfName);
    }
    printf("max number of surface vertex=%d\n", maxNumSurfTris);

    // allocate max memory for vertex crds
    if ( surface == NULL && maxNumSurfTris>0 )      
      surface=new SurfTriangle[maxNumSurfTris+1];   // add 1 triangle at end end for loop terminating
    vertexfp=fopen(InParm.surfName, "rb+");         // open vertex file for reading, changed from rb to rb+
  }
  currentmol = cmol;
  readSurface();                         // read surface vertexes into memory for current structure
  if ( !msg )
  {
    printf("\nCMOIL: Surfaces generated by A. Varshney's SURF which was modified by R. Gillilan in XTM at ftp://bridgedec.chess.cornell.edu/pub/XTM-src.tar.gz .\n");
    msg=true;
  }
  getSubArea();

}

// get surface for single structure in aptr->list: run surf.exe and save to file *.vtx
void AtomArray::getSurface(int isurf)
{
  const int 	nTextLinePerTri=4;
  int     		i ;
  Atom   		*aptr;
  char   		 textLine[LINE_MAXLEN];
  char   		*infile="surf_in_tmp", *logfile="surf_in_tmp.log",  *outfile="surf_in_tmp.tri" ;
  SurfTriangle    sbuff;      

  printf("CMOIL: Generating surface #%d....\n", isurf+1);

  // compose the input for "surf.exe"
  //FILE *surfIn = fopen(infile, "wt");              // input data file for  
  FILE *surfIn = fopen(infile, "w");              // input data file for  
  if ( surfIn != NULL)
  {
    for (i=0, aptr=atomList; i<numAtom; i++, aptr++) 
    {
      if ( aptr->skipsurf == 0)
      {
        // surf.exe input format:  <atom_id> <atom_radius> <atom_center_x> <atom_center_y> <atom_center_z>
        fprintf(surfIn, "%d %.3f %.3f %.3f %.3f\n", i, aptr->radius, aptr->x, aptr->y, aptr->z);
      }
    }
    fclose(surfIn);
  } 
  else 
  {
    fprintf(stderr, "CMOIL: unable to open file, %s,  to write !\n",infile );
    return;
  }

  // execute surf.exe program
  sprintf(textLine, "%ssurf%s -R %f -E %f -W 1 %s 2>%s", InParm.exedir, BINSURFFIX, InParm.surfprobe, InParm.surfprobe, infile, logfile);
  printf ("surf cmd=%s\n", textLine);
 
  system(textLine);                     // ==surfCmd

  // get surf outputs into vertexOut
  int nSurfTris;
  FILE *surfLog= fopen(logfile, "rt");
  while ( fgets(textLine, LINE_MAXLEN, surfLog) != NULL ) 
  {
    if ( sscanf(textLine, "Total Triangles %d", &(nSurfTris)) >= 2 )      //numSurfaceTriagles
      break;
  }
  fclose(surfLog);

  surfHeader.numSurfTris=nSurfTris;      // init  
  surfHeader.outerProced=0; 
  surfHeader.cavityNormalReversed=0;

  fwrite(&surfHeader, SURF_HDR_LEN, 1, vertexfp);   // first elemnt of the block:  number of vertexes

  FILE *surfTri= fopen(outfile, "rt");
  int ii;

  sbuff.out=0;   // read from surf.exe output and write into binary vertex file, no outer surface processing yet
  sbuff.edge=0;
  i=0;
  bool notBuried=true;
  int  nOutputTris=0;    // actual tris written to vertex file
  int  edge=-1, lastAtom=-1, lastEdge=-1, outEdge=-1;
  while (i<nSurfTris*nTextLinePerTri && fgets(textLine, LINE_MAXLEN, surfTri) != NULL  ) 
  { 
      ii = (i+3)%nTextLinePerTri;       // map to i=[0,numOfTriLines] to  3, 0, 1, 2
      if (ii==nTextLinePerTri-1) 
      {
        sscanf(textLine, "%d%d%d", &(sbuff.atomNum), &edge, &(sbuff.nextAtom));         // get atom num 

        // Fix small H radius bury bug in surf.exe output. It causes H displays "mushroom" like isolated surfaces.
        // (must for cavity display)
        notBuried=true;
        if ( sbuff.atomNum >= 0 ) 
        {
          Atom *a=atomList+sbuff.atomNum;      
          if ( ( *(a->atomName)=='H' || *(a->atomName+1)=='H' ) && a->numNeighbors>0 )
          {
            Atom *b = atomList+a->neighbors[0]; 
            Coordinate c1(*a), c2(*b);
            if ( a->radius + (c1-c2).len() < b->radius )
              notBuried=false;
          }
          if ( notBuried )
          {
            if ( sbuff.nextAtom < 0 )
              sbuff.nextAtom = sbuff.atomNum;         // change -1 to link itself

            // reorder edge number 
            if ( sbuff.atomNum != lastAtom )
            {
              lastAtom=sbuff.atomNum;
              outEdge=0;
              lastEdge=edge;
            } 
            else if ( edge != lastEdge )
            {
              outEdge++;
              lastEdge=edge;
            }
          }
        }
      }
      else
        sscanf(textLine, "%f%f%f%f%f%f", sbuff.vt[ii],sbuff.vt[ii]+1, sbuff.vt[ii]+2, sbuff.nm[ii],sbuff.nm[ii]+1, sbuff.nm[ii]+2);
      if (ii==2 && notBuried && sbuff.atomNum >=0 )
      {
        nOutputTris++;
        sbuff.edge = sbuff.atomNum*MaxNumSurfEdgePerAtom+outEdge;    // change to global unique edge_id=atom#*25+face_id
        fwrite(&sbuff, SURF_TRI_LEN, 1, vertexfp);
      }
      i++;
  }
  fflush(vertexfp);

  if ( nOutputTris != nSurfTris )
  {
    surfHeader.numSurfTris=nOutputTris;
    fseek(vertexfp, (long)(-(SURF_TRI_LEN*nOutputTris+SURF_HDR_LEN)), SEEK_CUR);  
    fwrite(&surfHeader,SURF_HDR_LEN, 1, vertexfp); 
    fseek(vertexfp, (long)(SURF_TRI_LEN*nOutputTris), SEEK_CUR); 
    fflush(vertexfp);
  }
 
  if ( nOutputTris > maxNumSurfTris )
    maxNumSurfTris=nOutputTris;                 // maximum number of vertexes 

  fclose(surfTri);
  printf("CMOIL: ... surface done!\n\n");
}

// read next surface vertexes from file into memory
void AtomArray::readSurface()
{
  if (surface == NULL)
    return;
 
  if ( (currentsurf>=InParm.structno) ||  (currentsurf>=currentmol) ) 
  {
    currentsurf = 0;
    if (maxNumSurfTris>0 && surfIndex ) 
    {
      delete[] surfIndex;                                     // free resource
      surfIndex=NULL;
    }
  }

  if ( vertexfp!=NULL )
  {
    if (currentsurf==0) 
    {
      if ( maxNumSurfTris >  0 ) 
        fseek(vertexfp, (long)(sizeof(int)), SEEK_SET);     			// skip first integer of maxNumSurfTris
      else 
        fread(&maxNumSurfTris, sizeof(int), 1, vertexfp);   	// read first integer to get maxNumSurfTris 
    }

    currentsurf++; 
    while (currentsurf<currentmol)
    {     
      if ( fread(&surfHeader, sizeof(surfHeader), 1, vertexfp)==0 )
        surfHeader.numSurfTris=0 ;  // read header structure to get numSurfTris
      else
        fseek(vertexfp, (long)(surfHeader.numSurfTris*SURF_TRI_LEN), SEEK_CUR); // skip following block
      currentsurf++;
    }

    // read in surface for current structure
    if ( fread(&surfHeader, SURF_HDR_LEN, 1, vertexfp)==0)
      surfHeader.numSurfTris=0  ;       // read first integer to get numSurfTris
    else if (surface != NULL)    // read in "nmlXYZ, XYZ, color_type" block
    {
      fread(surface, (surfHeader.numSurfTris*SURF_TRI_LEN), 1,  vertexfp); 
      surface[surfHeader.numSurfTris].atomNum=-1;
      surface[surfHeader.numSurfTris].edge=-100;      // as terminate condition on edge loopping

      // get min/maxmum charge on surface
      int i, atomNum;
      SurfTriangle *sptr=surface;
      minAtomCharge=maxAtomCharge=atomList[sptr->atomNum].charge;
      for (i=1, sptr++; i<surfHeader.numSurfTris; i++, sptr++)
      {
        atomNum = sptr->atomNum;
        if ( atomList[atomNum].charge < minAtomCharge )
          minAtomCharge=atomList[atomNum].charge;
        if ( atomList[atomNum].charge > maxAtomCharge )
          maxAtomCharge=atomList[atomNum].charge;
      }
      markCavity();
    }
  }
}

//************ for cavity display *******************************
const int SurfEdgesEmpty=-1000000;     // no edges exist from this edge id on
const int SurfEdgeInvalid=-1;          // the edge deleted

// the  bottle neck of the performance for cavity finding :
// check all triangles in both edges for connecting information
// no restrict to atom#s, edge1 and edge2 can belong to different atoms
bool AtomArray::surfEdgeConnected(SurfEdgeIndex *edge1, SurfEdgeIndex *edge2)
{
  SurfTriangle *sptr1, *sptr2;
  int i, j, m, n;
  float *a, *b; 
  bool notConnected=true;  

  for (i=edge1->startTri, sptr1=surface+i; notConnected && i<edge1->startTri+edge1->triCnt; i++, sptr1++)
    for (j=edge2->startTri, sptr2=surface+j; notConnected && j<edge2->startTri+edge2->triCnt; j++, sptr2++)
      // surf.exe output triagles are connected by common vertexes
      for (m=0, a=(float*)sptr1->vt; notConnected && m<3; m++, a+=3 )
        for (n=0, b=(float*)sptr2->vt; notConnected && n<3; n++, b+=3 )
        {
          if ( *a==*b && *(a+Y)==*(b+Y) && *(a+Z)==*(b+Z) )      // same coordinate
             notConnected=false;
        }
  return (!notConnected);

} // surfEdgeConnected()

// connect the outerEdge to another atom if any connection exists
void AtomArray::surfAtomEdgeConnect(int outerEdge, int atomIndex)
{
  if ( atomIndex >= 0 )
  {
    int i,  ei0=atomIndex*MaxNumSurfEdgePerAtom;    	// first edge index of the atom
    SurfEdgeIndex *ii;   					// first tri of edge i
    for (i=0, ii=surfIndex+ei0; (ii->startTri!=SurfEdgesEmpty) && i<MaxNumSurfEdgePerAtom; i++, ii++)
      if (ii->startTri>=0)                  
      {
        if (surfEdgeConnected(surfIndex+outerEdge,ii) )
          buildSurfOutConnectivity(atomIndex, ei0+i);
      }
  }
}


// put all connected edges into same edge id for each atom
// modify surfIndex for the edge line number
// after all done then writeback to the surface file
//
void AtomArray::setSurfAtomEdgeSelfGroup(int *maxzEdge)
{
  for (int i=0; i<numAtom; i++)
  {  
    SurfEdgeIndex *ii=surfIndex+i*MaxNumSurfEdgePerAtom;             // first index of atom i
    if (  ii->startTri==SurfEdgesEmpty  )  
      continue;  

    bool changed=true; 
    int j, k;
    SurfEdgeIndex *jj, *kk;                             // point to edge j, k
    while (changed)
    {
      changed=false;
      for (j=0, jj=ii; (jj->startTri!=SurfEdgesEmpty) && j<MaxNumSurfEdgePerAtom; j++, jj++)
      {
        if ( jj->startTri<0 ) 
          continue;    
      
        for (k=j+1, kk=ii+k; (kk->startTri!=SurfEdgesEmpty) && k<MaxNumSurfEdgePerAtom; k++, kk++)
        {
          if (  kk->startTri < 0 )
            continue; 

          if ( surfEdgeConnected(jj, kk) )
          {
            if ( (surface+kk->startTri)->edge == *maxzEdge )    
              *maxzEdge = (surface+jj->startTri)->edge ; 	// update maxzEdge if it's merged 
            //printf("i=%d, j=%d, k=%d\n", i, j, k);
            reorderEdgePair(j, k, jj, kk);
            changed=true;
          }
        } 
      }
    }
  }
}  // setSurfAtomEdgeSelfGroup()

// reorder edges put same edge# block together
// j & k must belong to same atom i and j<k: [0,24], and
//       jj=surfIndex+25i+j, kk=surfIndex+25i+k
void AtomArray::reorderEdgePair(int j, int k, SurfEdgeIndex *jj, SurfEdgeIndex *kk) 			
{
  const int surfBuffSize=200;
  SurfTriangle sbuff[surfBuffSize];

  int sj=jj->startTri, sk=kk->startTri;      // tri start positions
  int cj=jj->triCnt,   ck=kk->triCnt;        // tri counts

  if ( sj<0 || sk<0)
    return;

  // change edge2# to edge1#
  int n, edge1=(surface+sj)->edge;
  SurfTriangle *nn;
  for (n=sk, nn=surface+sk; n<sk+ck; n++, nn++)
    nn->edge = edge1;    // set to same edge#  if two edges belong to the same atom

  // merge edge k into edge j in an atom
  // must:  jj < kk  (j<k),  jj=ii+j, kk=ii+k
  // ck:  total block size( number of triangles) to move 
  for (int blockLEN=ck; blockLEN>0 ; blockLEN -=surfBuffSize )          // loop in case #triangles>50
  {
    int blocklen = ((blockLEN>=surfBuffSize)?surfBuffSize:blockLEN);	// size(in tris) for each copy command
    int blockbytes=blocklen*sizeof(SurfTriangle);            		// size(in bytes) for each copy
    int insertIndex =sj+cj;               					// start of the block after j. 
    int targetIndex =sk+ck-blocklen;               				// beginnig of the segement to move in blcok k

    if ( insertIndex < targetIndex)						// don't do any thing if move to same position
    {
      memcpy(sbuff, surface+targetIndex, blockbytes);         	// copy k block to buff 
      for (int m=targetIndex; m>insertIndex; m -= blocklen)
      {
        if ( m<blocklen)
          memcpy(surface+blocklen, surface, m);              	// at very begining of surface
        else
          memcpy(surface+m, surface+m-blocklen, blockbytes);    	// move all blocks down from i
      }
      memcpy(surface+insertIndex, sbuff, blockbytes);
    }
  }

  //modify indexes
  jj->triCnt += ck;  			// kk moved into jj
  kk->startTri=SurfEdgeInvalid;  	// delete kk info

  SurfEdgeIndex *mm=jj+1;                   		// update all edges between jj-kk
  for (n=j+1; n<k; n++, mm++)
  {
    if ( mm->startTri >=0 )
      mm->startTri += ck;     	// update start position, keep lengh at mm->triCnt
  }
} // reorderEdgePair()

//  recursive to build outer surface connectivity
//  outerEdge, outerEdge :  known atom-edge on the outer surface, 
// then build connectivity to all triangles in this edge 
void AtomArray::buildSurfOutConnectivity(int outerAtom, int outerEdge)
{  
  if ( outerAtom<0 || outerEdge < 0 ) 
    return; 

  int firstLine =  surfIndex[outerEdge].startTri;
  if (firstLine < 0 )
    return;  

  SurfTriangle *sptr=surface+firstLine; 
  if (sptr->out != 0 )
    return;                 // whole edge already set to outer surface

  surfIndex[outerEdge].out=true;     // mark whole edge is on outer surface

  int iatom, lastAtom=-1;
  for (int i=0; i<surfIndex[outerEdge].triCnt; i++, sptr++)
  {
    sptr->out =1;		              // put self to outer surface first to stop processing outerAtom again from nested calls
    iatom = sptr->nextAtom;           // atom to connect
    if ( iatom != outerAtom && iatom != lastAtom )    // not connect to self, and check another atom once
    {
       lastAtom=iatom;
       surfAtomEdgeConnect(outerEdge, iatom);		// call buildSurfOutConnectivity(iatom, edge) from it
    } 
  }
} //buildSurfOutConnectivity()

//find all surface atoms (no cavity)
//surface array should not null,  set aptr->skipsurfdisp
void AtomArray::markCavity()
{
  if (surface==NULL || InParm.showCavity==false || surfHeader.outerProced != 0 || vertexfp == NULL )
    return;

  printf("\nCMOIL: Processing surface data to identify cavities....\n");

  numSurfEdges=numAtom*MaxNumSurfEdgePerAtom;
  // get all surf atom index, triagles are in atom# order
  if ( surfIndex == NULL )
  {
      surfIndex=new SurfEdgeIndex[numSurfEdges];      // for atom edge information, store offset for each atom/edge
                                                      // for start/stop line#
  }
  if ( surfIndex == NULL)
     return;
     
  int 	i,j;
  SurfEdgeIndex *ii;
  for ( i=0, ii=surfIndex; i<numSurfEdges; i++, ii++) // init
  { 
     ii->startTri = SurfEdgesEmpty ;                  // init to empty list
     ii->out=false;
  }

  int 	edgeI=-1, maxzEdge=-1;
  int       atomI=-1, oldAtom=-1, oldEdge=-1, maxzAtom=-1;
  SurfTriangle  *sptr=surface;
  float     maxZ,z;
  
  // find triangle with max(z) value, which must on the outer surface
  // -- several surface clusters are not supported
  maxZ = sptr->vt[0][Z];      
  for (i=0; i<surfHeader.numSurfTris; i++, sptr++)
  {
    if (edgeI != sptr->edge )   // a new edge
    {
       oldEdge=edgeI;        
       atomI = sptr->atomNum ;
       edgeI = sptr->edge;
       surfIndex[edgeI].startTri= i;                  // store index
       if (oldEdge>=0)
         surfIndex[oldEdge].triCnt=i-surfIndex[oldEdge].startTri;     // store length
    }
    for (j=0; j<3; j++)
    {
      z = sptr->vt[j][Z];             // find first outer surface triangle
      if ( z > maxZ )
      {
         maxZ=z;
         maxzAtom=atomI;
         maxzEdge=edgeI;
      }
    }
  }
  if (edgeI>=0)
     surfIndex[edgeI].triCnt=i-surfIndex[edgeI].startTri;

  setSurfAtomEdgeSelfGroup(&maxzEdge);
  buildSurfOutConnectivity(maxzAtom, maxzEdge);
  fixSurfConnectivity();

  // reverse cavity normal vectors, and write back to file
  for (i=0, sptr=surface; i<surfHeader.numSurfTris; i++, sptr++)
  {
    if ( sptr->out==0)
    {
       float *f = (float*)sptr->nm;
       for ( j=0; j<9; j++, f++)
         *f = - *f; 
    }
  }

  // save modified code back to vertex file
  surfHeader.outerProced=1;
  surfHeader.cavityNormalReversed=1;
  i = surfHeader.numSurfTris*SURF_TRI_LEN ;   // surface block length
  if ( fseek(vertexfp, -(long)(SURF_HDR_LEN+i),SEEK_CUR) != 0 )
  {
    fprintf(stderr, "file operation error!\n");
    return;
  }
  fwrite(&surfHeader,SURF_HDR_LEN,1,vertexfp);
  fwrite(surface,i,1,vertexfp);
  fflush(vertexfp);

  printf("\nCMOIL: .... cavities done!\n");

} // markCavity()

// connectivity recheck: fix surf.exe output problem
void AtomArray::fixSurfConnectivity()
{
  SurfEdgeIndex *jj, *mm;
  SurfTriangle  *sptr, *sptr2;
  int i, j, k, l;
  int theAtom;

  for (i=0, theAtom=0; i<numSurfEdges; i+=MaxNumSurfEdgePerAtom, theAtom++)     // init 
  {
    if ( surfIndex[i].startTri == SurfEdgesEmpty)
      continue;         // no surface atoms
   
    for ( j=i,jj=surfIndex+i; (jj->startTri!=SurfEdgesEmpty) && !(jj->out) && j<i+MaxNumSurfEdgePerAtom; j++, jj++);
    if ( (jj->startTri!=SurfEdgesEmpty) && j<i+MaxNumSurfEdgePerAtom )
      continue;         // found a surface outer edge for theAtom
    
    // no edge on outer surface, check if any linked atoms on outer surface
    for ( j=i, jj=surfIndex+i; (jj->startTri!=SurfEdgesEmpty) && j<i+MaxNumSurfEdgePerAtom; j++, jj++)   	// loop on edges of theAtom
    {
      if ( jj->startTri<0 )                   // edge not exist
        continue;   
      
      int linkedAtom=-1, lastLinked=-1;
      for (k=0, sptr=surface+jj->startTri; k<jj->triCnt; k++, sptr++)  // loop on tris of the edge
      {
        // Fix mismatched atom link pair bug in surf.exe output
        linkedAtom=sptr->nextAtom; 
        if (linkedAtom != lastLinked && linkedAtom>=0 && linkedAtom != sptr->atomNum )
        {
          lastLinked = linkedAtom;       // linked atom, check if this atom is also linked back to it
          int m = linkedAtom*MaxNumSurfEdgePerAtom; 
          for (l=0, mm=surfIndex+m; (mm->startTri!=SurfEdgesEmpty) && l<MaxNumSurfEdgePerAtom && mm->startTri<0 ; l++)
             mm++;
          if ( mm->startTri==SurfEdgesEmpty || l>=MaxNumSurfEdgePerAtom )
             continue;

          sptr2 = surface+mm->startTri;
          while ( sptr2->atomNum==linkedAtom && sptr2->nextAtom!=theAtom )
            sptr2++;
          if ( sptr2->atomNum==linkedAtom)
            continue;                 // theAtom is linked from the lastLinked, thus already checked

          // found unpaired atom link, fix it
          for (l=0, mm=surfIndex+m; (mm->startTri!=SurfEdgesEmpty) && l<MaxNumSurfEdgePerAtom; l++, mm++)  
          {
            if ( mm->out)
            {
              //printf("outedge=%d(atom#%d), check atom #%d\n", m+l, linkedAtom, theAtom);
              //if (theAtom != 1516 || m+1!= 21551) 
              surfAtomEdgeConnect(m+l, theAtom);
            }
	    }
        }
      }
    }
  }
} //fixSurfConnectivity


// get the sub-area of current structure
void AtomArray::getSubArea()
{

  if (surface==NULL || surfHeader.numSurfTris<=0 || InParm.subsurf==NULL )   // the structure should be read in already
  {
    totalSubArea=0.;
    numSubAreaAtoms=0;
    return;
  }

  int lastAtom=-1, atomcnt=0, i, j, atomNum;
  double areasum=0.0, ar=0, arsum=0;

  // loop init
  Coordinate pt1, pt2, pt3;
  SurfTriangle  *sptr;

  // calculate area of a triangle each time
  if (InParm.showSubSurf)
    printf("\nSurface area by atom (prob=%.2f) : \n----------\n", InParm.surfprobe);
  for (i=0, j=0, sptr=surface; i<surfHeader.numSurfTris; i+=3 )
  {
    // get atomNum and check if it's on the selected sub-area
    atomNum = sptr->atomNum; 
    if ( inPick(atomNum, InParm.subsurf) )       // is picked sub-surface?
    {
      pt1.set(sptr->vt[0]);	// 1st vertex vector
      pt2.set(sptr->vt[1]);	// 2nd vertex vector
      pt3.set(sptr->vt[2]);	// 3rd vertex vector

      // caculate area
      ar = 0.5*((pt1-pt3)*(pt2-pt3)).len();   // area  
      if ( lastAtom != atomNum  )
      {                                       // output previous atom 
         if ( InParm.showSubSurf && (lastAtom>=0) ) 
         {
           atomcnt++;
           printf("atom#%d (%4s, r=%#9.5lf):\t%#9.5lf\n", lastAtom+1, atomList[lastAtom].atomName, atomList[lastAtom].radius, arsum);
         }
         arsum=0.;
         lastAtom = atomNum;
      }
      arsum += ar;
      areasum += ar;

    }    
  }
  if ( InParm.showSubSurf && (lastAtom>=0)  )
  {
    atomcnt++;
    printf("atom#%d (%4s, r=%#9.5lf):\t%#9.5lf\n", lastAtom+1, atomList[lastAtom].atomName, atomList[lastAtom].radius, arsum);
  }
  numSubAreaAtoms=atomcnt;
  totalSubArea=areasum;
  if ( InParm.showSubSurf)
    printf("----------\nTotal area : %lg\n", areasum );

} //getSubArea()

float *pointsPlane(int nPoints, Coordinate *c, char cutaxe='Z');

// set AtomArray::Coeff to the cutoff plane
// if isAxe=1:  use InParm.cutoffaxe=X,Y,Z as cutoff axe
//    isAxe=0:  use last picked three atomsm, if no pick, use default axe
//
void AtomArray::setSurfCutoffPlane(int isAxe)
{
  Coordinate pickedCrds[3];
  if ( isAxe==0 )
    nSurfPlanePoints = State::PickQ.getPickedCrds(pickedCrds);
  if (nSurfPlanePoints <=0 || isAxe)
  {
    InParm.cutoffaxe= 'X' + (InParm.cutoffaxe-'X'+1)%3;   // move to next axe
    pickedCrds[0].set(0., 0., 0.);   // use CoMass for cutoff instead
    nSurfPlanePoints=1;
  }
  memcpy(Coeff, pointsPlane(nSurfPlanePoints, pickedCrds, InParm.cutoffaxe), sizeof(float)*4);   // set Coeff
  InParm.cutoffval=0.;          // reset cutoff
}

void AtomArray::flipSurfCutoffPlane()
{
  for (int i=0; i<4; i++)
    Coeff[i] = -Coeff[i];
  InParm.cutoffval = -InParm.cutoffval;
}

// display surface for PickSphere Cut
void AtomArray::SphereNet(int colorFilter) 
{
  GLUquadricObj *qobj=gluNewQuadric();
  gluQuadricDrawStyle(qobj, GLU_LINE);
  glEnable(GL_LIGHTING);
  glEnable(GL_LINE_STIPPLE);
  glLineWidth(0.1);
  glLineStipple(2, 0xAAAA);
  tRGBA c(75,75,0);
  glColor3ub(Color2R(colorFilter,&c), Color2G(colorFilter,&c), Color2B(colorFilter,&c)); 
  gluSphere(qobj, InParm.pickradius, 20, 20);
  glDisable(GL_LINE_STIPPLE);

}

// display surface image from memory    
#define GrayRGB  160
#define RGBup(charge)    (GrayRGB+(255-GrayRGB)*(charge>0?charge/maxAtomCharge:charge/minAtomCharge))
#define RGBdown(charge)  (GrayRGB*(charge>0?1-charge/maxAtomCharge:1-charge/minAtomCharge))
#define surfG(charge)    ((int)(charge!=0?RGBdown(charge):GrayRGB))
#define surfR(charge)    ((int)(charge>0?RGBup(charge):surfG(charge)))
#define surfB(charge)    ((int)(charge<0?RGBup(charge):surfG(charge)))

void AtomArray::surf(int colorFilter)
{
  if (surface == NULL)
    return;

  float   cutPlane[4]={Coeff[0],Coeff[1],Coeff[2],Coeff[3]+InParm.cutoffval};  // move along Normal

  int    *aacolor, gg;
  int     numTriangles=0, picked=0;
  GLubyte  transparent = InParm.transparent;
  Coordinate crd1, crd2, crd3;
  SurfTriangle *sptr;

  glEnable(GL_LIGHTING);
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  aacolor=InParm.AAcolor[AminoAcids::AA_UNI];
  tRGBA c(aacolor[0],aacolor[1],aacolor[2]);
  glColor4ub(Color2R(colorFilter,&c),Color2G(colorFilter,&c),Color2B(colorFilter,&c),transparent);  
  if ( InParm.showSurfInMesh )
  {
    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE );
    glLineWidth(0.1);
  }
  else
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL );

  glBegin(GL_TRIANGLES);

  for (int i=0; i<surfHeader.numSurfTris; i++ ) 
  {
    sptr  = surface + i;      // note: 1) sptr changed in loop
                              //       2) in 3rd atom# position to store flag for surface (=-1) 

    // check triagle for Cutoff
    int surfAtom = sptr->atomNum;  // point to atomNum
    if (InParm.showSubSurf && inPick(surfAtom, InParm.subsurf) )
    {
       numTriangles=1;        // display any way
       picked=1;
    } 
    else if (InParm.showCavity && ( sptr->out!=0 ) )
    {
      numTriangles=0;
      picked=0;
    }
    else
    {
       numTriangles=triangleToDraw(&sptr, cutPlane);  // sptr could be changed in the subroutine
       picked=0;
    }

    if ( numTriangles != 0 )                      	// draw one/two triangles
    {    
      GLubyte  mytrans=picked?255:transparent;  	// highlight selected area
      if ( picked ) 
	  c.set( 255,255,0);     	// yellow
      else if ( InParm.AAcolorGroups != SURF_COLOR_UNI ) 
      { 
        if ( InParm.AAcolorGroups == SURF_COLOR_RES ) {
          aacolor = InParm.AAcolor[AminoAcidsObj.GetType((char*)(resName[atomList[surfAtom].rn]))]; 
          c.set(aacolor[0],aacolor[1],aacolor[2]);
        } else {
          float chrg = atomList[surfAtom].charge;
          gg = surfG(chrg);
          if (gg <0)  
            gg=0;
          c.set(surfR(chrg),gg,surfB(chrg)); 
        }
      } 
      else if ( InParm.showSubSurf )  
        c.set(aacolor[0],aacolor[1],aacolor[2] ); 
      glColor4ub(Color2R(colorFilter,&c),Color2G(colorFilter,&c),Color2B(colorFilter,&c),transparent);  
     
      while ( numTriangles>0 )
      {
        for (int j=0; j<3; j++)            // draw three vertex 
        {
          glNormal3fv(sptr->nm[j]);        
          glVertex3fv(sptr->vt[j]);        
        }
        sptr++;
        numTriangles--;
      }
    } 
  }
  glEnd();

  if ( InParm.showSubSurf == true )
  {
    if (InParm.AAcolorGroups == SURF_COLOR_UNI  && InParm.subsurf != NULL )
    {
      char msg[200];
      sprintf(msg, "%d selected atoms on surface. area=%g", numSubAreaAtoms, totalSubArea );
      MsgBrd.set(msg);
    }
    else
      MsgBrd.erase();
  }
  glDisable(GL_BLEND);
}

// check the surface triangle to see if the cutoff plane goes through it.
// if the plane crosses the triangle, then split the triangle to draw a
// smoth cutoff edge.
// vertexcptr : points to 3 lines of surface vertexes ( a triangle ) 
//             with the size of TriBlockBytes
//
const int TriCrdBytes   = sizeof(float)*3;

int AtomArray::triangleToDraw(SurfTriangle **surfPtr, float *cutoffPlane)
{
  static SurfTriangle vertex[3];              	// 2 triangles (6 vertexes) and a highlight line(not implemented yet)
  Coordinate coeff(cutoffPlane);                // coeff, constD : the cutoff plane   [coeff]^X=constD
  Coordinate crd[3];
  bool  cut[3]={false, false, false};
  int   cnt=0, i, j1, j2;                  	// cnt: number of points to cut
  float crossPt1[3], crossPt2[3];         	// two cross points

  // apply cutoff
  if (InParm.cutoffval < 1.e+10 )		// cutoff not set if (cutoffval>=1.e+10)
  {
    for (i=0; i<3; i++)
    {
      crd[i].set( (*surfPtr)->vt[i]) ;   // triangle crds 
      cut[i]= ((coeff^crd[i])>cutoffPlane[3]);    // draw those only Ax+By+Cz<=constD 
      if ( cut[i] ) 
        cnt++;                   		//cuted triangle count
    }
  }

  switch(cnt)	  // map "cutted triangle count" to "output triangle count"
  {							
  case 0: 
    cnt=1;		  // keep the triangle;   switch cnt meaning: output triangle count
    break;
  case 1:           // cutoff one vertex,  result in 2 output triangles
    for ( i=0; (i<3)&& !cut[i]; i++);     // i is the vertex to drop
    j1 = (i+1)%3;                         // j1, j2, vertexes to keep
    j2 = (i+2)%3;
    if ( linePlaneCrossPoint(cutoffPlane, crd[i], crd[j1], crossPt1) &&
         linePlaneCrossPoint(cutoffPlane, crd[i], crd[j2], crossPt2) )
    {
       // split into two triangle
       // copy vertex block
       memcpy(vertex,   *surfPtr, SURF_TRI_LEN);
       memcpy(vertex+1, *surfPtr, SURF_TRI_LEN); 
   
       // first triangle : substitue i vertex crd to crossPt1
       memcpy(vertex[0].vt[i], crossPt1, TriCrdBytes); 

       // second triangle: add [crossPt2,corsPt1,j2] triangle (keep the j2)
       memcpy(vertex[1].vt[i], crossPt2, TriCrdBytes); 
       memcpy(vertex[2].vt[j1], crossPt1, TriCrdBytes);
       *surfPtr=vertex;
       cnt=2;                          // output trangle count
    } 
    else
      cnt=1;          			// should never be here
    break;
  case 2:  // cutoff two vertex, result in 1 output triangle
    cnt = 1;
    for ( i=0; (i<3) && cut[i]; i++);     // i is the vertex to keep
    j1=(i+1)%3;                           // vertexes to drop
    j2=(i+2)%3;
    if ( linePlaneCrossPoint(cutoffPlane, crd[i], crd[j1], crossPt1) &&
         linePlaneCrossPoint(cutoffPlane, crd[i], crd[j2], crossPt2) )
    {
      memcpy(vertex,   *surfPtr, SURF_TRI_LEN);
      // substitue the i+1 and i+2 vertexes in the first triangle
      memcpy( vertex->vt[j1], crossPt1, TriCrdBytes); 
      memcpy( vertex->vt[j2], crossPt2, TriCrdBytes); 
      *surfPtr=vertex;
    }
    break;
  case 3 : // cut off all vertexes, 0 output triangle
    cnt=0;
    break;
  default: // keep the triangle
    cnt=1;
    break;
  }
  return cnt;
}


// Ramachandran diagram of Dihedral angles to determine protein secondary structure.
// Do not consider H-bond.
//
// .secStruct is for Dihedral angles between Residue(i) & Residure(i+1), donot count the residue
// only first atom was assigned now
//                     res1           res2           listi
//                     |              |              |
// .----.--.--.--------.--.--.--------.--.--.--------.
// C    N  CA C        N  CA C        N  CA C        N
// 
// secStruct = 0x[phi<0][psi<0]
//
const double MaxCADist=5.0;        // max distance between CAs in single chain
const double MaxNOdist=3.5;        // max H-bond distance
const int    MaxICount=5;
void AtomArray::setSecStructStep1()  // first atom of residue has secStructType
{
  Coordinate N1, CA1, C1, N2, CA2, C2, O2;   // Atom Crds
  int        foundCA1(0), foundN1(0), foundC1(0);
  int        foundCA2(0),  foundN2(0), foundC2(0), foundO2(0);  // save previous two residue info
  int        numResRead=0;      
  Atom       *listi;
  Atom       *res2=NULL, *res1=NULL; 
  double     phi, psi;     // torsional angles
  int        i;

  foundO2=0;

  for (i=0,listi=atomList; i<numAtom; i++, listi++)
    listi->secStruct=STRUCT_UNKNOWN;              // reset to unknown

  for (i=0,listi=atomList,res1=NULL,res2=atomList; i<=numAtom; i++, listi++)  // i==numAtom for last residue
  {
     // check for a new residue, no read in atom info yet
    if (i==numAtom || res2->rn != listi->rn )     // get a new residue, check previous residue info 
    {
      if (foundCA2==0 )
      {
        if ( res1 != NULL)           
          res1->secStruct |=  STRUCT_CHAINEND; 
        res2->secStruct |=  STRUCT_CHAINEND; 
      }
      else if (res1!=NULL)
      { 
        if ( foundCA1==0)
          res1->secStruct |= STRUCT_CHAINEND;       // read in first CA of the chain
        else if  ( ((CA1-CA2).len())>=MaxCADist || res1->chainIndex!=res2->chainIndex )    // must (foundCA1&foundCA2)  
          res1->secStruct |= STRUCT_CHAINEND; 
      }
      else
       foundCA1=0;

      if (foundC1==1 && foundN2 && foundCA2 && foundC2 )
      { 
         phi=angleT(C1, N2, CA2, C2);
         if (phi <= 0.0 )    
         {
           res2->secStruct |= STRUCT_PHI;
           if (res1==atomList)
             res1->secStruct |= STRUCT_PHI;
         }

        //if ( res2->rn> 90 && res2->rn< 100 ) 
        // printf("atom=%d res=%d phi=%g\n", i, res2->rn, phi);

      }
      if ( foundN1 && foundC1==1 && foundN2 )
      {
         psi=angleT( N1, CA1, C1,  N2);
         if (psi >= -170. && psi <= 0.0 )     // [-170, 0]
           res1->secStruct |= STRUCT_PSI;
         if (i==numAtom)
           res2->secStruct |= STRUCT_PSI;

         //if ( res1->rn> 90 && res1->rn< 100 ) 
         //printf("atom=%d res=%d psi=%g\n", i, res1->rn, psi);
      }

      if ( foundO2)
         res2->secStruct |= STRUCT_OFOUND;

      if ( i==numAtom)
      {
        res2->secStruct |= STRUCT_CHAINEND; 
        return;
      }

      if ( sIsChainEnd(res2->secStruct) )
      {
        res1=NULL;
        foundCA1=0;
        foundC1=0;
        foundN1=0;
      } else {
        res1=res2;
        foundCA1=foundCA2;
        foundC1=foundC2;
        foundN1=foundN2;
      }
      foundCA2=0;
      foundC2=0;
      foundN2=0;
      res2=listi;
      foundO2=0;
    } //     if (res1->rn != listi->rn )     // get a new residue

    if( !strncmp(listi->atomName, "CA", 2))           // consider CAE, CAx
    {
      CA1 = CA2;
      CA2 = *listi;
      foundCA2=1;  
      //printf("foundCA2 %d\n", listi->rn);
    } 
    else if ( !strcmp(listi->atomName, "N" ))  
    {
      N1=N2;
      N2= *listi;
      foundN2=1;
    }   
    else if ( !strcmp(listi->atomName, "C" ))  
    {
      C1=C2;
      C2 = *listi;
      foundC2=1;
    }   
    else if ( !strcmp(listi->atomName, "O" ))  
    {
      O2 = *listi;
      foundO2=1;
    }   
  } // for(i)  

}  // setSecStructStep1()

// check input pick list to see if the atom is in the pickList
// return 1 if is otherwise return 0
// if getColor=true,  enter path to get color attribute, set picked color to inPickColor
// note:  input pick list always use order(start from 1) instead of index(start from 0)
//
bool AtomArray::inPick(int atomindex, char *pickList, bool getColor)
{
  int currInPick=0;
  int lastInPick=0;
  int inPickSum=0;
  char op='|';
  char *cptr=pickList;
  int  state=1000;
  int  nlen=0;           /* total name length */
  int  wpos=0;           /* wild card * position */
  int  nBegin, nEnd;
   
  inPickColor = InParm.ribbonColorBase;
  if (cptr!=NULL && atomindex>=0 && atomindex<numAtom)
  {
    int rn= atomList[atomindex].rn;    // res from index to order

    while (1)
    {
      switch (state)
      {
      case 1000:    // set op
        while (op==' ')
        {
          
          if (*cptr==',' || *cptr=='|' || *cptr=='&' || *cptr=='\0')
          {
            op = *cptr;
            break;
          }
          cptr++;
        }
        state=1001;
        break;
      case 1001:
        if (lastInPick==1 && op=='|' )
        { 
          while (*cptr!='\0' && *cptr!='&' && *cptr!=',')
             cptr++;
          op=*cptr;
          cptr++;
          state=2000;
        }
        else
          state=1002;
        break;
      case 1002:
        if (getColor==false && lastInPick==1 && op==',')
        {
          state=1000;
          op='\0';       //ending
        }
        else
          state=2000;
        break;
      case 2000:
        if (strncmp(cptr, "mono ", 5)==0)
          state=2001;
        else if (strncmp(cptr, "prtc ", 5)==0)
          state=2002;
        else if (strncmp(cptr, "#mon ", 4)==0)
          state=2003;
        else if (strncmp(cptr, "#prt ", 4)==0)
          state=2004;
        else
          state=2000;
        if ( *cptr == '\0')
           state = 1000;  // ending
        else if (state!=2000)
          cptr += 4;
        else
          cptr++;
        break;
      case 2001:     // chem mono XXXX
        nlen=0;
        wpos=-1;
        while (!isspace(*cptr) && *cptr!='\0' && *cptr!=',' && *cptr!='|' && *cptr!=',' && nlen<MAXLEN_NAMES)
        { 
          if ( *cptr=='*' && wpos==-1 )          
            wpos=nlen;            // consider widecard *
          cptr++;
          nlen++;
        }
        if ( wpos==-1 ) 
           wpos=nlen;
        //strncpy(name, cptr-nlen, wpos);
        if (wpos==0 || strncmp(resName[rn], cptr-nlen, wpos)==0 )
          currInPick=1;
        else 
          currInPick=0;
        state=3000;
        break;
      case 2002:    // chem prtc XXXX
        nlen=0;
        wpos=-1;
        while (!isspace(*cptr) && *cptr!=',' && *cptr!='|' && *cptr!=',' && nlen<MAXLEN_NAMES)
        { 
          if ( *cptr=='*' && wpos==-1 )
            wpos=nlen;
          cptr++;
          nlen++;
        }
        if (wpos==-1 ) 
           wpos=nlen;
        //strncpy(name, cptr-nlen, wpos);
        if (wpos==0 || strncmp(atomList[atomindex].atomName, cptr-nlen, wpos)==0 )
          currInPick=1;
        else
          currInPick=0;
        state=3000;
        break;
      case 2003:    // #mono ddd ddd
        if ( sscanf(cptr, "%d%d", &nBegin, &nEnd) == 2 )
        {
          int rnPlus = rn+1;                      // res order to res index
          if ( nBegin<=rnPlus && rnPlus<=nEnd )       
            currInPick=1;
          else 
            currInPick=0;
          while (isdigit(*cptr) || isspace(*cptr)) 
            cptr++;
          state=3000;
        }
        else
          state=2003;
        break;
      case 2004:   // #prt ddd ddd
        if ( sscanf(cptr, "%d%d", &nBegin, &nEnd) == 2 )
        {
          int indexPlus = atomindex+1;
          if ( nBegin<=indexPlus && indexPlus<=nEnd)
            currInPick=1;
          else 
            currInPick=0;
          while (isdigit(*cptr) || isspace(*cptr)) 
            cptr++;
          state=3000;
        }
        else
          state=2003;
        break;
      case 3000:
        if ( getColor==true && *cptr == '<' )
          state=4000;
        else
          state=5000;
        break;
      case 4000:
        if ( currInPick == 1 )
        { 
          if (sscanf(cptr, "<%d-%d-%d>", &inPickColor.r, &inPickColor.g, &inPickColor.b) == 3 )
          {
            state = 5000;
            inPickColor.setGrey();
          }
          else 
            state=4000;
        }
        else
          state=5000;
        cptr=strchr(cptr, '>');
        cptr++;
        break;
      case 5000:
        if ( op == '&' )
           lastInPick &= currInPick;
        else
           lastInPick |= currInPick;
        if ( (op==',' || op=='\0' || isspace(op) ) && lastInPick==1 && inPickSum==0 )
           inPickSum=1;
        state=1000;
        op=' ';
        break;
      default:
        op='\0';
      } //switch()
      if (state==1000 && (cptr==NULL || op == '\0' || *cptr=='\0') )
        break;  
      while ( isspace(*cptr) ) 
        cptr++;        // skip blanks
    } // while()
  } // if() 

  if ( (inPickSum|lastInPick) == 0  ||  *(atomList[atomindex].atomName) == '\0' )
    return false;
  else 
    return true; 

  //return inPickSum|lastInPick;

} // inPick

int AtomArray::GetInputParameters(char *inputfile, int crdIndex, char *cmoil) 
{
  return InParm.ReadInput(inputfile, crdIndex, cmoil);
}

void AtomArray::SubtractCoMass(double cX, double cY, double cZ)
{
  CoMX=cX;
  CoMY=cY;
  CoMZ=cZ; 
 
  // keep coordinates for XYZ format 
  if ( CoMX != 0.0 || CoMY !=0.0 || CoMZ !=0.0)
  {
    Atom *aptr = atomList;
    for (int i=0; i<numAtom; i++, aptr++)
    {
      aptr->x -= CoMX;
      aptr->y -= CoMY;
      aptr->z -= CoMZ;
    }
  }
}

void AtomArray::SubtractCoMass()
{
  double sumX(0), sumY(0), sumZ(0);
  Atom *aptr;
  int counter;

  aptr = atomList;
  for (counter=0; counter<numAtom; counter++, aptr++)
  {
    sumX=sumX+aptr->x;
    sumY=sumY+aptr->y;
    sumZ=sumZ+aptr->z;
  }
  CoMX=sumX/double(numAtom);
  CoMY=sumY/double(numAtom);
  CoMZ=sumZ/double(numAtom); 
 
  // keep coordinates for XYZ format 
  if ( CoMX != 0.0 || CoMY !=0.0 || CoMZ !=0.0)
  {
    aptr = atomList;
    for (counter=0; counter<numAtom; counter++, aptr++)
    {
       aptr->x-=CoMX;
       aptr->y-=CoMY;
       aptr->z-=CoMZ;
    }
  }
}

void AtomArray::ResetCrd()
{
  Atom *aptr = atomList;
  for (int i=0; i<numAtom; i++, aptr++)
  {
    aptr->x = aptr->orgX;
    aptr->y = aptr->orgY;
    aptr->z = aptr->orgZ;
  }
}

void AtomArray::FindmaxAngle()
{
  angle=(40.0/180)*PI;
}
void AtomArray::FindmaxZ()
{
  Atom *aptr=atomList;
  int counter;
  
  if ( aptr == NULL)
    return;

  maxZ=0;
  double tempz(0.0);
  for (counter=1; counter<numAtom; counter++, aptr++)
  {
    if (aptr->atomName[0]!='H' && aptr->atomName[1]!='H' )   // no lowercase
    {
      tempz=sqrt(aptr->x*aptr->x+aptr->y*aptr->y+aptr->z*aptr->z);
      //tempz=(sqrt((aptr->x)*(aptr->x)+(aptr->y)*(aptr->y))/(tan(angle)))+aptr->z;
      if(tempz>maxZ)
	  maxZ=tempz;
    }
  }
  minZ = -maxZ ;
}

bool AtomArray::checkNeigh(int no, char *a,int *sn)
{
  int i, nn;
  for (i=0; i<atomList[no].numNeighbors; i++)
    {
      nn=atomList[no].neighbors[i];
      if(strcmp(atomList[nn].atomName, a)==0)
	{
	  *sn = nn;
	  return 1;
	}
    }
  return 0;
}

bool AtomArray::sameatom(char *a, char *b) {
  size_t i,j;
  for(i=0;i<strlen(b);i++) {
    if (a[i] != b[i]) {
      return false;
    }
  }	
  for(j=i;j<strlen(a);j++) {
    if (!isdigit(a[j])) {
      return false;
    }
  }
  return true;
}

void AtomArray::setColor()
{
  int i;
  char tn[MAXLEN_NAMES];
  char *settings=NULL;
  tRGBA *c;
  int   counter(0);
  GLubyte whiteBG=0 ;
  GLubyte colorC;   // =whiteBG*0.3 + 0.7*(1-whiteBG); 
  GLubyte colorH;
  GLubyte colorUnknown;

  if ( InParm.bg.r>165 && InParm.bg.g > 165 && InParm.bg.b > 165 )
     whiteBG=255;
  colorC=(GLubyte)(whiteBG*0.45+0.70*(255-whiteBG));    // adjust upon bg
  colorH=(GLubyte)(whiteBG+(255-whiteBG)*0.60); 
  colorUnknown=(GLubyte)(whiteBG*0.55+0.65*(255-whiteBG)); 

  settings=InParm.colorsettings;
  for(i=0; i<numAtom; i++)
  {
      strcpy(tn,atomList[i].atomName);
      c = &atomList[i].color;
      if (settings!=NULL && inPick(i, settings, true) )
      {
	  c->r=inPickColor.r;  
	  c->g=inPickColor.g;   
	  c->b=inPickColor.b;   
      }
      else if ( *tn == 'C' || 
	  sameatom(tn, "C") ||
	  sameatom(tn, "CA")||
	  sameatom(tn, "CB")||
	  sameatom(tn, "ME")||
	  sameatom(tn, "CG")||
	  sameatom(tn, "CD")||
	  sameatom(tn, "CH")||
	  sameatom(tn, "CE")||
	  sameatom(tn, "CZ")||
	  sameatom(tn, "CAE")||
	  sameatom(tn, "CBE")||
	  sameatom(tn, "CN") )
	{
	  c->r=colorC;   // 0.83;    carbon grey
	  c->g=colorC;   // 0.83;
	  c->b=colorC;   // 0.83;
	} 
      else if ( *tn == 'N' || sameatom(tn, "N")||
	       sameatom(tn, "NE") ) 
	{
	  c->r=170; // 0.68;
	  c->g=210; // 0.85;
	  c->b=230; // 0.90;
	} 
      else if (sameatom(tn, "P")) 
	{
	  c->r=255; // 1.0;
	  c->g=165; // .65;
	  c->b=0;   // 0.0;
	}
      else if ( *tn == 'O' || sameatom(tn, "O")||
	       sameatom(tn, "OE")||
	       sameatom(tn, "OA")||
	       sameatom(tn, "OS")||
	       sameatom(tn, "ON")||
	       sameatom(tn, "OH") )
	{
	  c->r=255; // 1.0;
	  c->g=0;   // 0.0;
	  c->b=0;   // 0.0;
	}
      else if ( *tn == 'H' || ( isdigit(*tn) && *(tn+1)=='H') || sameatom(tn, "H")||
	       sameatom(tn,"HOE")||
	       sameatom(tn, "HE") )
	{
	  c->r=colorH;
	  c->g=colorH;
	  c->b=colorH;
	}
      else if (sameatom(tn, "NA"))
	{
	  c->r=100;  // .39;
	  c->g=150;  // .58;
	  c->b=240;  // .93;
	}
      else if (sameatom(tn, "CL"))
	{
	  c->r=0;   // 0.0;
	  c->g=255; // 1.0;			
	  c->b=0;   //0.0;			
	}
      else 
	{
	  //printf("I do not know what type of atom %s is\n",tn);
          c->r=colorUnknown;
	  c->g=colorUnknown;
	  c->b=colorUnknown;
	}
      c->setGrey();
  }  // for()

  if ( ribbonColor != NULL )
  {
    int rn, currRn=-1;
    settings=InParm.ribboncolor;
    c = ribbonColor;
    for (i=0; i<numAtom; i++)
    {
      rn = atomList[i].rn;
      if ( currRn != rn )                // assign color to a new residue
      {
        if ( inPick(i, settings, true) )
        {
          c->r = inPickColor.r;
          c->g = inPickColor.g;
          c->b = inPickColor.b;
        }
        else 
          *c=InParm.ribbonColorBase;
        c->setGrey();
        currRn = rn;
        c++;
      } 
    }
  } 
}

int AtomArray::DispMode()
{
  return InParm.dispmode;
}

int AtomArray::BondMode()
{
  return InParm.bondmode;
}

// get operator from input 
char* AtomArray::NextOp(char *astring)
{
  if (astring == NULL)
    return NULL;

  while ( *astring != '\0' && *astring != '|' && *astring != '&' && *astring != '<' )
     astring++;

  if ( *astring == '\0' )
    astring = NULL;

  return astring;
}

// set align field upon input overlap settings
//
void AtomArray::setAlign()
{
  if (InParm.alignsettings != NULL )
  {
    for (int m=0; m<numAtom; m++)
    {
      if ( atomList[m].skip==0 && inPick(m, InParm.alignsettings) )
      {
        atomList[m].align=true; 
      }
    }
  }
  else			// no alignment input, default to backbone alignment
  {
    for (int i=0; i<numAtom; i++)
      if ( !strcmp(atomList[i].atomName, "CA" ) )
        atomList[i].align=true;		// default is false
  }
}

// if chainID=='-' output all chains otherwise output selected chain
void AtomArray::OutputCrdInPDB(FILE *crdf, int chainIndexID, char *extraTxt)
{
  if (crdf == NULL)
    return;

  fprintf(crdf, "HEADER    %-40.40s00-XXX-00   %-4.4s\n", "CMOIL PDB", moleculeName);
  fprintf(crdf, "CRYST1%9.3f%9.3f%9.3f%7.2f%7.2f%7.2f %-11.11s%4d\n",1.,1.,1.,90.,90.,90.,"P 1",1);
  fprintf(crdf, "REMARK   1 This file is automatically generated by CMOIL\n");

  Atom *atom=atomList, *catom=NULL;
  char chain=' ';
  int  lastChain=-2;   // for the chainid, default is -1
  int  k, lastResNum=-1, serNum=0, numPrintRes=0, lastRn=-1, numChainRes=0, lastChainRes=-1, chainindex;

  // print SEQRES section
  for(k=0; k<numAtom; k++, atom++)
  {
    if ( *(atom->atomName) == '\0' )
      continue;

    chainindex = atom->chainIndex; 
    if ( chains != NULL) 
       chain=(chains+chainindex)->id;

    if ( chainIndexID<0 || chainIndexID==chainindex)
    {
       if ( atom->rn != lastRn )          // output a new res
       {
         if ( chainindex != lastChain )        // start a new chain
         {
            lastChain = chainindex;
            numPrintRes=0;
            serNum=0;

            // count total number of residues in a chain
            catom=atom;
            numChainRes=0;
            lastChainRes=-1;
            for (int i=k; i<numAtom && (catom->chainIndex==chainindex); i++, catom++)
            {
               if ( catom->rn != lastChainRes )
               {
                  numChainRes++;
                  lastChainRes=catom->rn;
               }
            } 
         }
         if ( strncmp(resName[atom->rn], "HOH", 4) && strncmp(resName[atom->rn], "TIP3", 4) )
         {     
           if ( (numPrintRes%13) == 0  )  {      // 13 residues per line
             if (atom != atomList )
               fprintf(crdf, "\n");
             fprintf(crdf, "SEQRES  %2d %c %4d ", ++serNum, chain, numChainRes);
           }
           fprintf(crdf, "%4.4s", resName[atom->rn]);
           lastRn=atom->rn;
           numPrintRes++;
         }
       }
    }
  }
  fprintf(crdf, "\n");    // end of SEQRES

  chain=' ';
  atom=atomList;
  for(k=0; k<numAtom; k++, atom++)
  {
    if ( *(atom->atomName) == '\0' )
      continue;

    chainindex = atom->chainIndex; 
    if ( chains != NULL) 
       chain=(chains+chainindex)->id;

    if ( chainIndexID<0 || chainIndexID==chainindex)
    {
      if ( strlen(atom->atomName)<4)
      {
        fprintf(crdf,"ATOM  %5d  %-3.3s %-4.4s%c%4d    %8.3f%8.3f%8.3f%26.26s\n",
                      k+1, atom->atomName, resName[atom->rn], chain, atom->rn+1, atom->orgX, atom->orgY, atom->orgZ, " "); 
      } else {
        fprintf(crdf,"ATOM  %5d %-4.4s %-4.4s%c%4d    %8.3f%8.3f%8.3f%26.26s\n", 
                      k+1, atom->atomName, resName[atom->rn], chain, atom->rn+1, atom->orgX, atom->orgY, atom->orgZ, " "); 
      }
    }
  }
}

// CRD format from fortran:
//	  write(ucrd,103,err=7)i,poimon(i),char1,char2,xtmp
//     1		,ytmp,ztmp,char3,char4,moretmp
//103	  format(i5,i5,1x,a4,1x,a4,3(f10.5),1x,a4,1x,a4,f10.5)	 
//
void AtomArray::WriteOrgCrd(char *extraTxt)
{
  // take the current crd name as base name
  char  tmpname[FILENAME_MAXLEN], *cptr;
  strcpy(tmpname, InParm.crdName);
  cptr=strrchr(tmpname, '.');
  if (cptr != NULL)
  {
    if (strchr(cptr, '\\') == NULL && strchr(cptr, '/')==NULL )
      *cptr = '\0';            // delete the suffix
  }
  cptr = tmpname + strlen(tmpname);
  sprintf(cptr, "_ml_%d", currentmol);

  if (InParm.filetag == fPDB )      // compose the filename to save
    strcat(tmpname, ".pdb");
  else 
    strcat(tmpname, ".crd");

  FILE *crdf = fopen(tmpname, "w");
  if (crdf == NULL)
  {
    fprintf(stderr, "ERROR: Unable to open file for writing: %s\n", tmpname);
    return;
  }

  if (InParm.filetag == fPDB) 
  {     // write as PDB if original file is PDB
    OutputCrdInPDB(crdf, -1, extraTxt);
  } 
  else 
  {    // otherwise write current structure to a CRD file
    fprintf(crdf, "* This file is automatically generated by CMOIL\n");
    fprintf(crdf, "* %.20s\n", moleculeName);
    fprintf(crdf, "* %s\n", extraTxt);
    fprintf(crdf, "%5d\n", numAtom);

    Atom *atom=atomList;
    for(int k=0; k<numAtom; k++, atom++)
    {
      fprintf(crdf,"%5d%5d %-4.4s %-4.4s%10.5f%10.5f%10.5f %4.4s FREE%10.5f\n", 
         k+1, atom->rn+1, resName[atom->rn], atom->atomName, atom->orgX, atom->orgY, atom->orgZ, moleculeName, 0.);
    }
  }
  fclose(crdf);

  DisplayMsgWrtCrd(tmpname);
}

void AtomArray::DisplayMsgWrtCrd(char *filename)
{
  MsgBrd.set(Ov.structs[1]);
  MsgBrd.append(" coordinates saved to ");
  MsgBrd.append(filename);
  glutPostRedisplay();
}

void AtomArray::DisplayStructure(int &imgIndex, int colorFilter, int lineWidth) 
{
  if ( (dispmodeComb&(STICK|SPACEBALL|STICKBALL)) == 0 ) 
    return;  // do not display basic structure

  glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

  int i, typeno=-1;
  Atom *listi; 
 
  for (i=imgIndex,listi=atomList; i<numAtom+imgIndex; i++, listi++)
  {
    if (listi->skip!=0)
      continue;
    if (InParm.filetag== fXYZ )
      typeno=1;
    else if ( PickMode==1 && (listi->dispmode&UNITBALL)!=0 )
    {
      typeno=-1;      // picked, don't show real atom sphere
    }
    else 
    {
      if  ( (listi->dispmode&SPACEBALL)!=0 ) 
        typeno =int( listi->radius*100+.5); // Get the atom type that we want to draw
      else if ( (listi->dispmode&STICKBALL)!=0 )
        typeno=1;
      else
        typeno=-1;
    }
       
    if ( typeno >= 0 )           // draw a ball
    {
      glEnable(GL_LIGHTING);
      glColor3ub( Color2R(colorFilter, &(listi->color)),
                  Color2G(colorFilter, &(listi->color)),
                  Color2B(colorFilter, &(listi->color)) ); 
	glTranslatef(listi->x, listi->y, listi->z);	// Move to the proper location
      glCallList(typeno); 					// Render sphere display list
	glTranslatef(-listi->x, -listi->y, -listi->z);// Move back
    } //if
       
    if ( ( (listi->dispmode&(STICK|STICKBALL))!=0 && (listi->dispmode&SPACEBALL)==0 )  )   // Stick mode
    {
      Atom *listnn;

      if (DrawBondAsLine)
        glDisable(GL_LIGHTING);
      else
        glEnable(GL_LIGHTING);

      //glLineWidth((float)lineWidth);
      glLineWidth(InParm.stickWidthFactor);

	for (int j=0; j<listi->numNeighbors; j++)
      {
        listnn = &(atomList[listi->neighbors[j]]);
        glColor3ub( Color2R(colorFilter, &(listi->color)),
                    Color2G(colorFilter, &(listi->color)),
                    Color2B(colorFilter, &(listi->color)) ); 

        if (DrawBondAsLine)
        {
          glBegin(GL_LINES);
          glVertex3d(listi->x, listi->y, listi->z);
          glVertex3d((listi->x+listnn->x)*0.5, (listi->y+listnn->y)*0.5, (listi->z+listnn->z)*0.5);
          glEnd();
        }
        else
          Stick::Draw(listi, listnn);
	}
    }
  } // for(i)
  imgIndex += i;    // for overlap 
}

// only for dispmode
// if useOverlapMatrix==true, then apply OverlapMatrix to org data
// This method is call for 
// 1.  rebuild structural info
// 2.  reread crd to orgXYZ
//
void AtomArray::ReProcessData()
{
  if ( (InParm.filetag != fXYZ ) && (InParm.dispmode & SURF) != 0 )
  {
    getSurfaces();         // for SURF command, in ProcessData, it's called in ReadCoordinate();
  }
  if (orgCrdModified)
  {
    if ( currentmol >=1 )
      currentmol--;
    ReadNextStruct();
  }
 
  setDefaultDispMode();    // set all atom dispmodes to the overall displayMode

  setRibbon();
  setSpace();
  setStick();
}

// read coordinates for both atom position and surface vertex
void AtomArray::ReadNextStruct()
{
  if ( InParm.structend <= 0 )
  {
    ReadCoordinate();       // read current structure
    if (  (InParm.dispmode&SURF) != 0 )
      readSurface();        // read surface for currentmol
  }
  else
  {
    for (int i=0; i<InParm.structend; i++)
    {
       ReadCoordinate();       // read current structure
       if (  (InParm.dispmode&SURF) != 0 )
         readSurface();        // read surface for currentmol
    }
    InParm.structend=-1;
  }
}

void AtomArray::ProcessData(void)
{
  if ( InParm.filetag != fXYZ && InParm.filetag != fPDB )
    ReadConn(); 
 
  if ( (InParm.filetag != fXYZ ) && (InParm.dispmode & SURF) != 0 )
  {
    getSurfaces();        // get surface vertexes for all structures first
    currentmol=0;
  }

  ReadNextStruct();       // read in coordinates

  if (BondMode()==1)
    CalcNeighbors();

  FindmaxAngle();
  FindmaxZ();
  setColor();
  setSphere();
  setSkip();
  setAlign();
  setDefaultDispMode();       // set all atom dispmodes to the overall displayMode

  setRibbon();
  setSpace();
  setStick();
}
}

