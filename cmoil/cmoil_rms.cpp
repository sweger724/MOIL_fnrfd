//*************************************************************************************************
//*  Filename:   cmoil_rms.cpp
//*
//*  Description: 
//*    Compute new overlap coordinates for two structures.
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Jul 25, 2001	Baohua Wang	Extracted from nloopp code and add wrap function and fix plan problem
//*
//*************************************************************************************************
#include "cmoil.h"
#include "cmoil_rms.h"

namespace CMOIL {
OverlapAlign::OverlapAlign() 
{
  ovtype=' ';
}

//========================================================================
//========================================================================
// perform trans1->rotate->trans2 on aa's original coordinates orgX,Y,Z
// rot[3][3] is row major , trans1[3], trans2[3]
void OverlapAlign::OverlapMoveStruct(AtomArray &aa)
{
  Atom *myatom;
  double f1, f2, f3;
  double *trans1=&(OverlapMatrix[0][0]), *trans2=&(OverlapMatrix[1][0]), *rot=&(OverlapMatrix[2][0]);
  int   k;

  // translate1 for all atoms in aa
  for (k=0, myatom=aa.atomList; k<aa.numAtom; k++, myatom++)		
  {
    myatom->orgX += trans1[0];
    myatom->orgY += trans1[1];
    myatom->orgZ += trans1[2];
  }

  // rotate for all atoms in aa 
  for (k=0, myatom=aa.atomList; k<aa.numAtom; k++, myatom++) 
  {
    f1 = myatom->orgX*rot[0] + myatom->orgY*rot[1] + myatom->orgZ*rot[2];
    f2 = myatom->orgX*rot[3] + myatom->orgY*rot[4] + myatom->orgZ*rot[5];
    f3 = myatom->orgX*rot[6] + myatom->orgY*rot[7] + myatom->orgZ*rot[8];  
    myatom->orgX=f1;
    myatom->orgY=f2;
    myatom->orgZ=f3;    
  }

  // traslate2 for all atoms in aa
  for ( k=0, myatom=aa.atomList; k<aa.numAtom; k++, myatom++)
  {		                                              
    myatom->orgX += trans2[0];
    myatom->orgY += trans2[1];
    myatom->orgZ += trans2[2];
  }
  aa.orgCrdModified=true;
}

// reread the coordinates and perform overlap on the orgX,Y,Z directly
//
// translates[2] for aa1 and aa2 translates and rotate2 is aa2 ratation after translate
// returen number of overlapped atoms and rmsd
// 

//copy aa crd to match structure
int OverlapAlign::cpyAA2Match(AtomArray &aa, COORDINATES &match, bool cpyOrg)
{
  Atom *aptr;
  int n=0, i;

  aptr=aa.atomList; 
  for (i=0, n=0; (i<aa.numAtom)&&(n<=MAX_SEQ); i++, aptr++)
  {
    if ( aptr->align && (  ovtype=='D' || 
                           !strncmp(aptr->atomName, "CA",4) ||
                           ( ovtype=='E' && ( !strncmp(aptr->atomName, "C",4)  || 
                                              !strncmp(aptr->atomName, "N",4)  || 
                                              !strncmp(aptr->atomName, "O",4)   )  ) ) ) 
    {
      if ( cpyOrg)
      {
        match.x[n]=aptr->orgX;            // copy original crds
        match.y[n]=aptr->orgY;
        match.z[n]=aptr->orgZ;
      } else {
        match.x[n]=aptr->x;               // copy aliged crds
        match.y[n]=aptr->y;
        match.z[n]=aptr->z;
      }
      n++;
    }
  }
  return n;
}

//decode the overlap error code into text
void OverlapAlign::setErrMsg()
{
   switch (numAlignedAtoms)
   {
     case 0  : sprintf(ErrMsg, "%s,%s: -- unable to overlap", structs[0], structs[1]) ;
               break;
     case -1 : case -2 : 
               sprintf(ErrMsg, "%s,%s: -- unable to overlap: struct#%d '%s' has less atoms", 
                           structs[0], structs[1], -numAlignedAtoms,  structs[-numAlignedAtoms-1]);
               break;
     case -3 :
               sprintf(ErrMsg, "%s,%s: -- unable to overlap: None alignment input should be empty", 
                           structs[0], structs[1], -numAlignedAtoms,  structs[-numAlignedAtoms-1]);

   }
}


int OverlapAlign::overlap(AtomArray &aa1, AtomArray &aa2)
{
  COORDINATES match1, match2;
      
  printf("CMOIL: overlapping....\n");
  memset(OverlapMatrix, 0, sizeof(OverlapMatrix[5][3]));
  strncpy(structs[0],aa1.moleculeName,4);
  strncpy(structs[1],aa2.moleculeName,4);
  numAlignedAtoms = -1; 
  switch (ovtype)
  {
    case 'A':                             // auto alignment by CE
      autoAlign.overlap(aa1,aa2);
      rmsd = autoAlign.rmsd;
      return 1;
    case 'B': case 'C':                   // auto alignment by CE, then, align Ca or backbone
      autoAlign.matches(aa1,aa2,match1,match2, ovtype);
      numAlignedAtoms=autoAlign.numAlignedAtoms;
      sprintf(structs[0]+4, ":%c", aa1.displayChainid==-1?' ':(aa1.chains+aa1.displayChainid)->id);
      sprintf(structs[1]+4, ":%c", aa2.displayChainid==-1?' ':(aa2.chains+aa2.displayChainid)->id);
      break;
    case 'D': case 'E': case 'F':                 // align by selected atoms
      if ( aa1.InParm.alignsettingsLen>0 && aa2.InParm.alignsettingsLen>0 )
      {                                           // has alignment string from both structure
        int n1 = cpyAA2Match(aa1, match1); 
        int n2 = cpyAA2Match(aa2, match2); 
        if ( n1<n2 )
          numAlignedAtoms=-1;
        else if ( n1>n2)
          numAlignedAtoms=-2;
        else 
          numAlignedAtoms=n1; 
      } 
      else 
        numAlignedAtoms=-3;
      autoAlign.numAlignedMonos=-1;            // set flag of not using autoalign
      break;
  }
  if ( numAlignedAtoms <0 )
  {
    setErrMsg();
    return numAlignedAtoms;
  }

  int rtn;
  if ( (rtn=rms_dist(&match1, &match2, numAlignedAtoms)) >= 0 )
  {
    OverlapMoveStruct(aa2);

    // copy orgXs to x,y,z for display
    aa1.ResetCrd();   /* should be reread from orgX,Y,Z to x,y,z*/
    aa2.ResetCrd();

    // modified for display, aa2 takes aa1's mass center
    aa1.SubtractCoMass();
    aa2.SubtractCoMass(aa1.CoMX, aa1.CoMY, aa1.CoMZ);
    // no need to set structural info again

    // compute rmsd (Root Mean Square Deviation)
    if ( numAlignedAtoms > 0 )
    {
      double ff,rms=0.0;
      if ( ovtype <= 'C' )
      {
        autoAlign.alignAA2Match(aa1, match1, 1);
        autoAlign.alignAA2Match(aa2, match2, 2);
      } else {
        cpyAA2Match(aa1, match1); 
        cpyAA2Match(aa2, match2); 
      }
      for (int k=0;  k<numAlignedAtoms; k++) {
        ff=match1.x[k]-match2.x[k];
        rms += ff*ff;
        ff=match1.y[k]-match2.y[k];
        rms += ff*ff;
        ff=match1.z[k]-match2.z[k];
        rms += ff*ff;
        //printf("%g,%g,%g <-> %g,%g,%g\n", match1.x[k],match1.y[k],match1.z[k], match2.x[k],match2.y[k],match2.z[k]  ); 
      }
      rmsd = (float) sqrt(rms/numAlignedAtoms); 
    }
  }
  return rtn;

}  // overlap()

int OverlapAlign::rms_dist(COORDINATES *match1, COORDINATES *match2, int n)
{
 int     k,ii,jj,kk;
 double  kabsch[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double  kab2[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double  b[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double  rot[3][3]= {{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double  eigval[3]= {0.0,0.0,0.0};
 double  eigvec[3][3];
 double  check[3][3]={{0.0,0.0,0.0},{0.0,0.0,0.0},{0.0,0.0,0.0}};
 double  norm,xcm_prt,ycm_prt,zcm_prt,xcm_nat,ycm_nat,zcm_nat;
 int     find_eigen;
 double  mysmall = 1.e+20 ;
 int     ismall;
 double  det; 
 int     set = 0;      			/* Take C_alpha  posoition */
 double  epsilon = 1.e-6 ;

 if ( n <= 2 )
   return 0;          // unable to overlap

 xcm_nat=0.0;  ycm_nat=0.0;  zcm_nat=0.0;
 xcm_prt=0.0;  ycm_prt=0.0;  zcm_prt=0.0;
 
 for (k=0;  k<n; k++) {   
        xcm_nat += match2->x[k];
        ycm_nat += match2->y[k];
        zcm_nat += match2->z[k];
 
        xcm_prt += match1->x[k];
        ycm_prt += match1->y[k];
        zcm_prt += match1->z[k];
 }

  norm=1.0/(double)n;
  xcm_nat=xcm_nat*norm;
  ycm_nat=ycm_nat*norm;
  zcm_nat=zcm_nat*norm;
  xcm_prt=xcm_prt*norm;
  ycm_prt=ycm_prt*norm;
  zcm_prt=zcm_prt*norm;

 /*---------      subtract geometric center of native */
  for ( k=0; k< n; k++){
      match2->x[k] -= xcm_nat;
      match2->y[k] -= ycm_nat;
      match2->z[k] -= zcm_nat;
  }

  /*---------      subtract geometric center of prot */
  for ( k=0; k<n; k++){
      match1->x[k] -= xcm_prt;
      match1->y[k] -= ycm_prt;
      match1->z[k] -= zcm_prt;
  }

 /*Construct the Kabsch matrix  R  ------- */

 for (k=0;  k<n; k++) {
      kabsch[0][0] +=  match2->x[k] * match1->x[k];
	kabsch[1][0] +=  match2->x[k] * match1->y[k];;
	kabsch[2][0] +=  match2->x[k] * match1->z[k];;

	kabsch[0][1] +=  match2->y[k] * match1->x[k] ;
	kabsch[1][1] +=  match2->y[k] * match1->y[k];
	kabsch[2][1] +=  match2->y[k] * match1->z[k];
		 
	kabsch[0][2] +=  match2->z[k] * match1->x[k];
	kabsch[1][2] +=  match2->z[k] * match1->y[k];
	kabsch[2][2] +=  match2->z[k] * match1->z[k];
  }
 /*----------multiply kabsch by its transpose: R^T R  ---------------- */
 for (ii=0;ii<3;ii++) 
   for (jj=0;jj<3;jj++) 
     for (kk=0;kk<3;kk++)
       kab2[ii][jj] += kabsch[kk][ii] * kabsch[kk][jj];

 for (ii=0;ii<3; ii++)
   for (jj=0;jj<3; jj++)
     eigvec[ii][jj] = kab2[ii][jj];

 /*get normailzed eigenvectors stored as coloumns */
 find_eigen = diag_kabsch(eigvec,eigval);  
 if (find_eigen == NO) 
	 return 0 ; //	(-1.0);
 
 /* Find smallest eigenvalue */
 for (ii=0; ii<3;ii++){
   if (eigval[ii] < mysmall){
      ismall = ii;
      mysmall  = eigval[ii];
   }
 }
        
 /*compute normalized  b vectors stored as coloumns  b=R^T.a.norm--------- */
  for (jj=0;jj<3;jj++) {
     norm=1.0/sqrt(eigval[jj]);   
      for (ii=0;ii<3;ii++) {
       for (kk=0;kk<3;kk++)
	 b[ii][jj] += kabsch[ii][kk] * eigvec[kk][jj] * norm;
     }
  }

  // fix for all atoms in a single plane
  if ( eigval[0] <= RMS_ZERO || eigval[0] <= RMS_ZERO )
    return 0;

  if (eigval[2] <= RMS_ZERO)
  {
    eigval[2] = 0.0;
    eigvec[0][2]=eigvec[1][0]*eigvec[2][1] - eigvec[1][1]*eigvec[2][0];
    eigvec[1][2]=eigvec[2][0]*eigvec[0][1] - eigvec[2][1]*eigvec[0][0];
    eigvec[2][2]=eigvec[0][0]*eigvec[1][1] - eigvec[0][1]*eigvec[1][0];

    b[0][2]=b[1][0]*b[2][1] - b[1][1]*b[2][0];
    b[1][2]=b[2][0]*b[0][1] - b[2][1]*b[0][0];
    b[2][2]=b[0][0]*b[1][1] - b[0][1]*b[1][0];
  }
  
  /*Find if det U is negative that is if det R is negative */
  det = kabsch[0][0]*kabsch[1][1]*kabsch[2][2] - kabsch[0][0]*kabsch[1][2]*kabsch[2][1]
      - kabsch[0][1]*kabsch[1][0]*kabsch[2][2] + kabsch[0][1]*kabsch[1][2]*kabsch[2][0] 
      + kabsch[0][2]*kabsch[1][0]*kabsch[2][1] - kabsch[0][2]*kabsch[1][1]*kabsch[2][0];
  
  if(det < 0) {
    for (ii=0;ii<3; ii++)
        b[ii][ismall] *=-1.0 ;
  }

 /* rotation matrix  U = bXa^T where R^TR=I------ */
 for (ii=0;ii<3;ii++)     
   for (jj=0;jj<3;jj++) 
     for (kk=0;kk<3;kk++)
       rot[ii][jj] += b[ii][kk] * eigvec[jj][kk];

 // Output:  copy over translate and rotation matrices
 double *overlapmatrix= (double*)OverlapMatrix;

 *(overlapmatrix++) = -xcm_nat;             // translate for all atoms in aa2
 *(overlapmatrix++) = -ycm_nat;
 *(overlapmatrix++) = -zcm_nat;

 *(overlapmatrix++) = xcm_prt;              // don't translate aa1, move aa2 instead after rotate
 *(overlapmatrix++) = ycm_prt;
 *(overlapmatrix++) = zcm_prt;

 memcpy(overlapmatrix, rot, sizeof(double)*9);   // rotate for all atoms in aa2 , row major
 
 return  n;   

}	//  rms_dist()

/* =====================================================================
  The next functions were taken from Numerical Recipes 
  flowers for the NR guys 
  note small modifications with respect to original code 
   ======================================================================
*/

int OverlapAlign::diag_kabsch(double kab[3][3], double *eigval){
  int /* i,*/ j,k, /*kk,l,ll, */ nrot;
 double a[3][3],b[3][3]  /*,c[3][3]*/; /* let c be a buffer for shift */
  double *d,**v,**e;
  int    find_eigen;

  for (j=0;j<3;j++) 
    for (k=0;k<3;k++) 
	a[k][j]=kab[k][j];

  /* allocate such that indices start from 1 instead of 0 */
  d=vector(1,3);
  e=convert_matrix(&a[0][0],1,3,1,3);
  v=convert_matrix(&b[0][0],1,3,1,3);
  find_eigen = jacobi(e,3,d,v,&nrot);
  if (find_eigen == NO ) return(find_eigen);
  eigsrt(d,v,3);
  
  /* destroy now the original matrix kab */
  for (j=1;j<=3;j++) *(eigval+j-1)=d[j];
  for (j=1;j<=3;j++) 
    for (k=1;k<=3;k++) 
      kab[k-1][j-1]=v[k][j];

  free_convert_matrix(e,1,3,1,3);
  free_convert_matrix(v,1,3,1,3);
  free_vector(d,1,3);
  return(1);
}

#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);\
	a[k][l]=h+s*(g-h*tau);


/* ========================================================= 
   diagonalize matrix using Jacobi method 
   =========================================================
*/
int  OverlapAlign::jacobi(double **a, int n, double d[], double **v, int *nrot)
{
  int j,iq,ip,i;
  double tresh,theta,tau,t,sm,s,h,g,c,*b,*z;
  int find_eigen = YES;

  b=vector(1,n);
  z=vector(1,n);
  for (ip=1;ip<=n;ip++) 
    {
      for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
      v[ip][ip]=1.0;
    }
  for (ip=1;ip<=n;ip++) 
    {
      b[ip]=d[ip]=a[ip][ip];
      z[ip]=0.0;
    }
  *nrot=0;
  for (i=1;i<=50;i++) 
    {
      sm=0.0;
      for (ip=1;ip<=n-1;ip++) 
	{
	  for (iq=ip+1;iq<=n;iq++)
	    sm += fabs(a[ip][iq]);
	}
      if (sm == 0.0) 
	{
	  free_vector(z,1,n);
	  free_vector(b,1,n);
	  return (find_eigen);
	}
      if (i < 4) tresh=0.2*sm/(n*n);
      else tresh=0.0;
      for (ip=1;ip<=n-1;ip++) 
	{
	  for (iq=ip+1;iq<=n;iq++) 
	    {
	      g=100.0*fabs(a[ip][iq]);
	      if (i > 4 && (double)(fabs(d[ip])+g) == (double)fabs(d[ip])
		  && (double)(fabs(d[iq])+g) == (double)fabs(d[iq]))
		a[ip][iq]=0.0;
	      else if (fabs(a[ip][iq]) > tresh) 
		{
		  h=d[iq]-d[ip];
		  if ((double)(fabs(h)+g) == (double)fabs(h))
		    t=(a[ip][iq])/h;
		  else 
		    {
		      theta=0.5*h/(a[ip][iq]);
		      t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
		      if (theta < 0.0) t = -t;
		    }
		  c=1.0/sqrt(1+t*t);
		  s=t*c;
		  tau=s/(1.0+c);
		  h=t*a[ip][iq];
		  z[ip] -= h;
		  z[iq] += h;
		  d[ip] -= h;
		  d[iq] += h;
		  a[ip][iq]=0.0;
		  for (j=1;j<=ip-1;j++) 
		    {
		      ROTATE(a,j,ip,j,iq)
		    }
		  for (j=ip+1;j<=iq-1;j++) 
		    {
		      ROTATE(a,ip,j,j,iq)
		    }
		  for (j=iq+1;j<=n;j++) 
		    {
		      ROTATE(a,ip,j,iq,j)
		    }
		  for (j=1;j<=n;j++) 
		    {
		      ROTATE(v,j,ip,j,iq)
		    }
		  ++(*nrot);
		}
	    }
	}
      for (ip=1;ip<=n;ip++) 
	{
	  b[ip] += z[ip];
	  d[ip]=b[ip];
	  z[ip]=0.0;
	}
    }
  find_eigen = NO;
  return(find_eigen);
}

#undef ROTATE

/* ===================================================
  sort eigenvectors in descending eigenvalue order 
   ===================================================
*/
void OverlapAlign::eigsrt(double d[], double **v, int n)
{
  int k,j,i;
  double p;

  for (i=1;i<n;i++) 
    {
      p=d[k=i];
      for (j=i+1;j<=n;j++)
	if (d[j] >= p) p=d[k=j];
      if (k != i) 
	{
	  d[k]=d[i];
	  d[i]=p;
	  for (j=1;j<=n;j++) 
	    {
	      p=v[j][i];
	      v[j][i]=v[j][k];
	      v[j][k]=p;
	    }
	}
    }
}

#define NR_END 1


/*===================================================
  allocate a float vector with subscript range v[nl..nh] 
  =====================================================
*/
double *OverlapAlign::vector(long nl, long nh)
{
  double *v;

  v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
  // ASSERT(v ,printf("allocation failure in vector()\n"));
  return v-nl+NR_END;
}


/* ========================================================
  free a float vector allocated with vector() 
   ========================================================
*/
void OverlapAlign::free_vector(double *v, long nl, long nh)
{
  free((char* ) (v+nl-NR_END));
}

/* ======================================================================
   ======================================================================
*/
double **OverlapAlign::convert_matrix(double *a, long nrl,long nrh, long ncl, long nch)
{
  long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
  double **m;

  /* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
     declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
     and ncol=nch-ncl+1. The routine should be called with the address
     &a[0][0] as the first argument. */

  /* allocate pointers to rows */
  m=(double **) malloc((unsigned int) ((nrow+NR_END)*sizeof(double*)));
  // ASSERT(m,printf("allocation failure in convert_matrix()"));
  m += NR_END;
  m -= nrl;

  /* set pointers to rows */
  m[nrl]=a-ncl;
  for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
  /* return pointer to array of pointers to rows */
  return m;
}


/* =============================================================
   free a matrix allocated by convert_matrix() 
   =============================================================
*/

void OverlapAlign::free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch)
{
  free((char* ) (b+nrl-NR_END));
}

#undef NR_END

}
