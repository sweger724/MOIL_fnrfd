#define	EXTERN	
#include "surf.h"
#undef EXTERN
#define byte BYTE
#define Vector VECTOR
#include <dx/dx.h>


#define RSQR(x,y,z)  ((x)*(x)+(y)*(y)+(z)*(z))

static Error DoSurf(Object o, float probe, float gridres, int colorflag);

int refinter, hashdim, tri_count, nohash, *hashtable, *hash_index_list, natoms;
int colorflag,xmax,xmin,ymax,ymin,zmax,zmin;

float *pos,*poss,*nrml, *vdw_radii,margin,tol2;
int *conn,*ndata;
float *colors,*ncolors; 

m_Surf(Object *in, Object *out)
{
    Object o = NULL;
    float gridres, probe;
    int maxpolys;

    if (!in[1]) { probe = 1.4;}             /* radius of probe sphere */
    else { DXExtractFloat(in[1],&probe); }

    if (!in[2]) { gridres = 1.2;}           /* maximum triangle edge */
    else { DXExtractFloat(in[2],&gridres); }

    if (!in[3]) { maxpolys = 10;}           /* expected polys in thousands */
    else { DXExtractInteger(in[3],&maxpolys); }

    if (!in[4]) { hashdim =  50;} /* size of hash grid */
    else { DXExtractInteger(in[4],&hashdim); }

    if (!in[5]) { nohash =  0;} /* hash vertices or not */
    else { DXExtractInteger(in[5],&nohash); }

    if (!in[6]) { colorflag =  0;} /* skip hashing */
    else { DXExtractInteger(in[6],&colorflag); }

    tol2 = 1.0e-20; /* closeness criterion for hashing */
   
    tol2 *= tol2;

    margin = 1.0;

 /* set global variables for Varshney's code */

   Checks_On = FALSE;
   Write_Option = 2;
   Num_atoms = 0;
   Probe_radius = probe;
   Max_Tess_Len = gridres;
   Max_Gp_Polys = 1000*maxpolys;

 /* copy the input structure */

    if (!in[0]) DXErrorGoto(ERROR_BAD_PARAMETER, "missing object");

    o = DXCopy(in[0], COPY_STRUCTURE);

    if (!o) goto error;
  
    /* Do the recursive traversal.  */

    if (!DoSurf(o,probe,gridres,colorflag)) goto error;

    /* successful return */

    out[0] = o;

    return OK;

error:
    DXDelete(o);
    return ERROR;
  }

static Error DoSurf(Object o,float probe, float gridres,int colorflag) {
    Object oo;

    Array a,b,c,d;
    Array na,nb,nc,nd,ne;

    float hx,hy,hz,x,y,z;
    float *corners;
    int i,j,k,l,n,npts;

    int ind0,ind1,ind2;

    int ip;
    float cxmin,cxmax,cymin,cymax,czmin,czmax;
    float dist,distmin,distmina,distminb,tmp;
    float res;

    float rtmp,gtmp,btmp;

    int s,il,jl,kl,iv,jv,kv,ic,jc,kc;
    int index,corner,iatom,iatoma,iatomb;

    int isubdiv,itmp;

    /* determine the class of the object */

    switch (DXGetObjectClass(o)) {
    
    case CLASS_FIELD:

	/* 
         * Extract component arrays from the input field
         */

	a = (Array) DXGetComponentValue((Field)o, "positions");

	if (!a)
	    DXErrorReturn(ERROR_MISSING_DATA, "input has no positions");
	if (!DXTypeCheck(a, TYPE_FLOAT, CATEGORY_REAL, 1, 3))
	    DXErrorReturn(ERROR_BAD_TYPE, "positions are not 3D floating point");
	if (!DXGetArrayInfo(a, &natoms, NULL, NULL, NULL, NULL))
	    return ERROR;
	pos = (float *) DXGetArrayData(a);
	if (!pos)
	    return ERROR;

	b = (Array) DXGetComponentValue((Field)o, "vdw_radii");

	if (!b)
	    DXErrorReturn(ERROR_MISSING_DATA, "input has no vdw_radii");
	if (!DXTypeCheck(b, TYPE_FLOAT, CATEGORY_REAL, 0, 0))
	    DXErrorReturn(ERROR_BAD_TYPE, "problem with vdw_radii array");
	if (!DXGetArrayInfo(b, &natoms, NULL, NULL, NULL, NULL))
	    return ERROR;
	vdw_radii = (float *) DXGetArrayData(b);
	if (!vdw_radii)
	    return ERROR;

	c = (Array) DXGetComponentValue((Field)o, "box");

	if (!c)
	    DXErrorReturn(ERROR_MISSING_DATA, "input has no box");
	if (!DXTypeCheck(c, TYPE_FLOAT, CATEGORY_REAL, 1, 3))
	    DXErrorReturn(ERROR_BAD_TYPE, "problem with box array");
	if (!DXGetArrayInfo(c, &i, NULL, NULL, NULL, NULL))
	    return ERROR;
	corners = (float *) DXGetArrayData(c);
	if (!corners) return ERROR;

        if (colorflag) {

           d = (Array) DXGetComponentValue((Field)o, "colors");

    	   if (!d)
	       DXErrorReturn(ERROR_MISSING_DATA, "input has no colors");
	   if (!DXTypeCheck(d, TYPE_FLOAT, CATEGORY_REAL, 1, 3))
	       DXErrorReturn(ERROR_BAD_TYPE, "problem with colors array");
	   colors = (float *) DXGetArrayData(d);
	   if (!colors) return ERROR;

	 }

        /* will create new positions, connections, colors */

	DXChangedComponentStructure((Field)o, "connections");
	DXChangedComponentStructure((Field)o, "positions");
        DXChangedComponentStructure((Field)o, "normals");

	if (colorflag) DXChangedComponentStructure((Field)o, "colors");

/* set external variables and call Varshney's code */

   Num_atoms = natoms;

/* transfer atoms and data to Varshney's interal arrays */

   ainput(natoms,pos,vdw_radii);   

   DXMessage("Constructing solvent-accessible surface ..");

   /* start timing the molecular surface computations here */
   START

   /* initialize and compute the molecular surface */

   init_and_compute(); 

   /* stop timing the molecular surface computations here */
   STOP

   if (Write_Option == 1)
     DXMessage("Surface construction + writing time %4.2f seconds", et);
   else
     DXMessage("Surface construction time %4.2f seconds", et);

   DXMessage(">>> %d polygons found <<<",Num_polys);

   DXMessage("---- Patches ----");
   DXMessage("%d case I,II toriodal",torus_count_I);
   DXMessage("%d case I,II convex",convex_count_I);
   DXMessage("%d case I,II concave",concave_count_I);
   DXMessage("%d case III toroidal",torus_count_III);
   DXMessage("%d case III convex",convex_count_III);
   DXMessage("%d case IV spherical",sphere_count_IV);
   DXMessage(" ");
   DXMessage("%d patches total",torus_count_I+convex_count_I+
              concave_count_I+torus_count_III+convex_count_III+sphere_count_IV);
   DXMessage("--------------");

        /* create new components */

        /* this allocates space but does not set number of items yet */
        /* DXTrim will remove unneeded elements  later */

        na =   DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3);

               DXAllocateArray(na,3*Num_polys);  

        poss = (float *) DXGetArrayData(na);


        nb =   DXNewArray(TYPE_INT,CATEGORY_REAL,1,3);

               DXAllocateArray(nb,Num_polys); 

        conn = (int *) DXGetArrayData(nb);


        nc =  DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3);

               DXAllocateArray(nc,3*Num_polys); 

        nrml = (float *) DXGetArrayData(nc);

        if (colorflag) {

           nd =  DXNewArray(TYPE_FLOAT,CATEGORY_REAL,1,3);
                 DXAllocateArray(nd,3*Num_polys); 

           ncolors = (float *) DXGetArrayData(nd);
	 }

           ne =  DXNewArray(TYPE_INT,CATEGORY_REAL,1,1);
                 DXAllocateArray(ne,3*Num_polys); 

           ndata = (int *) DXGetArrayData(ne);
	 

        refinter = 0; /* initialize new position count */ 

/* extract actual box dimensions of molecule */

 xmax = corners[12];  xmin = corners[0];
 ymax = corners[7 ];  ymin = corners[1];
 zmax = corners[5 ];  zmin = corners[2];

/* add enough length to ensure that all atomic     */
/* radii fit inside the new box including probe    */


for (i=0; i<natoms; ++i) {
    if (vdw_radii[i] > margin) margin = vdw_radii[i]+.1;
}

xmax += margin+probe;   xmin -= margin+probe;
ymax += margin+probe;   ymin -= margin+probe;
zmax += margin+probe;   zmin -= margin+probe;

/* hash vertices to eliminate redundancies */

if (nohash) { transfer_vertices(colorflag);}

else { DXMessage("Hashing (removing redundant vertices)");
       START
       hash_vertices();
       STOP
       DXMessage("finished hashing in %4.2f sec",et); }

/* free pointers used by Varshney's code */

 FREEMAT(verts, Max_Gp_Polys, 3); 

 free(atoms);  

/* transfer colors */

if (colorflag) { for (i=0; i<refinter; ++i) {

   ncolors[3*i  ] = colors[3*ndata[i]  ];
   ncolors[3*i+1] = colors[3*ndata[i]+1];
   ncolors[3*i+2] = colors[3*ndata[i]+2];

 }
	       }

/* set actual number of positions  for trim */

  DXAddArrayData(na,0,refinter,NULL);
  DXAddArrayData(nb,0,tri_count,NULL);
  DXAddArrayData(nc,0,refinter,NULL);
  if (colorflag) {DXAddArrayData(nd,0,refinter,NULL);}
  DXAddArrayData(ne,0,refinter,NULL);

  DXSetComponentValue((Field)o,"positions",(Object)na);
  DXSetComponentValue((Field)o,"connections",(Object)nb);
  DXSetComponentValue((Field)o,"normals",(Object)nc);
  if (colorflag) {DXSetComponentValue((Field)o,"colors",(Object)nd);}
  DXSetComponentValue((Field)o,"data",(Object)ne);

  DXSetComponentAttribute((Field)o,"connections","ref",
                       (Object)DXNewString("positions"));

  DXSetComponentAttribute((Field)o,"connections","element type",
                       (Object)DXNewString("triangles"));

  DXSetComponentAttribute((Field)o,"normals","dep",
                       (Object)DXNewString("positions"));

  if (colorflag) {
  DXSetComponentAttribute((Field)o,"colors","dep",
                       (Object)DXNewString("positions"));

  DXSetComponentAttribute((Field)o,"data","dep",
                       (Object)DXNewString("positions"));
}

   /* finalize field */

	if (!DXEndField((Field)o))
	    return ERROR;
	break;

    case CLASS_GROUP:

	/* recursively traverse groups */
	for (i=0; oo=DXGetEnumeratedMember((Group)o, i, NULL); i++)
	    if (!DoSurf(oo,probe,gridres,colorflag))
		return ERROR;
	break;

    default:
	DXErrorReturn(ERROR_BAD_PARAMETER, "object must be a group or field");

  }
  
    /* successful return */
    return OK;
   }

/* allocate transfer atoms to Varshney's struct */

ainput(int natoms, float *pos, float *vdw_radii)

{
  int i;

    ALLOCN(atoms, Gp_Atom, natoms);

    for(i = 0; i < natoms; i++) {
       atoms[i].radius = vdw_radii[i];
       atoms[i].center[X] = pos[3*i];
       atoms[i].center[Y] = pos[3*i+1];
       atoms[i].center[Z] = pos[3*i+2];
       atoms[i].type = 0;
    }

}

/* direct transfer of vertices without hashing */

transfer_vertices(colorflag)
int colorflag;
{

int i,dgen_count;
float tmp0,tmp1,tmp2;

tri_count = 0;
dgen_count = 0;

for (i=0; i<Num_polys; ++i) {


/* check for degenerate triangles */

  tmp0 = (verts[i][0].Coord[X] - verts[i][1].Coord[X])*(verts[i][0].Coord[X] - verts[i][1].Coord[X])+
         (verts[i][0].Coord[Y] - verts[i][1].Coord[Y])*(verts[i][0].Coord[Y] - verts[i][1].Coord[Y])+
         (verts[i][0].Coord[Z] - verts[i][1].Coord[Z])*(verts[i][0].Coord[Z] - verts[i][1].Coord[Z]);

  tmp1 = (verts[i][0].Coord[X] - verts[i][2].Coord[X])*(verts[i][0].Coord[X] - verts[i][2].Coord[X])+
         (verts[i][0].Coord[Y] - verts[i][2].Coord[Y])*(verts[i][0].Coord[Y] - verts[i][2].Coord[Y])+
         (verts[i][0].Coord[Z] - verts[i][2].Coord[Z])*(verts[i][0].Coord[Z] - verts[i][2].Coord[Z]);

  tmp2 = (verts[i][2].Coord[X] - verts[i][1].Coord[X])*(verts[i][2].Coord[X] - verts[i][1].Coord[X])+
         (verts[i][2].Coord[Y] - verts[i][1].Coord[Y])*(verts[i][2].Coord[Y] - verts[i][1].Coord[Y])+
         (verts[i][2].Coord[Z] - verts[i][1].Coord[Z])*(verts[i][2].Coord[Z] - verts[i][1].Coord[Z]);

if ((tmp0 > 1.0e-10) && (tmp1 > 1.0e-10) && (tmp2 > 1.0e-10)) {

  poss[9*tri_count  ] = verts[i][0].Coord[X];
  poss[9*tri_count+1] = verts[i][0].Coord[Y];
  poss[9*tri_count+2] = verts[i][0].Coord[Z];

  poss[9*tri_count+3] = verts[i][1].Coord[X];
  poss[9*tri_count+4] = verts[i][1].Coord[Y];
  poss[9*tri_count+5] = verts[i][1].Coord[Z];

  poss[9*tri_count+6] = verts[i][2].Coord[X];
  poss[9*tri_count+7] = verts[i][2].Coord[Y];
  poss[9*tri_count+8] = verts[i][2].Coord[Z];

  nrml[9*tri_count  ] = verts[i][0].Normal[X];
  nrml[9*tri_count+1] = verts[i][0].Normal[Y];
  nrml[9*tri_count+2] = verts[i][0].Normal[Z];

  nrml[9*tri_count+3] = verts[i][1].Normal[X];
  nrml[9*tri_count+4] = verts[i][1].Normal[Y];
  nrml[9*tri_count+5] = verts[i][1].Normal[Z];

  nrml[9*tri_count+6] = verts[i][2].Normal[X];
  nrml[9*tri_count+7] = verts[i][2].Normal[Y];
  nrml[9*tri_count+8] = verts[i][2].Normal[Z];

  ndata[3*tri_count  ] = atom_type[i];
  ndata[3*tri_count+1] = atom_type[i];
  ndata[3*tri_count+2] = atom_type[i];

  ++tri_count;

} else { ++dgen_count;
/*
  DXMessage("0:%f %f %f\n",verts[i][0].Coord[X],verts[i][0].Coord[Y],verts[i][0].Coord[Z]);
  DXMessage("1:%f %f %f\n",verts[i][1].Coord[X],verts[i][1].Coord[Y],verts[i][1].Coord[Z]);
  DXMessage("2:%f %f %f\n",verts[i][2].Coord[X],verts[i][2].Coord[Y],verts[i][2].Coord[Z]);
*/
       }

}

/* fill in connections component */

  for (i=0; i< tri_count; ++i) {
   conn[3*i] = 3*i; conn[3*i+1] = 3*i+1; conn[3*i+2] = 3*i+2;
}

if (dgen_count > 0) {DXMessage("%d degenerate triangles removed",dgen_count);}

refinter = 3*tri_count; /* number of vertices */

}

/* function to hash vertices */

hash_vertices()

{

int atom_number,i,ncon,ind1,ind2,ind3;
int dgen_count;
float x,y,z,nx,ny,nz;
float tmp0,tmp1,tmp2;

int Hash(float, float, float, float, float, float, int);

hashtable = (int *) DXAllocate(sizeof(int)*hashdim*hashdim*hashdim);
hash_index_list = (int *) DXAllocate(sizeof(int)*3*Num_polys);

for (i=0; i<hashdim*hashdim*hashdim; ++i) { hashtable[i] = -1; }
for (i=0; i<3*Num_polys; ++i) { hash_index_list[i] = -1; }

refinter = 0;
dgen_count = 0;
tri_count = 0;

for (i=0; i<Num_polys; ++i) {

/* check for degenerate triangles */

  tmp0 = (verts[i][0].Coord[X] - verts[i][1].Coord[X])*(verts[i][0].Coord[X] - verts[i][1].Coord[X])+
         (verts[i][0].Coord[Y] - verts[i][1].Coord[Y])*(verts[i][0].Coord[Y] - verts[i][1].Coord[Y])+
         (verts[i][0].Coord[Z] - verts[i][1].Coord[Z])*(verts[i][0].Coord[Z] - verts[i][1].Coord[Z]);

  tmp1 = (verts[i][0].Coord[X] - verts[i][2].Coord[X])*(verts[i][0].Coord[X] - verts[i][2].Coord[X])+
         (verts[i][0].Coord[Y] - verts[i][2].Coord[Y])*(verts[i][0].Coord[Y] - verts[i][2].Coord[Y])+
         (verts[i][0].Coord[Z] - verts[i][2].Coord[Z])*(verts[i][0].Coord[Z] - verts[i][2].Coord[Z]);

  tmp2 = (verts[i][2].Coord[X] - verts[i][1].Coord[X])*(verts[i][2].Coord[X] - verts[i][1].Coord[X])+
         (verts[i][2].Coord[Y] - verts[i][1].Coord[Y])*(verts[i][2].Coord[Y] - verts[i][1].Coord[Y])+
         (verts[i][2].Coord[Z] - verts[i][1].Coord[Z])*(verts[i][2].Coord[Z] - verts[i][1].Coord[Z]);

if ((tmp0 > 1.0e-10) && (tmp1 > 1.0e-10) && (tmp2 > 1.0e-10)) {

   ind1 = Hash(verts[i][0].Coord[X],verts[i][0].Coord[Y],verts[i][0].Coord[Z],
               verts[i][0].Normal[X],verts[i][0].Normal[Y],verts[i][0].Normal[Z],atom_type[i]);

   ind2 = Hash(verts[i][1].Coord[X],verts[i][1].Coord[Y],verts[i][1].Coord[Z],
               verts[i][1].Normal[X],verts[i][1].Normal[Y],verts[i][1].Normal[Z],atom_type[i]);

   ind3 = Hash(verts[i][2].Coord[X],verts[i][2].Coord[Y],verts[i][2].Coord[Z],
               verts[i][2].Normal[X],verts[i][2].Normal[Y],verts[i][2].Normal[Z],atom_type[i]);

   conn[3*tri_count] = ind1; conn[3*tri_count+1] = ind2; conn[3*tri_count+2] = ind3;
   ++tri_count;

 } else { ++dgen_count;

/*
  DXMessage("0:%f %f %f\n",verts[i][0].Coord[X],verts[i][0].Coord[Y],verts[i][0].Coord[Z]);
  DXMessage("1:%f %f %f\n",verts[i][1].Coord[X],verts[i][1].Coord[Y],verts[i][1].Coord[Z]);
  DXMessage("2:%f %f %f\n",verts[i][2].Coord[X],verts[i][2].Coord[Y],verts[i][2].Coord[Z]);
*/
	}

}

if (dgen_count > 0) {DXMessage("%d degenerate triangles removed",dgen_count);}

DXFree(hashtable);
DXFree(hash_index_list);

}

/* return index of position using hash table */

int Hash(float x, float y, float z, float nx, float ny, float nz, int aind)

{

/* uses the following globals: tol2, hashdim */

float dist2;
int i,j,k,index,old_index,loop;

/* compute integer location of point within cube */

i = (int) ((hashdim-1)*(x-xmin)/(xmax-xmin));
j = (int) ((hashdim-1)*(y-ymin)/(ymax-ymin));
k = (int) ((hashdim-1)*(z-zmin)/(zmax-zmin));

if (i > hashdim -1) DXMessage("i:%d",i);
if (j > hashdim -1) DXMessage("j:%d",j);
if (k > hashdim -1) DXMessage("k:%d",k);

/* look up index of first position in hashtable */

index = hashtable[hashdim*hashdim*k+hashdim*j+i]; 

if (index > Num_polys*3) {DXMessage("error"); exit(0);}

if (index >= 0) {  

  dist2 = RSQR(x-poss[3*index],y-poss[3*index+1],z-poss[3*index+2]);

  if (dist2 <= tol2) { 

     return(index);}  


  while ((dist2 > tol2)&&(index >=0)) { 

        old_index = index;

        index = hash_index_list[index];

        if (index >= 0) {               

           dist2 = RSQR(x-poss[3*index],y-poss[3*index+1],z-poss[3*index+2]);

	 }
        else { 
           poss[3*refinter] = x;  poss[3*refinter+1] = y;  poss[3*refinter+2] = z; 
           nrml[3*refinter] = nx; nrml[3*refinter+1] = ny; nrml[3*refinter+2] = nz; 
           ndata[refinter] = aind;

           hash_index_list[old_index] = refinter;
           ++refinter;  

           return(refinter-1);

	   }
	   }

  return(index);  

} else {

   poss[3*refinter] = x;  poss[3*refinter+1] = y;  poss[3*refinter+2] = z; 
   nrml[3*refinter] = nx; nrml[3*refinter+1] = ny; nrml[3*refinter+2] = nz; 
   ndata[refinter] = aind;

   hashtable[hashdim*hashdim*k+hashdim*j+i]  = refinter;

   ++refinter;

   return(refinter-1);

 }

}



