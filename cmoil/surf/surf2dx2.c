#include <stdio.h>

#define NTRI 6466
#define NPOS 3*NTRI

#define RSQR(x,y,z)  ((x)*(x)+(y)*(y)+(z)*(z))

float xmin,xmax,ymin,ymax,zmin,zmax;
float *pos,*nrm;
int refinter, hashdim, *hashtable, *hash_index_list;
float tol2;

main()

{

FILE *fin, *fout, *fout1, *fout2;

int atom_number,i,ncon,ind1,ind2,ind3;

float x,y,z,nx,ny,nz;

int Hash(float x, float y, float z, float nx, float ny, float nz);

fin = fopen("CA.sys.tri","r");
fout  = fopen("CA.sys.pos","w");
fout1 = fopen("CA.sys.con","w");
fout2 = fopen("CA.sys.nor","w");

/* allocate and initialize */

pos = (float *) malloc(3*sizeof(float)*NPOS);
nrm = (float *) malloc(3*sizeof(float)*NPOS);

hashdim = 50; tol2 = 1.0e-10;

hashtable = (int *) malloc(sizeof(int)*hashdim*hashdim*hashdim);
hash_index_list = (int *) malloc(sizeof(int)*NPOS);

for (i=0; i<hashdim*hashdim*hashdim; ++i) { hashtable[i] = -1; }
for (i=0; i<NPOS; ++i) { hash_index_list[i] = -1; }

refinter = 0;

/* bounding box */

xmin = 10.0;    xmax = 48.7;
ymin = -12.4;   ymax = 56.3;
zmin = 8.8;   zmax = 66.0;

/* read in verts and normals and hash  */

for (i=0; i<NTRI; ++i) {

fscanf(fin,"%d",&atom_number);

fscanf(fin,"%f %f %f %f %f %f",&x,&y,&z,&nx,&ny,&nz);
ind1 = Hash(x,y,z,nx,ny,nz);

fscanf(fin,"%f %f %f %f %f %f",&x,&y,&z,&nx,&ny,&nz);
ind2 = Hash(x,y,z,nx,ny,nz);

fscanf(fin,"%f %f %f %f %f %f",&x,&y,&z,&nx,&ny,&nz);
ind3 = Hash(x,y,z,nx,ny,nz);

fprintf(fout1,"%d %d %d\n",ind1,ind2,ind3);

}

/* write out the final non-redundant positions */

fprintf(stderr,"%d positions\n",refinter);

for (i=0; i<refinter; ++i) {

 fprintf(fout,"%f %f %f\n",pos[3*i],pos[3*i+1],pos[3*i+2]);
 fprintf(fout2,"%f %f %f\n",nrm[3*i],nrm[3*i+1],nrm[3*i+2]);

}

}

/* return index of position using hash table */

int Hash(float x, float y, float z, float nx, float ny, float nz)

{

/* uses the following globals: tol2, hashdim */

float dist2;
int i,j,k,index,old_index;

/* compute integer location of point within cube */

i = (int) ((hashdim-1)*(x-xmin)/(xmax-xmin));
j = (int) ((hashdim-1)*(y-ymin)/(ymax-ymin));
k = (int) ((hashdim-1)*(z-zmin)/(zmax-zmin));

/* look up index of first position in hashtable */

index = hashtable[hashdim*hashdim*k+hashdim*j+i]; 

if (index >= 0) {  

  dist2 = RSQR(x-pos[3*index],y-pos[3*index+1],z-pos[3*index+2]);
  
  if (dist2 <= tol2) { return(index);}  

  while ((dist2 > tol2)&&(index >=0)) { 

        old_index = index;
        index = hash_index_list[index];

        if (index >= 0) {               
           dist2 = RSQR(x-pos[3*index],y-pos[3*index+1],z-pos[3*index+2]);
	 }
        else { 
           pos[3*refinter] = x;  pos[3*refinter+1] = y;  pos[3*refinter+2] = z; 
           nrm[3*refinter] = nx; nrm[3*refinter+1] = ny; nrm[3*refinter+2] = nz; 
           hash_index_list[old_index] = refinter;
           ++refinter;
           return(refinter-1);
	   }
	   }

  return(index);  

} else {

   pos[3*refinter] = x;  pos[3*refinter+1] = y;  pos[3*refinter+2] = z; 
   nrm[3*refinter] = nx; nrm[3*refinter+1] = ny; nrm[3*refinter+2] = nz; 
   hashtable[hashdim*hashdim*k+hashdim*j+i]  = refinter;
   ++refinter;

   return(refinter-1);

 }

}
 

