#include <stdio.h>

main()

{

FILE *fin, *fout, *fout1, *fout2;

int atom_number,i,ncon;

float x,y,z,nx,ny,nz;

fin = fopen("ACh.sys.tri","r");
fout  = fopen("ACh.pos","w");
fout1 = fopen("ACh.con","w");
fout2 = fopen("ACh.nor","w");

ncon = 0;

for (i=0; i<2078; ++i) {

fscanf(fin,"%d",&atom_number);

fscanf(fin,"%f %f %f %f %f %f",&x,&y,&z,&nx,&ny,&nz);
fprintf(fout,"%f %f %f\n",x,y,z);
fprintf(fout2,"%f %f %f\n",nx,ny,nz);

fscanf(fin,"%f %f %f %f %f %f",&x,&y,&z,&nx,&ny,&nz);
fprintf(fout,"%f %f %f\n",x,y,z);
fprintf(fout2,"%f %f %f\n",nx,ny,nz);

fscanf(fin,"%f %f %f %f %f %f",&x,&y,&z,&nx,&ny,&nz);
fprintf(fout,"%f %f %f\n",x,y,z);
fprintf(fout2,"%f %f %f\n",nx,ny,nz);

fprintf(fout1,"%d %d %d\n",ncon,ncon+1,ncon+2);

ncon += 3;

}

}

