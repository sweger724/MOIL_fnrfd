#include <stdio.h>

const unsigned short int ALIGNSIZE=8;	

void SwapBite(char *s, int aSize)
{
  char t;

  aSize--; 		// change to index
  for (int i=0; i<= aSize/2; i++) 
  {
    t = s[i];
    s[i] = s[aSize-i];
    s[aSize-i] = t;
  }
}

void main(int argc, char **argv)
{
  if ( argc <= 1 ) 
  {
    fprintf(stderr, "\nswap8p.exe converts coordinate file in path format between big-endian and little-endian systems.\n");
    fprintf(stderr, "\nUsage: swap8p.exe [input pth|dcd filename] {[output bin filename]}\n\n");
    fprintf(stderr, "       {[output bin filename]} is optional, default is swap.out in working directory\n");   
    return;
  }

  char *fileName=argv[1];
  char *outputName="swap.out";
  if ( argc >= 3 )
  {
     outputName=argv[2];
  }
    
  FILE *crdfp = fopen(fileName, "rb");
  FILE *outfp = fopen(outputName, "wb");

  if ( crdfp==NULL || outfp==NULL )
  {
    fprintf(stderr, "Unable to open input/output file: %s, %s\n", fileName, outputName);
    return;
  }

  char dummy[ALIGNSIZE];
  for (int i=1; i<=3; i++)         // read three integers
  {
    if ( fread(&dummy,ALIGNSIZE/2,1,crdfp)==1 )
    {
      SwapBite(dummy, ALIGNSIZE/2);
      fwrite(&dummy,ALIGNSIZE/2,1,outfp);
    }
  }

  while ( !feof(crdfp) )          // read doubles
  {
    if ( fread(&dummy,ALIGNSIZE,1,crdfp)==1 ) {
      SwapBite(dummy, ALIGNSIZE);
      fwrite(&dummy,ALIGNSIZE,1,outfp);
    }
  }

  if ( fread(&dummy,ALIGNSIZE/2,1,crdfp)==1 )         // read another 4 bytes
  {
    SwapBite(dummy, ALIGNSIZE/2);
    fwrite(&dummy,ALIGNSIZE/2,1,outfp);
  }

  fclose(crdfp);
  fclose(outfp);

}
