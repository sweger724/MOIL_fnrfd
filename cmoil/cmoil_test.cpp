#include <cmoil_dib.h>
#include <LIST>
#include <VECTOR>
#include <Magick++.h>

void main (int argc, char **argv) 
{
  char imageTmpName[256];
  BYTE bytes8[8];
   
  if (argc == 1) return;

  FILE *in=fopen(argv[1], "r+b");
  if ( in == NULL )
    return;

  fread(bytes8, 1, 8, in);
  fread(bytes8, 1, 3, in);

  printf("%x\n", bytes8[2]);

}
