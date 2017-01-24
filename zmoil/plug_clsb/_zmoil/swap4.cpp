#include <stdio.h>

const unsigned short int ALIGNSIZE = 4;

void SwapBite( char *a )
{
	char t = *a;
	*a = a[3];
	a[3] = t;
	t = a[1];
	a[1] = a[2];
	a[2] = t;
}

void main( int argc, char **argv )
{
	if( argc <= 1 ) {
		fprintf( stderr,
				 "\nswap4.exe swaps bytes for each 4 bytes alignment boundaries.\nIt converts binary file format between big-endian and little-endian.\n" );
		fprintf( stderr, "\nUsage: swap4.exe [input pth|dcd filename] {[output bin filename]}\n\n" );
		fprintf( stderr, "       {[output bin filename]} is optional, default is swap.out in working directory\n" );
		return;
	}

	char *fileName = argv[1];
	char *outputName = "swap.out";
	if( argc >= 3 ) {
		outputName = argv[2];
	}

	FILE *crdfp = fopen( fileName, "rb" );
	FILE *outfp = fopen( outputName, "wb" );

	if( crdfp == NULL || outfp == NULL ) {
		fprintf( stderr, "Unable to open input/output file: %s, %s\n", fileName, outputName );
		return;
	}

	char dummy[ALIGNSIZE];
	while( !feof( crdfp ) ) {
		if( fread( &dummy, ALIGNSIZE, 1, crdfp ) == 1 ) {
			SwapBite( dummy );
			fwrite( &dummy, ALIGNSIZE, 1, outfp );
		}
	}
	fclose( crdfp );
	fclose( outfp );

}
