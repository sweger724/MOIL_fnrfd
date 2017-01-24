//*************************************************************************************************
//*  Filename:   cmoil_parm.cpp
//*
//*  Description: 
//*    Read in and process input file to get all parameters.
//*
//*  Modification History:
//*  
//*  Date           Developer   Description
//*  ------------   ----------- ---------------------------------------------------------------------
//*  Aug. 2000  ??      Initial Development
//*  July 2001  Baohua Wang Seperated from cmoil.cpp and reshape to support mutiple structures
//*  Sept 2007  Thomas Blom   Update for OSX/Linux
//*************************************************************************************************
#include "cmoil_parm.h"
#include "wingl.h"
#include "GL/gl.h"
#include "GL/glu.h"
#include <stdlib.h>
#include <ctype.h>

#include "cmoil_globals.h"

extern void trace( char *fmt, ... );

namespace CMOIL {

const float MOIL_COLORSCALE = ( float ) 0.004;	// map GUI color range [0,500] to gl range [0,2]

// default AminoAcid surface color in RGB color [0,255]
const int AA_DEF_COLOR[NUM_AA_TYPE][3] = { 126, 126, 126,	// grey    hydrophobic
	50, 126, 50,			// green   polar
	126, 30, 30,			// red     charged
	50, 50, 150,			// blue    negtively charged
	110, 60, 180,			//         unknown
	110, 60, 180			//         uni-color(single color)
};

InputParameters::InputParameters(  ) {
	fileName[0] = '\0';
	bmpName[0] = '\0';
	surfName[0] = '\0';
	connName[0] = 0;
	filetag = fNONE;
	sleepseconds = 0.0;
	structno = 1;			// tfb: this was uninitialized; default to 1 structure
	structend = -1;			// default to non-stop
	maxframes = 100;		// save max 100 frames for movie in default
	strcpy( imagetype, "gif" );	// default movie saving image format
	gifloops = 500;
	framedistance = 1;		// print all structures
	structpstart = 0;		// first structure to be printed:
	//  0:current structure; 1,..,i,..n: ith structure
	dispmode = STICK;		
	bondmode = 0;
	swapmode = 0;

	colorsettings = NULL;
	colorsettingsLen = 0;
	alignsettings = NULL;
	alignsettingsLen = 0;
	ribboncolor = NULL;
	ribboncolorLen = 0;
	stick = NULL;
	stickLen = 0;
	space = NULL;
	spaceLen = 0;
	ribbon = NULL;
	ribbonLen = 0;
	skip = NULL;
	skipLen = 0;
	skipsurf = NULL;
	skipsurfLen = 0;
	picksphere = NULL;
	picksphereLen = 0;
	pickradius = 60;		// default to 60A
	rotateDegree = ( float ) ( 5.0 * AngleToRadius );	// 5.0 degree for each button click by default
	subsurf = NULL;
	subsurfLen = 0;
	transparent = 255;		// 
	uniqueSurfaceColors=0;

	ribbonColorBase.set( 0, 204, 0 );	// default to green
	ribbonpolymode = GL_FILL;	//GL_LINE;
	ribbonWidthFactor = 1.0;
	stickWidthFactor = 1.0;

	bg.set( 0, 0, 0 );		// background default to black

	cutoffaxe = 'Y';		// surface display defaults set to Y->'Z'
	cutoffval = 0;
	cutoffstep = 2.0;
	surfprobe = 1.4f;		// default radius : water=1.4
	useExistingVtxFile = true;
	showSubSurf = false;
	showCavity = false;
	showSphereNet = true;
	showSurfInMesh = false;

	AAcolorGroups = SURF_COLOR_UNI;	// print in single color by default
	memcpy( AAcolor, AA_DEF_COLOR, sizeof( AAcolor ) );

	exedir[0] = '\0';

	//printf("InputParameters: dispmode=%d\n", dispmode );
}

InputParameters::~InputParameters(  ) {
	DropText( &colorsettings, &colorsettingsLen );
	DropText( &ribboncolor, &ribboncolorLen );
	DropText( &alignsettings, NULL /*&alignsettingsLen */  );	//alignsettingsLen==0 will trig auto_alignment
}

void InputParameters::SwapBite( void *a, int aSize ) {
	unsigned char *s = ( unsigned char * ) a;
	char t;

	aSize--;				// change to index
	for( int i = 0; i <= aSize / 2; i++ ) {
		t = s[i];
		s[i] = s[aSize - i];
		s[aSize - i] = t;
	}
}

// read in parameters from a input instruction file
int InputParameters::ReadInput( char *inputfile, int crdIndex, char *cmoil ) {
	// @TODO: resolve the ambiguity with State::ExeDir and InParm.exedir.  The
	// former is the one to be trusted in ZLAB version (zmoil)

	dispmode = 0;

	char s[LINE_MAXLEN];	//, c[LINE_MAXLEN];
	int k( 0 ), j( 0 ), n( 0 ), q( 0 );
	FILE *fp;
	int currCrd = -1;
	int currWcon = -1;

	if( cmoil != NULL )		// get cmoil's dir path
	{
		char *cptr;
		strcpy( exedir, cmoil );
		cptr = strrchr( exedir, DIRDELIM );
		if( cptr == NULL )
			exedir[0] = '\0';
		else {
			*( cptr + 1 ) = '\0';
		}

		if( ExeDir == NULL )
			ExeDir = exedir;	// assign to global
	}

	fp = fopen( inputfile, "r" );
	if( fp == NULL ) {
		trace( "Unable to open file : %s\n", inputfile );
		exit( -1 );
	}

	int l;

	strcpy( fileName, inputfile );

	while( fgets( s, LINE_MAXLEN, fp ) != NULL ) {
		l = strlen( s ) - 1;
		while( l >= 0 && isspace( s[l] ) ) {
			s[l] = '\0';	// chop '\n'
			l--;
		}
		if( strstr( s, "~" ) ) {
			;				// fgets(c, LINE_MAXLEN, fp);
		}
		else if( strstr( s, "connfile=" ) ) {
			currWcon++;
			if( currWcon == crdIndex ) {
				sscanf( s, "connfile=%s", connName );
				SetValidPath( connName );
			}
		}
		else if( strstr( s, "pthfile=" ) ) {
			currCrd++;
			if( currCrd == crdIndex ) {
				sscanf( s, "pthfile=%s", pthName );
				SetValidPath( pthName );
				filetag = fPTH;
			}
		}
		else if( strstr( s, "crdfile=" ) ) {
			currCrd++;
			if( currCrd == crdIndex ) {
				sscanf( s, "crdfile=%s", crdName );
				SetValidPath( crdName );
				filetag = fCRD;
			}
		}
		else if( strstr( s, "dcdfile=" ) ) {
			currCrd++;
			if( currCrd == crdIndex ) {
				sscanf( s, "dcdfile=%s", dcdName );
				SetValidPath( dcdName );
				filetag = fDCD;
			}
		}
		else if( strstr( s, "xyzfile=" ) )	// coordinate file  in  "x y z"  format
		{
			currCrd++;
			if( currCrd == crdIndex ) {
				sscanf( s, "xyzfile=%s", crdName );
				SetValidPath( crdName );
				filetag = fXYZ;
			}
		}
		else if( strstr( s, "pdbfile=" ) )	// coordinate file  in  "x y z"  format
		{
			currCrd++;
			if( currCrd == crdIndex ) {
				sscanf( s, "pdbfile=%s", pdbName );
				SetValidPath( pdbName );
				filetag = fPDB;
			}
		}
		else if( strstr( s, "bmpfile=" ) ) {
			sscanf( s, "bmpfile=%s", bmpName );
			SetValidPath( bmpName );
		}

		else if( strstr( s, "dispmode" ) || strstr( s, "dispmode" ) )	// dispmode1=* or dispmode2=*
		{

			// @TODO: when moving to zlab, I added the zMsgQueue functions so that
			// the UI will correctly reflect the selected states.  Thus the setting
			// here of displaymode etc should be unnecessary, but subsequnet calls to 
			// ProcessData rely on dispmode already set, so for now I'll leave it;
			// instead we could place calls to ProcessData in the Cmoil_DisplayMode 
			// message handler.  (tfb)


			if( strstr( s + 10, "space" ) ) {
				dispmode |= SPACEBALL;
			}
			else if( strstr( s + 10, "stickball" ) ) {
				dispmode |= STICKBALL;
			}
			else if( strstr( s + 10, "stick" ) ) {
				dispmode |= STICK;
			}
			else if( strstr( s + 10, "ribbon" ) ) {
				dispmode |= RIBBON;
			}
			else if( strstr( s + 10, "backbone" ) ) {
				dispmode |= BACKBONE;
			}
			else if( strstr( s + 10, "structure" ) ) {
				dispmode |= STRUCT2nd;
			}
			else if( strstr( s + 10, "surface" ) ) {
				dispmode |= SURF;
			}
			else {
				dispmode |= STICK;	//default
			}
		}
		else if( strstr( s, "bondbuild=" ) )
			sscanf( s, "bondbuild=%d", &bondmode );
		else if( strstr( s, "swapmode=" ) )
			sscanf( s, "swapmode=%d", &swapmode );
		else if( strstr( s, "skipres=" ) )
			AddText( &skip, s + 8, &skipLen, NULL );
		else if( strstr( s, "skipsurf=" ) )
			AddText( &skipsurf, s + 9, &skipsurfLen, NULL );
		else if( strstr( s, "picksphere=" ) )
			AddText( &picksphere, s + 11, &picksphereLen, NULL );
		else if( strstr( s, "pickradius=" ) )
			sscanf( s, "pickradius=%f", &pickradius );
		else if( strstr( s, "rotatedegree=" ) ) {
			sscanf( s, "rotatedegree=%f", &rotateDegree );
			rotateDegree *= ( float ) AngleToRadius;
		}
		else if( strstr( s, "cavityonly=" ) ) {
			showCavity = atoi( s + 11 ) == 0 ? false : true;
			cutoffval = 1.e+20f;	// no cutoff by default
		}
		else if( strstr( s, "subsurf=" ) )
			AddText( &subsurf, s + 8, &subsurfLen, NULL );
		else if( strstr( s, "colorsetting=" ) )
			AddText( &colorsettings, s + 13, &colorsettingsLen, ", " );	// color separated by ","
		else if( strstr( s, "alignsetting=" ) )
			AddText( &alignsettings, s + 13, &alignsettingsLen, NULL );	// single line pick
		else if( strstr( s, "ribboncolor=" ) )
			AddText( &ribboncolor, s + 12, &ribboncolorLen, ", " );	// color separated by ","
		else if( strstr( s, "stickmode=" ) )
			AddText( &stick, s + 10, &stickLen, NULL );
		else if( strstr( s, "spacemode=" ) )
			AddText( &space, s + 10, &spaceLen, NULL );
		else if( strstr( s, "ribbonmode=" ) )
			AddText( &ribbon, s + 11, &ribbonLen, NULL );
		else if( strstr( s, "nostruct=" ) )
			sscanf( s, "nostruct=%d", &structno );
		else if( strstr( s, "struct2stop=" ) ) {
			sscanf( s, "struct2stop=%d", &structend );	// the order of pth or dvd structure 
			if( structend > 0 )
				sleepseconds = 0;
		}
		else if( strstr( s, "sleepseconds=" ) ) {
			sscanf( s, "sleepseconds=%g", &sleepseconds );
			sleepseconds *= 1000;	// transferm to ms units
		}
		else if( strstr( s, "maxframes=" ) )
			sscanf( s, "maxframes=%d", &maxframes );
		else if( strstr( s, "imagetype=" ) )
			sscanf( s, "imagetype=%3s", imagetype );
		else if( strstr( s, "gifloops=" ) )
			sscanf( s, "gifloops=%d", &gifloops );
		else if( strstr( s, "framedistance=" ) )
			sscanf( s, "framedistance=%d", &framedistance );
		else if( strstr( s, "structpstart=" ) )
			sscanf( s, "structpstart=%d", &structpstart );
		else if( strstr( s, "ribboncolorbase=" ) ) {
			int a, b, c;
			sscanf( s + 17, "%d-%d-%d", &a, &b, &c );
			ribbonColorBase.r = a;
			ribbonColorBase.g = b;
			ribbonColorBase.b = c;
			ribbonColorBase.setGrey(  );
		}
		else if( strstr( s, "ribbonfill" ) )
			ribbonpolymode = GL_FILL;
		else if( strstr( s, "ribbonline" ) )
			ribbonpolymode = GL_LINE;
		else if( strstr( s, "bg=" ) ) {
			int a, b, c;
			sscanf( s + 4, "%d-%d-%d", &a, &b, &c );
			bg.r = a;
			bg.g = b;
			bg.b = c;
			bg.setGrey(  );
		}
		else if( strstr( s, "stereomode=" ) ) {
			StereoMode = atoi( s + 11 ) + NONE_STEREO;
			//const char * names[3] = { "NONE_STEREO", "HARDWARE_STEREO", "ANAGLYPH_STEREO" };
			//zMsgQueue( "type=Cmoil_StereoMode mode=%s", names[ atoi(s+11) ] );
		}
		else if( strstr( s, "cutoffvalue=" ) )
			sscanf( s + 12, "%g", &cutoffval );
		else if( strstr( s, "cutoffaxe=" ) )
			cutoffaxe = s[10];
		else if( strstr( s, "surfprobe=" ) ) {
			sscanf( s + 10, "%g", &surfprobe );
			printf( "%f\n", surfprobe );
		}
		else if( strstr( s, "regeneratevtxfile=" ) )
			useExistingVtxFile = ( atoi( s + 18 ) == 0 );
	}

	if( !IsBinaryCrd( filetag ) )
		maxframes = 1;

	if( structno < 1 && currCrd > 0 )
		structno = 1;

	if( dispmode == 0 ) {
		dispmode |= SKIP;
		// default 
	}

	if( structend > 0 && structno < structend )
		structend = structno;

	if( bmpName[0] == '\0' && crdIndex == 0 )	// set default bmp filename
	{
		if( *connName && filetag != fXYZ && filetag != fPDB )
			strncpy( bmpName, connName, FILENAME_MAXLEN );
		else
			strncpy( bmpName, crdName, FILENAME_MAXLEN );
		bmpName[FILENAME_MAXLEN - 1] = '\0';

		char *pChr = strrchr( bmpName, '.' );
		strcpy( pChr, ".bmp" );
	}

	skipresno = k;
	stickresno = j;
	spaceresno = n;
	ribbonresno = q;
	fclose( fp );

	if( crdIndex > currCrd || ( crdIndex > currWcon && filetag != fXYZ && filetag != fPDB && filetag != fCRD) )	// missing crd or wcon file
		return -1;
	else
		return crdIndex;

}							// ReadInput()

// add string the color var storage, lenVar is allocated memory length fro allocStr 
char *InputParameters::AddText( char **allocStr, char *newStr, unsigned int *lenVar, char *separator ) {
	if( allocStr == NULL || newStr == NULL || lenVar == NULL )
		return NULL;

	if( *allocStr == NULL )	// recognize multi-lines
	{
		*allocStr = ( char * ) malloc( LINE_MAXLEN );
		*lenVar = LINE_MAXLEN;
		( *allocStr )[0] = '\0';
	}
	while( strlen( *allocStr ) + strlen( newStr ) + 1 >= *lenVar )	// null terminated
	{
		*lenVar += LINE_MAXLEN;
		*allocStr = ( char * ) realloc( *allocStr, *lenVar );
	}
	if( *allocStr == NULL ) {
		fprintf( stderr, "CMOIL: Unable to allocate memory!\n" );
		exit( -1 );
	}
	if( ( *allocStr )[0] == '\0' )
		strcpy( *allocStr, newStr );
	else {
		if( separator != NULL )
			strcat( *allocStr, separator );
		strcat( *allocStr, newStr );
	}

	return *allocStr;
}

void InputParameters::DropText( char **allocStr, unsigned int *lenVar ) {
	if( allocStr != NULL && *allocStr != NULL )
		free( *allocStr );
	if( lenVar != NULL )
		*lenVar = 0;
}

}
