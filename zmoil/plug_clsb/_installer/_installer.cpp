// @ZBS {
//		*WIN32_LIBS_DEBUG opengl32.lib glu32.lib
//		*WIN32_LIBS_RELEASE opengl32.lib glu32.lib
//		*SDK_DEPENDS freeimage pcre
//		*MODULE_DEPENDS zuitexteditor.cpp zuifilepicker.cpp zwinutil.h  
// }

// OPERATING SYSTEM specific includes:
#include "wingl.h"
#include "GL/gl.h"
#include "GL/glu.h"

#ifdef WIN32
#undef APIENTRY
#undef WINGDIAPI
#endif


// STDLIB includes:
#include "stdio.h"
#include "stdlib.h"

// MODULE includes:

// ZBSLIB includes:
#include "zhashtable.h"
#include "zplugin.h"
#include "zui.h"
#include "ztmpstr.h"
#include "zwinutil.h"

#include "string.h"		// strlen, debug	

//extern void trace( char *, ... );
extern char * getUserLocalFilespec( char *basename, int bMustExist );
extern ZHashTable options;
	// all imported from main.cpp

// The following values will be replaced post-compile with ascii numeric values
// that will be read by this progam to locate portions of the file to extract
// to the payload.zip and unzip.exe files which are used during installation.
char * payloadBegin = "PAYLOADB";
char * payloadLen   = "PAYLOADL";
#ifdef WIN32
char * unzipBegin   = "UNZIPBEG";
char * unzipLen     = "UNZIPLEN";
#endif


//==================================================================================================
// Forward declares (outside plugin namespace); utility functions; see end of this file for impl.

char * getVersionString();
	// get the version string that is embedded in the TITLE macro

//==================================================================================================
// Plugin

ZPLUGIN_BEGIN( installer );

void setMOILHOME( char *moilHome ) {
	int hresult = zWinUtilSetRegString( "HKEY_LOCAL_MACHINE\\System\\CurrentControlSet\\Control\\Session Manager\\Environment\\MOIL_HOME\0", moilHome );
}

void setMOILPath( char *guiPath, char *exePath ) {
	// add the given paths to the environment, replacing any existing moil-related paths.
	// One or both arguments may be NULL, and thus can be used to remove all moil-related paths.

	// GET PATH ENV & split into individual paths
	char buffer[8192];
	int hresult = zWinUtilGetRegString( "HKEY_LOCAL_MACHINE\\System\\CurrentControlSet\\Control\\Session Manager\\Environment\\Path\0", buffer, 8191 );
	buffer[8191]=0;
	ZStr *paths = zStrSplitByChar( ';', buffer );

	// WALK LIST, REPLACE first occurance of \moil.gui with new guiPath, remove remaining,
	// adding to end if no \moil occurance is found.
	buffer[0]=0;
	int found = 0;
	int count = zStrCount( paths );
	for( int i=0; i<count; i++ ) {
		char *p = paths->getS( i );
		if( strstr( p, "\\moil\\" ) || strstr( p, "\\moil.gui" ) || strstr( p, "\\moil.exe" ) || strstr( p, "\\moil.source\\exe" ) ) {
			if( !found ) {
				if( guiPath ) {
					strcat( buffer, guiPath );
					strcat( buffer, ";" );
				}
				if( exePath ) {
					strcat( buffer, exePath );
					strcat( buffer, ";" );
				}
				found = 1;
			}
		}
		else {
			strcat( buffer, p );
			strcat( buffer, ";" );
		}
	}
	if( !found ) {
		if( guiPath ) {
			strcat( buffer, guiPath );
			strcat( buffer, ";" );
		}
		if( exePath ) {
			strcat( buffer, exePath );
		}
	}

	int len = strlen( buffer );
	while( len && buffer[--len] == ';' ) {
		buffer[len] = 0;
	}
	
	// SET new path
	hresult = zWinUtilSetRegString( "HKEY_LOCAL_MACHINE\\System\\CurrentControlSet\\Control\\Session Manager\\Environment\\Path\0", buffer );
}

int writePayload( int offset, int len, char *filespec ) {
	if( offset && len ) {
		FILE *fw = fopen( filespec, "w" );
		FILE *fr = fopen( options.getS( "argv0" ), "r" );
		if( fw && fr ) {
			char buffer[4096];
			int bytesRead;
			while( ( bytesRead=fread( buffer, 1, 4096, fr ) ) != 0 ) {
				fwrite( buffer, 1, bytesRead, fw );
			}
			fclose( fr );
			fclose( fw );
			return 1;
		}
	}
	return 0;
}

void startup() {
	printf( "startup...\n");

	// LOAD textures
	// chevronTex = zglGenTextureRGBAZBits( "_kin/chevron2.png", 0, 1 );

	// LOAD fonts
	// zglFontLoad( "verdana12", "verdana.ttf", 12 );

	// LOAD the ui definition
	ZUI::zuiExecuteFile( ZTmpStr( "_%s/_%s.zui", thisPluginName, thisPluginName ) );

	char buffer[8192];
	int hresult = zWinUtilGetRegString( "HKEY_LOCAL_MACHINE\\System\\CurrentControlSet\\Control\\Session Manager\\Environment\\Path\0", buffer, 8191 );
	buffer[8191]=0;
	printf( "\n\nPath is:\n%s\n\n", buffer );

	printf( "Will now set moil-related path to d:\\thomas\\moil\\moil.gui...\n\n" );
	setMOILPath( "d:\\thomas\\moil\\moil.gui", "d:\\thomas\\moil\\moil.exe" );

	hresult = zWinUtilGetRegString( "HKEY_LOCAL_MACHINE\\System\\CurrentControlSet\\Control\\Session Manager\\Environment\\Path\0", buffer, 8191 );
	buffer[8191]=0;
	printf( "Path is now:\n%s\n\n", buffer );

	printf( "Now clearing the moil path...\n\n" );
	setMOILPath( 0, 0 );

	hresult = zWinUtilGetRegString( "HKEY_LOCAL_MACHINE\\System\\CurrentControlSet\\Control\\Session Manager\\Environment\\Path\0", buffer, 8191 );
	buffer[8191]=0;
	printf( "Path is now:\n%s\n\n", buffer );

	hresult =  zWinUtilGetCommonProgramsFolderPath( buffer );
	printf( "the Progarms folder path is:\n%s\n\n", buffer );
}

void shutdown() {
}

void update() {
}

void render() {
}


// TEMPLATE (generic) Message Handlers
//----------------------------------------------------------------------------------

ZMSG_HANDLER( ZLAB_Shutdown ) {
	shutdown();
	exit( 0 );
}

ZMSG_HANDLER( TutorialOK ) {
}

ZPLUGIN_EXPORT_PROPERTY( shadowGardenInterface, "1" );
ZPLUGIN_EXPORT_SYMBOL( startup );
ZPLUGIN_EXPORT_SYMBOL( shutdown );
ZPLUGIN_EXPORT_SYMBOL( render );
ZPLUGIN_EXPORT_SYMBOL( update );


ZPLUGIN_END;

char * getVersionString() {
	// get the version string that is embedded in the TITLE macro
#ifdef TITLE
	static ZRegExp versionStr( "\\S+ Version [0-9]\\.[0-9]\\.[0-9]+" );
	versionStr.test( TITLE );
	return versionStr.get( 0 );
#else
	static char *unk = "Unknown";
	return unk;
#endif
}



