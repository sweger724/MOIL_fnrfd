// OPERATING SYSTEM specific includes:
#include "wingl.h"
#include "GL/gl.h"
#include "GL/glu.h"
// SDK includes:
// STDLIB includes:
// MODULE includes:
// ZBSLIB includes:
#include "zvars.h"
#include "zplugin.h"
#include "zmsg.h"
#include "ztime.h"
#include "zgltools.h"
#include "zglmatlight.h"
#include "ztmpstr.h"

//#include "GL/gl.h"
//#include "GL/glu.h"
#include "GL/glfw.h"

#include "stdlib.h"
#include "string.h"

extern void trace( char *fmt, ... );
extern char * getUserLocalFilespec( char *basename, int bMustExist );


ZPLUGIN_BEGIN( gltest );

ZTime timeStartup;
ZTime timeElapsedSinceStartup;
ZTime timeLast = -1.0;
ZTime timeAccum;
GLFWvidmode desktopMode;
char hostname[255];
char username[255];


void sendmail( char *subj, char *body, char *from="thomas@ices.utexas.edu", char *to="blomcode@gmail.com" ) {
	char *tmpfile = getUserLocalFilespec( "mail.txt", 0 );
	FILE *f = fopen( tmpfile, "wt" );
	if( f ) {
		fprintf( f, "to: %s\n", to );
		fprintf( f, "subject: %s\n", subj );
		fprintf( f, "from: %s\n\n", from );
		fprintf( f, body );
		fclose( f );
		#ifndef WIN32
			char buf[255];
			sprintf( buf, "sendmail -t < %s", tmpfile );
			system( buf );
		#endif
	}
}

void render() {
	static int frameCounter = 0;
	static int sentMail = 0;

	static float rotateX = 0;
	static float rotateY = 0;
	static float rotateZ = 0;
	static int w=1,h=1;

	gluPerspective( 30.0, (float)w/(float)h, 4.0, 100.0 );
	gluLookAt( 0.0, 0.0, -5.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0 ); 

	glEnable(GL_COLOR_MATERIAL);
	//glEnable(GL_DEPTH_TEST);
	//glEnable( GL_NORMALIZE );
	glEnable( GL_BLEND );
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

	glEnable(GL_LIGHTING);
//	float ambLight[4] = { 1.f, 1.f, 1.f, 1.f };
	//glLightModelfv( GL_LIGHT_MODEL_AMBIENT, ambLight );

	glEnable(GL_LIGHT0);
	ZGLLight _light0;
	_light0.setLightNumber( 0 );
	_light0.makePositional();
	_light0.pos[0] = .1f;
	_light0.pos[1] = .1f;
	_light0.pos[2] = -.1f;
	_light0.diffuse[0] = 1.f;
	_light0.diffuse[1] = 1.f;
	_light0.diffuse[2] = 1.f;
	_light0.active = 1;
	glEnable( GL_LIGHT0 );
	_light0.setGL();

	rotateX += 1.f;
	rotateY += 2.f;
	rotateZ += 3.f;
	glRotatef( rotateX, 1.f, 0.f, 0.f );
	glRotatef( rotateY, 0.f, 1.f, 0.f );
	glRotatef( rotateZ, 0.f, 0.f, 1.f );
	glColor4f( .8f, .1f, .1f, .7f );
	//glScalef( 3.0, 3.0, 3.0 );
	zglCube( 1.f );
	// A cube centered on the origin with length of side equal to dim

	ZTime timeNow = zTimeNow();
	if( timeLast < 0 ) {
		timeLast = timeNow;
	}
	timeElapsedSinceStartup = timeNow - timeStartup;
	ZTime t = timeAccum;
	if( 0 && !sentMail && timeElapsedSinceStartup > 20 ) {
		zMsgQueue( "type=QuitApp" );
	}
	else if( timeAccum > 3.0 && !sentMail ) {
		// display fps
		double fps = frameCounter / timeAccum;
		trace( "[ %04d x %04d ] FPS: %.3lf\n", w, h, fps );
		frameCounter = 0;

		if( fps < 150.0 ) {
			sendmail( "gltest", ZTmpStr( "%s, %s: gltest avg render over last 3 seconds was %.4lf\n", hostname, username, fps ) );
			trace( " *Sent an email about low framerate condition.\n" );
			sentMail = 1;
		}
		else {
			// switch window pos/size
			int x = rand() % 50;
			int y = rand() % 50;
			glfwSetWindowPos( x, y );
			w = rand() % ( desktopMode.Width - x );
			h = rand() % ( desktopMode.Height - y );
			glfwSetWindowSize( max(w,10), max(h,10) );
			timeAccum = 0.0;
		}
	}
	timeAccum += timeNow - timeLast;
	timeLast = timeNow;

	frameCounter++;
}

void startup() {
	timeStartup = zTimeNow();
	srand( (unsigned int )timeStartup );
	glfwGetDesktopMode( &desktopMode );

	char *hostVar, *userVar;
#ifdef WIN32
	hostVar = "COMPUTERNAME";
	userVar = "USERNAME";
#else
	hostVar = "HOSTNAME";
	userVar = "USER";
#endif
	char *h = getenv( hostVar );
	char *u = getenv( userVar );
	if( !h ) {
		h = "unknown host";
	}
	if( !u ) {
		u = "unknown user";
	}
	strcpy( hostname, h );
	strcpy( username, u );
}

void shutdown() {
}

void handleMsg( ZMsg *msg ) {
}

ZPLUGIN_EXPORT_PROPERTY( shadowGardenInterface, "1" );
ZPLUGIN_EXPORT_SYMBOL( startup );
ZPLUGIN_EXPORT_SYMBOL( shutdown );
ZPLUGIN_EXPORT_SYMBOL( render );
ZPLUGIN_EXPORT_SYMBOL( handleMsg );

ZPLUGIN_END;

