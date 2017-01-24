// @ZBS {
//		*INTERFACE gui
//		*WIN32_LIBS_DEBUG opengl32.lib glu32.lib
//		*WIN32_LIBS_RELEASE opengl32.lib glu32.lib
//		*SDK_DEPENDS 
//		*MODULE_DEPENDS _zmoil.h zhash3d.cpp tbutil.cpp modelList.cpp zmat.cpp zuiradiobutton.cpp zuicheck.cpp zuifilepicker.cpp zuivar.cpp zuitexteditor.cpp zuislider.cpp zuitabbutton.cpp zuitabpanel.cpp zglgfiletools.cpp cmoil.cpp cmoil_aminoacids.cpp cmoil_atom.cpp cmoil_bspline.cpp cmoil_drawshape.cpp cmoil_msgboard.cpp cmoil_parm.cpp cmoil_pick.cpp cmoil_util.cpp trackball.cpp
//		*EXTRA_DEFINES SFFAST NATIVE_FILE_DIALOGS
// }

// TFB: some original cmoil modules have been removed because we prefer to implement their functionality using
// zlab, zbslib, or functionality in FreeImage etc...  These include  cmoil_dib, cmoil_menu, cmoil_movieprint,
// gifmerge, and possibly others. 

// OPERATING SYSTEM specific includes:
#include "wingl.h"
#include "GL/gl.h"
#include "GL/glu.h"
// SDK includes:
#include "FreeImage.h"
// STDLIB includes:
#include <cmath>
#include "sys/stat.h"
#ifdef WIN32
	#include <malloc.h>
	#include "direct.h"
#endif

// CMOIL includes
#include "_zmoil/cmoil.h"
#include "_zmoil/cmoil_atom.h"
#include "_zmoil/cmoil_pick.h"

// MODULE includes:
#include "_zmoil.h"


// ZBSLIB includes:
#include "zvars.h"
#include "zplugin.h"
#include "zmsg.h"
#include "zui.h"
#include "zwildcard.h"
#include "zfilespec.h"
#include "zmousemsg.h"
#include "zconfig.h"
#include "zglfont.h"
#include "ztmpstr.h"
#include "zbitmapdesc.h"
#include "zgraphfile.h"
#include "zgltools.h"
#include "ztime.h"
#include "zmat.h"
#include "zregexp.h"
#include "mainutil.h"
#include "tbutil.h"

extern int zMouseShift;
extern int zuiKeyEvent;


#include "modelList.h"

#ifdef __USE_GNU
#define stricmp strcasecmp
#endif

#define _ZVAR(type,name,val) type name = (type)val
#define _ZVARB(type,name,val,minb,maxb) type name = (type)val
	// I define these so that I can easily remove a variable from the editable list by
	// prepending an underscore to the macro name.  

ZVARB( float, Zmoil_BondLenMax, 4.f, 1.f, 100.f );
ZVARB( int, Zmoil_PeptidePlaneAtoms, 0, 0, 1 );
_ZVARB( int, Zmoil_PeptidePlaneFlipping, 1, 0, 1 );
ZVARB( int, Zmoil_PeptidePlaneVectors, 0, 0, 1 );
ZVARB( int, Zmoil_PeptidePlaneWireframe, 0, 0, 1 );
ZVARB( int, Zmoil_Perspective, 1, 0, 1 );
ZVARB( int, Zmoil_PthDisplayStructs, 1, 1, 64 );
ZVARB( int, Zmoil_PthDisplayAlphaLinear, 1, 0, 1 );
_ZVARB( int, Zmoil_DisplayWater, 0, 0, 1);
_ZVARB( int, Zmoil_DisplayHBonds, 0, 0, 1);
ZVARB( int, Zmoil_DisplayHBondsWater, 0, 0, 1);
ZVARB( int, Zmoil_HBondThickness, 2, 1, 5);
_ZVARB( float, Zmoil_HBondDistance, 2.f, 1.f, 4.f );
_ZVARB( float, Zmoil_SpaceFillingScale, 1.f, .1f, 10.f );
_ZVARB( int, Zmoil_VelocityDisplay, 0, 0, 1 );
ZVARB( int, Zmoil_VelocityCAOnly, 0, 0, 1 );
ZVARB( float, Zmoil_VelocityScale, 4.f, 0.1f, 10.f );
ZVARB( float, Zmoil_VelocityMin, 0.01f, 0.01f, 100.f );



ZVARB( int, Zmoil_UIdirtyRects, 0, 0, 1 );
ZVARB( int, Zmoil_UIdirtyRectsDebug, 0, 0, 1 );
ZVARB( int, Zmoil_DrawSheets, 1, 0, 1 );
ZVARB( float, Zmoil_eyeZ, 50.f, 1.f, 500.f );
ZVARB( float, Zmoil_CavityMinSurfArea, 40.f, 0.f, 500.f );

ZVARB( int, Zmoil_MFMView, 0, 0, 3 );
_ZVARB( int, Zmoil_gridColor, 1, 0, 1 );
ZVAR( double, Zmoil_gridColorMax, 700 );
ZVAR( double, Zmoil_gridColorMin, 200 );


struct CmoilModel {
	ZHashTable viewInfo;
		// Temporary information used by the view and zuis
	ZHashTable modeNameToValue;
		// lookup enum value of CMOIL::DisplayMode by name
	ZHashTable aminoNameToChar;
		// get one-letter character from 3 letter amino abbrev
} model;


#define UI_PANEL_WIDTH (360)

ZUI *g_pluginPanel = NULL;
ZUI *g_uiPanel = NULL;
int gEnableAxes=0;
int gEnableStructInfoDisplay=1;
//int gModelAlignMasterIndex=-1;
ZTime gLastActivityTime;
//float g_EyeZ;
int g_ShiftKey=0;
int g_AltKey=0;
double gTimeout = -1.0;

CMOIL::AtomArray *manualAlign=0;

char EXEFolder[255] = {0};
extern char exeName[255];

extern int useDirtyRects;
extern int drawDirtyRects;
	// from zui.cpp : controls use of dirty rect system; we control
	// these with Zmoil_UIdirtyRects, Zmoil_UIdirtyRectsDebug

//-----------------------------------------------------------------------------
// util functions, forward declare.  see end of file for implementation

void glPrint( int x, int y, char *text, char *font=0 );
	// handy print in screen coords
	// @TODO: replace this with simple zglFontPrint or some such

char tmpMessage[256] = {0,};
ZTime tmpMessageDisplayUntil;
unsigned char tmpMessageColor[3];
	// for displaying temporary messages onscreen

unsigned char * capturePluginImage( char *filename, bool notify=true );
	// take a screen shot of the plugin window

void updateImageCapture();
	// called by render for multiframe captures

int operationPending( int inc=0 );
	// a mechanism to track how many frames until critical opertaions
	// have completed; this is used for instance when loading an mfm, 
	// which takes ~5 frames due to various messaging etc, and some 
	// other operation, such as an automated screen capture, will want
	// to wait to capture until all opertaions are complete.  Pass the
	// number of frames needed by a current operation, or call with no 
	// args to get frames left until complete.

char * getVersionString() {
	// get the version string that is embedded in the TITLE macro
#ifdef TITLE
	static ZRegExp versionStr( "\\S+ Version [0-9]+\\.[0-9]+\\.[0-9]+" );
	versionStr.test( TITLE );
	return versionStr.get( 0 );
#else
	static char *unk = "(Unknown Version)";
	return unk;
#endif
}

void updateIndexSlider( bool labelOnly=true );

//-----------------------------------------------------------------------------
// We port CMOIL to zlab by implementing a ZUI that renders via the CMOIL
// openGL rendering code.  We use a ZUI instead of the plugin render() callback
// so that we can do things like requestExclusiveMouse, and potentially layout
// several such molecular views in future.

struct ZUIRenderCmoil : ZUIPanel {
	ZUI_DEFINITION(ZUIRenderCmoil,ZUIPanel);
	ZUIRenderCmoil() : ZUIPanel() { }
	virtual void factoryInit();
	virtual void update();
	virtual void render();
	virtual void handleMsg( ZMsg *msg );
};
ZUI_IMPLEMENT( ZUIRenderCmoil, ZUIPanel );

void ZUIRenderCmoil::factoryInit() {
	ZUIPanel::factoryInit();
	setUpdate( 1 );
}

void ZUIRenderCmoil::update() {
	dirty();
		// for now, always dirty the render window; we could be smarter
		// and only render when a display state is changed or the camera
		// is being adjusted.
}

int gEnableFog=1;
float cmoilPickX = -1.f;
float cmoilPickY = -1.f;
ZTime globalTime;
int globalTimerOn=0;
void ZUIRenderCmoil::render() {

	static ZTime startup = zTimeNow();
	globalTime = zTimeNow();
		
	if( !g_pluginPanel || w<=0 || h<=0 ) {
		return;
	}

	// Don't render molecules when fileopen dialog is up, so that the dialog is as snappy
	// as possible 
	//
	ZUI *filePicker = ZUI::zuiFindByName( "zuiFilePickerDialog" );
	if( (filePicker && !filePicker->getI( "hidden" )) || CMOIL::Main::NumImages < 1 ) {
		glClear( GL_COLOR_BUFFER_BIT );
		return;
	}

	// Some debugging/profile stuff...
	//
	INIT_TIMERS( 5 );
	if( ( nTimed < 5 || globalTimerOn ) ) {
		trace( "-------------------------------------- ZUIRenderCmoil::render()\n" );
	}
	START_TIMER( ZUIRender );
	START_TIMER( ZUIBegin );

	CMOIL::AtomArray *aaptr = CMOIL::Main::AAptrs[0];
	
	begin3d();
		// save opengl matrices/attribs; do viewport transform

	// TODO: this block might be more appropriately placed in CMOIL::Main::display
	////////////////////////////////////////////////////////////////////////////////
	
	float objectBoundRadius = 30; 
	if( aaptr && Zmoil_eyeZ < 1.f ) {
		objectBoundRadius = (float)aaptr->maxDisplacementXYZ * 1.33f;		
		Zmoil_eyeZ = objectBoundRadius / tanf( 45.f / 2.f );
	}
	END_TIMER( ZUIBegin );

	
	// SETUP perspective projection for molecule render; modelview matrix
	// is handled in call to Main::display() below
	START_TIMER( ZUITrans );
	glMatrixMode(GL_PROJECTION);
	if( Zmoil_Perspective ) {
		gluPerspective( 45.0, (float)w/(float)h , 1.0, 10000.0 );
	}
	else {
		glOrtho( -objectBoundRadius/2, objectBoundRadius/2, -objectBoundRadius/2, objectBoundRadius/2, 1, 10000 );
	}
	END_TIMER( ZUITrans );

	////////////////////////////////////////////////////////////////////////////////
	START_TIMER( ZUICallDisplay );
	if( ( nTimed < 5 || globalTimerOn ) )
	trace( " *** globalTime in zmoil: %9.8lf\n", globalTime );

	if( gEnableFog ) {
		GLfloat fogColor[4] = { 0, 0, 0, 1 };
		glFogi( GL_FOG_MODE, GL_LINEAR );
		glFogfv( GL_FOG_COLOR, fogColor );
		glHint( GL_FOG_HINT, GL_FASTEST );
		//glFogf( GL_FOG_START, objectBoundRadius*1.5 );
		//glFogf( GL_FOG_END, max( 1.f, objectBoundRadius*2.5 ) );
		float fogStart = Zmoil_eyeZ; //+  aaptr ? aaptr->centroid.z/2.0 : 0 );
		float fogEnd   = fogStart + objectBoundRadius;
		glFogf( GL_FOG_START, fogStart );
		glFogf( GL_FOG_END, fogEnd );
		glEnable(GL_FOG);
	}

	CMOIL::Main::display( (int)w, (int)h );
	if( gEnableFog ) {
		glDisable( GL_FOG );
	}
	static int slowRender = 0;
	if( ( nTimed < 5 || globalTimerOn ) ) {
		ZTime t = zTimeNow();
		double fps = 1.0 / ( t - globalTime );
		trace( " *** %9.8lf back in zmoil = %6.3lf fps (render only!)\n", t - globalTime, fps );
		if( nTimed > 3 && fps < 20.0 ) {
			//emailMsgTo( "Low framerate in zmoil!",  "blomcode@gmail.com" );
			int debug=1;
			#ifdef WIN32
			//_asm {
			//	int 3;
			//}
			#else
			//asm ( "int3" );
			// how to do the above without it causing the app to close immediately?
			// do we have to be attached with the debugger already?
			
			slowRender++;
			trace( "slowRender: %.4lf frames/sec\n", fps );
			#endif
		}
	}
	if( !slowRender && gTimeout>0  && ( globalTime - startup > gTimeout ) ) {
		zMsgQueue( "type=Zmoil_Shutdown" );
	}
	globalTimerOn=0;
	END_TIMER( ZUICallDisplay );

	// Requests to pick atom at location are made by setting coords to non-negative.
	// The pick is performed in the context of this render so that projection and
	// modelview matrices are intact for pickray projection.  It would be nicer
	// to pass matrices to the PickQ object and handle it in ZUIMouseClickOn.

	if( cmoilPickX >= 0.f ) {
		CMOIL::Main::PickQ.put( (int)cmoilPickX, (int)cmoilPickY );
		cmoilPickX = cmoilPickY = -1.f;
	}

	end3d();
		// restore opengl matrices/attribs

	if( Zmoil_VelocityDisplay && aaptr && aaptr->properties.has( "velocityMax" ) ) {
		spectralLegend( x+100, h-250, 30.f, 200.f, 0.f, (float)aaptr->properties.getD( "velocityMax"  ));
	}
	else if ( Zmoil_gridColor && aaptr ) {
		spectralLegend( x+10, h-250, 30.f, 200.f, /*aaptr->boxInfoMin*/ (float)Zmoil_gridColorMin, /*aaptr->boxInfoMax*/ (float)Zmoil_gridColorMax );
		static double gridColorMaxLast = 0;
		static double gridColorMinLast = 0;
		if( Zmoil_gridColorMax != gridColorMaxLast || Zmoil_gridColorMin != gridColorMinLast ) {
			Zmoil_gridColorMax = max( Zmoil_gridColorMax, .0001 );
			Zmoil_gridColorMin = max( Zmoil_gridColorMin, .0001 );
			gridColorMaxLast = Zmoil_gridColorMax ;
			gridColorMinLast = Zmoil_gridColorMin ;
			zMsgQueue( "type=Zmoil_BoxColor colors=1" );
		}
	}


	// tmalign resuylts printout:
	//
	float wx, wy;
	getWindowXY( wx, wy );
	if( ZUI::zuiFindByName( "AlignButton" )->getI( "selected" ) ) {
		int py = 30;
		glColor3ub( 200, 200, 200 );
		glPrint( (int)wx + 5, (int)h - py, "tmscore" );
		glPrint( (int)wx + 70, (int)h - py, "model" );
		py += 15;

		double avgScore=0.0,looppScore=0.0;
		int avgCount=0;
		int looppY;
		int alignedCount=0;
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ i ];
			if( aa->bDraw ) {
				ZFileSpec fs( aa->InParm.fileName );
				CMOIL::tRGBA &color = CMOIL::zmoilColors[ aa->properties.getI( "colorIndex", CMOIL::ZC_White ) ];
				glColor3ub( color.r, color.g, color.b );
				char *score = "master";
				ZTmpStr rank( "" );
				int lp8rank;
				if( lp8rank=aa->properties.getI( "lp8rank" ) ) {
					rank.set( " LP8=%d LP9=%d", lp8rank, aa->properties.getI( "lp9rank" ) ); 
				}
				if( aa->properties.getI( "tmalignMaster" ) ) {
					alignedCount++;
				}
				else if( aa->alignToModel ) {
					alignedCount++;
					double s = aa->properties.getD( "tmscore", 0.0 );
					score = formatFloat( s );

					char *type = aa->properties.getS( "type", "" );
					if( !strncmp( type, "model", 5 ) ) {
						// this is a model (not a template, or a native)
						if( strncmp( type, "model_hi", 8 ) )  {
							// this is type 'model', not one of the special models (e.g. loopp, or best model)
							avgScore += s;
							avgCount++;
						}
						else if( !strcmp( type, "model_hi1" ) ) {
							looppScore = s;
							looppY = py;
								// this is the loopp model (or best loopp model)
						}
					}
				}
				if( alignedCount ) {
					glPrint( (int)wx + 5, (int)h - py, score );
	
					glPrint( (int)wx + 70, (int)h - py, ZTmpStr( "%s%s%s", fs.getFile(), aa->properties.getI( "looppBestModel" ) ? "*" : "", rank.s ) );
					py += 15;
				}
			}
		}
		if( alignedCount ) {
			py = 30;
			glColor3ub( 200, 200, 200 );
			glPrint( (int)wx + 5, (int)h - py, "tmscore" );
			glPrint( (int)wx + 70, (int)h - py, "model" );
		}
//		if( avgCount > 0 ) {
//			glPrint( (int)x + 200 + 50, (int)h - looppY, ZTmpStr( "%s avg", formatFloat( looppScore - avgScore/avgCount ) ) );
//		}
	}

	// Hack temp messages here:
	glColor3ub ( tmpMessageColor[0], tmpMessageColor[1], tmpMessageColor[2] );
	if( *tmpMessage ) {
		if( tmpMessageDisplayUntil > zTimeNow() ) {
			glPrint( (int)wx + 5, (int)h - 15, tmpMessage );
		}
		else {
			tmpMessage[0] = 0;
		}
	}

	if( model.viewInfo.getI( "captureFrames" ) ) {
		updateImageCapture();
	}
	CMOIL::Main::animate();
		// in the case of "movie" mode, spinning, etc...

	TRACE_TIMERS();
}

void ZUIRenderCmoil::handleMsg( ZMsg *msg ) {
	gLastActivityTime = zTimeNow();

	// @TODO: might be nicest to use a ZViewpoint on this instead of
	// the custom quaternion stuff in cmoil.  For now just porting to zlab.

	if( zmsgIs(type,ZUIMouseClickOn) ) {
		requestExclusiveMouse( 1, 1 );
		if( zmsgIs(which,L) &&  zmsgIs(dir,D) ) {
			if( zmsgI(shift) ) {
				// pick
				cmoilPickX = zmsgF( localX );
				cmoilPickY = h - zmsgF( localY );
					// pick operation happens in context of render for transforms
			}
			else {
				// rotate
				CMOIL::Main::Spinning = 0;
				CMOIL::Main::Moving = 1;
				CMOIL::Main::BeginX = (int)zmsgF( localX );
				CMOIL::Main::BeginY = (int)(h - zmsgF( localY ));
			}
		}
		else if( zmsgIs(which,R) && zmsgIs(dir,D) ) {
			// translate
			CMOIL::Main::Translating = 1;
			CMOIL::Main::BeginX = (int)zmsgF( localX );
			CMOIL::Main::BeginY = (int)(h - zmsgF( localY ));
		}
		else if( zmsgIs(which,M) && zmsgIs(dir,D) ) {
			// no action for middle click: wheel zooms, though, see below
		}
		zMsgUsed();
	}
	else if( zmsgIs(type,ZUIExclusiveMouseDrag) ) {
		CMOIL::Main::motion( (int)zmsgF( localX ), (int)(h - zmsgF( localY )), w, h );	
		zMsgUsed();

	}
	else if( zmsgIs(type,ZUIMouseReleaseDrag) ) {
		if( zmsgIs(which,R) ) {
			CMOIL::Main::Translating = 0;
		}
		else if ( zmsgIs(which,L) ) {
			CMOIL::Main::Moving = 0;
		}
		zMsgUsed();
	}
	else if( (zmsgIs(type,Key) && zmsgIs(which,wheelforward)) ) {
		float localX = (float)zMouseMsgX;
		float localY = (float)zMouseMsgY;
		getLocalXY( localX, localY );

		if( containsLocal( localX, localY ) ) {
			Zmoil_eyeZ *= 1.15f;
			zMsgUsed();
		}
	}
	else if( (zmsgIs(type,Key) && zmsgIs(which,wheelbackward)) ) {
		float localX = (float)zMouseMsgX;
		float localY = (float)zMouseMsgY;
		getLocalXY( localX, localY );

		if( containsLocal( localX, localY ) ) {
			Zmoil_eyeZ *= .85f;
			zMsgUsed();
		}
	}
	else if( (zmsgIs(type,Key) && zmsgIs(which,ctrl_t)) ) {
		globalTimerOn=1;
		zMsgUsed();
	}


}

struct ZUIResidueText : ZUIPanel {
	#define MAX_MODELS 8
	
	int selectedTextIndex[ MAX_MODELS ];
	int selectedResidue[ MAX_MODELS ];
	float charWidth;
	float charHeight;


	ZUI_DEFINITION(ZUIResidueText,ZUIPanel);
	ZUIResidueText() : ZUIPanel() {}
	virtual void factoryInit() {
		ZUIPanel::factoryInit();
		charWidth = charHeight = -1.f;
		for( int i=0; i<MAX_MODELS; i++ ) {
			selectedTextIndex[ i ] = -1;
			selectedResidue[ i ] = -1;
		}
		setUpdate( 1 );
	}
	virtual void update();
	virtual void render();
	virtual void handleMsg( ZMsg *msg );
	void selectResidueFromMouseXY( int mX, int mY );
		// given mouse coordinates, determine the structure & residue index
};
ZUI_IMPLEMENT( ZUIResidueText, ZUIPanel );

void ZUIResidueText::update() {
	// See if we need to update the parent panel height.

	/*
	int count=0;
	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		if( CMOIL::Main::AAptrs[ i ]->bDraw ) {
			count++;
		}
	}
	float height = 13.f * count;
	if( h != height ) {
		parent->putF( "layout_forceH", height );
		parent->h = height;
		putF( "layout_forceH", height );
		h = height;
		layoutRequest();
	}
	*/
	dirty();
		// todo: optimize
}

void ZUIResidueText::render() {
	if( !g_pluginPanel || w<=0 || h<=0  ) {
		return;
	}
	float d;
	if( charWidth < 0 ) {
		zglFontGetDimChar( 'Q', charWidth, charHeight, d, "courier10" );
			// we use a fixed-pitch font, so which letter doesn't really matter...
	}

	ZUIPanel::render();
	begin3d();
	// PRINT RESIDE (amino acide) sequence
	static char residueText[4096];
	int scrollXHome = getI( "scrollXHome" );
	int i;
	int dispCount=0;
	for( i=0; i<CMOIL::Main::NumImages && dispCount < 2; i++ ) {
		residueText[0]=0;
		CMOIL::AtomArray *aaptr = CMOIL::Main::AAptrs[ i ];
		if( !aaptr->bDraw ) {
			continue;
		}
		if( !aaptr->alignSequence[0] ) {

			// GOOD GRIEF: there has to be a better way to scan the residues than atom by atom!
			char *code;
			int lastRes = -1;
			int resCount = 0;
			for( int j=0; resCount<4095 && j<aaptr->numAtom; j++ ) {
				if( aaptr->atomList[j].skip ) {
					continue;
				}
				int resN = aaptr->resNum( aaptr->atomList + j );
				if( resN != lastRes ) {
					code = model.aminoNameToChar.getS( aaptr->resName[ resN-1 ], 0 );
					if( !code ) {
						code = "*";
					}
					residueText[resCount++]=code[0];
					lastRes = resN;
				}
			}
			residueText[resCount]=0;
		}
		else {
			strcpy( residueText, aaptr->alignSequence );
		}
		glColor3ub( 200,200,200 );
		glPrint( (int)(scrollX - scrollXHome), (int)(h - charHeight*(dispCount+1)), residueText, "courier10" );
		putI( ZTmpStr( "modelForSequence%d", dispCount ), i );
		dispCount++;
	}
		
	// Hilight the selected residue
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
	zglPixelMatrixInvertedFirstQuadrant();
	for( i=0; i<CMOIL::Main::NumImages; i++ ) {
		if( selectedTextIndex[i] >= 0 ) {
			i==0 ? glColor4ub( 255, 127, 0, 96 ) : glColor4ub( 255, 0, 255, 96 );
			float x0 = selectedTextIndex[i] * charWidth + scrollX - scrollXHome;
			float y0 = i * charHeight;
			glBegin( GL_QUADS );
				glVertex2f( x0, y0 );
				glVertex2f( x0 + charWidth, y0 );
				glVertex2f( x0 + charWidth, y0 + charHeight );
				glVertex2f( x0,y0 + charHeight );
			glEnd();

		}
	}

	end3d();
}

void ZUIResidueText::handleMsg( ZMsg *msg ) {
	if( zmsgIs(type,ZUIMouseClickOn) ) {
		if( zmsgIs(which,L) && zmsgIs(dir,D) ) {
			selectResidueFromMouseXY( (int)zmsgF( localX ), (int)zmsgF( localY ) );
			zMsgUsed();
		}
	}

	if( !zMsgIsUsed() ) {
		ZUIPanel::handleMsg( msg );
	}
}

void ZUIResidueText::selectResidueFromMouseXY( int mX, int mY ) {
	assert( charWidth > 0 );
	int scrollXHome = getI( "scrollXHome" );
	mX += (int)(scrollXHome - scrollX);
	int mi = (int)(( h - mY ) / charHeight);  // model index (which line)
	int ti = (int)(mX / charWidth);		   // text index  (which char)

	selectedTextIndex[ mi ] = ti;
	selectedResidue[ mi ] = ti;

	int modelIndex = getI( ZTmpStr( "modelForSequence%d", mi) );

	CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ modelIndex ];
	if( aa && aa->alignSequence[0] ) {
		if( aa->alignSequence[ti] == '-' ) {
			selectedResidue[ mi ] = -1;
		}
		else {
			for( int i=0; i<ti; i++ ) {
				if( aa->alignSequence[i] == '-' ) {
					selectedResidue[mi]--;
				}
			}
		}
	}
	if( aa ) {
		aa->setPickedResidue( selectedResidue[mi] );
	}
}

void syncUIToInputParams() {
	// READ the CMOIL::Main::InParm struct and set the UI to match the settings
	// that were read from an input file.
	
	CMOIL::InputParameters *parms = CMOIL::Main::AAptrs[0] ? &CMOIL::Main::AAptrs[0]->InParm : 0;
	
	// Setup default basename for image capture
	char *capName = model.viewInfo.getS( "filename", 0 );
	if( !capName ) {
		capName = "zmoil";
	}
	ZFileSpec captureSpec( capName );
	model.viewInfo.putS( "captureFile", replaceExtension( capName, "png" ) );
	// default to png file
	ZUI::zuiFindByName( "captureFilename" )->putS( "text", model.viewInfo.getS( "captureFile" ) );
	
	// White Background
	if( parms && parms->bg.r == 0xFF ) {
		ZUI::zuiFindByName( "gen1" )->putI( "selected", 1 );		
	}
	
	// Shininess
	ZUI::zuiFindByName( "gen2" )->putF( "val", CMOIL::Main::Shine );
	
	// StereoMode
	char *buttonName = "st1";
	if( CMOIL::Main::StereoMode == CMOIL::ANAGLYPH_STEREO ) {
		buttonName = "st2";
		if( CMOIL::Main::ColorAnaglyph ) {
			buttonName = "st3";
		}
	}
	ZUI::zuiFindByName( buttonName )->putI( "selected", 1 );
	
	// DisplayMode
	buttonName = "dm0";
	if( parms && parms->dispmode & CMOIL::STICK ) {
		buttonName = "dm1";
	}
	if( parms && parms->dispmode & CMOIL::STICKBALL ) {
		buttonName = "dm2";
	}
	else if( parms && parms->dispmode & CMOIL::SPACEBALL ) {
		buttonName = "dm3";
	}
	ZUI::zuiFindByName( buttonName )->putI( "selected", 1 );
	ZUI::zuiFindByName( "dmhq" )->putI( "selected", CMOIL::Main::DrawBondAsLine==0 );
	ZUI *z = ZUI::zuiFindByName( "ballScale" );
	z->putS( "type", "float" );
	z->putP( "valPtr", &Zmoil_SpaceFillingScale );
	
	// Backbone/Ribbon/Secondary
	ZUI::zuiFindByName( "brs1" )->putI( "selected", parms && (parms->dispmode & CMOIL::BACKBONE) );
	ZUI::zuiFindByName( "brs2" )->putI( "selected", parms && (parms->dispmode & CMOIL::RIBBON) );
	ZUI::zuiFindByName( "brs3" )->putI( "selected", parms && (parms->dispmode & CMOIL::STRUCT2nd) );
	ZUI::zuiFindByName( "brs4" )->putI( "selected", parms && (parms->dispmode & CMOIL::RNABACKBONE) );
	ZUI::zuiFindByName( "brshq" )->putI( "selected", CMOIL::Main::DrawBackboneAsCurve==0 );
	
	// Structure/Movie
	if( parms ) {
		z = ZUI::zuiFindByName( "fps" );
		z->putI( "val", (int) min( 1000.f / max(parms->sleepseconds,1.f), z->getI( "rangeHigh" )) );
		updateIndexSlider();
	}
	
	// Surface
	ZUI::zuiFindByName( "sdm1" )->putI( "selected", parms && (parms->dispmode & CMOIL::SURF) );
	bool raviUI = false;
	if( parms ) {
		char * raviFilename = replaceExtension( parms->pdbName, "surface" );
		if( zWildcardFileExists( raviFilename ) ) {
			raviUI = true;
		}
	}
	ZUI::zuiFindByName( "sdm6" )->putI( "hidden", !raviUI );
	ZUI::zuiFindByName( "sdm7" )->putI( "hidden", !raviUI );
	ZUI::zuiFindByName( "sdm6" )->putI( "selected", parms && (parms->dispmode & CMOIL::RAVINORM) );
	ZUI::zuiFindByName( "sdm7" )->putI( "selected", parms && (parms->dispmode & CMOIL::RAVINORMCOARSE) );

	if( CMOIL::Main::NumImages < 2 ) {
		zMsgQueue( "type=ZUISet key=selected val=0 toZUI=AlignButton" );
		zMsgQueue( "type=ZUISet key=disabled val=1 toZUI=AlignButton" );
	}
	else {
		zMsgQueue( "type=ZUISet key=disabled val=0 toZUI=AlignButton" );
	}
	
	zMsgQueue( "type=ZUILayout" );
}

ZMSG_HANDLER( Zmoil_SyncUI ) {
	syncUIToInputParams();
}

bool modelLoad( char *filespec ) {
	
	// Discover the file format by looking at the extension, and read params.
	ZFileSpec fs( filespec );
	char *ext = fs.getExt();
	if( !stricmp( ext, "xyz" ) ) {
		CMOIL::Main::ReadFileParams( filespec, "-xyz", NULL );
	}
	else if( !stricmp( ext, "in" ) ) {
		CMOIL::Main::ReadFileParams( filespec, "-in", NULL );
	}
	else if( !stricmp( ext, "dcd" ) ) {
		CMOIL::Main::ReadFileParams( filespec, "-dcd", NULL );
	}
	else if( !stricmp( ext, "crd" ) ) {
		CMOIL::Main::ReadFileParams( filespec, "-crd", NULL );
	}
	else /*( !stricmp( ext, "pdb" ) || !stricmp( ext, "ent" ) || ! *ext )*/ {
		CMOIL::Main::ReadFileParams( filespec, "-pdb", NULL );
	}

	
	CMOIL::Main::AddSphereToList();
	// add sphere display lists for any unique sphere sizes for the
	// just loaded model
	
	// Recompute CoM since if a picksphere with cutoff distance was selected,
	// those skipped atoms were included in the first-pass CoM calculations.
	int index = CMOIL::Main::NumImages - 1;
	CMOIL::InputParameters *parms = CMOIL::Main::AAptrs[index] ? &CMOIL::Main::AAptrs[index]->InParm : 0;
	if( parms && parms->picksphere ) {
		zMsgQueue( "type=Zmoil_UpdateCoM when=now" );
	}
	
	// set default render quality based on number of atoms - such that
	// enormous molecules don't render painfully slowly.
	if( CMOIL::Main::AAptrs[index] ) {
		CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[index];
		int viewCount=0;
		for( int i=0; i<CMOIL::Main::AAptrs[index]->numAtom; i++ ) {
			if( !CMOIL::Main::AAptrs[index]->atomList[i].skip ) {
				viewCount++;
			}
		}
		CMOIL::Main::DrawBondAsLine = viewCount > 300;
		CMOIL::Main::DrawBackboneAsCurve = viewCount > 300;
		aa->connectWater();
	}
	
	// send message to reset pickQ
	zMsgQueue( "type=Zmoil_ProgramMode mode=pick" );
	
	// Testing looking at first model and last loaded in the case of
	// multiple models
	if( CMOIL::Main::NumImages > 1 ) {
		CMOIL::Main::ImageIndex = CMOIL::Main::NumImages;
		CMOIL::Main::Overlap1 = 0;
		CMOIL::Main::Overlap1 = CMOIL::Main::NumImages-1;
	}
	else {
		CMOIL::Main::ImageIndex = CMOIL::Main::NumImages-1;
	}
	
	
	
	return true;
}

ZMSG_HANDLER( Zmoil_SelectModel ) {
	void *model = STRING_TO_PTR( msg, model );
	ZUIList *dataList = (ZUIList*)ZUI::zuiFindByName( "dataList" );
	if( dataList ) {
		dataList->setSelectedUserItem( model );
	}
}

int mfmLoad( char *filespec ) {
	// An MFM "model family" file contains a list of related models (all pdb format)
	// along with some attributes that tell us how the models are related.  In
	// the initial impl we expect to see a native struct, a template used for
	// modeling this native struct, and models built using the template.
	// We'll color them so it's obvious what is what.
	int boneQuality = ZUI::zuiFindByName( "brshq" )->getI( "selected" );
		// preserve backbone quality setting across loads

	int autocolor_models = 0;
	int lastWasBestTmpl = 0;
	int templateLoaded = 0;
		// these are for special MFMView settings

	FILE *f = fopen( filespec, "r" );
	if( f ) {
		int count=0;
		ZFileSpec fs( filespec );
		char buf[256], name[64], type[32], tmpl[64], native[64], looppModelTmpl[8];
		float score;
		strcpy( looppModelTmpl, "xyzzy" );
		while( fgets( buf, 256, f ) ) {
			int fields = sscanf( buf, "%s %s", name, type );

			// check for commands before models
			if( !strcmp( name, "autocolor" ) ) {
				autocolor_models = 1;
				continue;
			}

			if( fields < 2 ) {
				continue; 
			}

			ZTmpStr modelFile( "%s%s/%s", fs.getDrive(), fs.getDir(), name );
			if( zWildcardFileExists( modelFile.s ) ) {

				// special view modes for the mfmset
				//
				// 
				if( Zmoil_MFMView == 1) {
					// special for Ron: only load the native, the best loopp, and the first submit

					if( !stricmp( type, "template" ) ) continue;	// no template
					if( strstr( type, "hi2" ) ) continue;			// no best model from other servers
					if( strstr( type, "hi3" ) ) continue;			// no model from best loopp template
					if( !stricmp( type, "model" ) && !strstr( name, "TS1" ) ) continue; // only top submit from loopp
				}
				if( Zmoil_MFMView == 2) {
					// special for thomas: only load the native, the best template, and model created from best template
					int want= 0;
					if( !strcmp( type, "native" ) ) want = 1;	
					else if( strstr( name, "BestTMPL" ) ) {
						want = 1;	
						lastWasBestTmpl = 1;
					}
					else if( lastWasBestTmpl ) {
						want = 1;
						lastWasBestTmpl = 0;
					}
					if( !want ) {
						continue;
					}
				}
				if( Zmoil_MFMView == 3) {
					// special for thomas: only load the native, the first template, and model generated by loopp
					// (for looking at badModelMFM)
					int want = 0;
					if( !strcmp( type, "native" ) ) want = 1;	
					else if( strstr( type, "hi1" ) ) want = 1;	
					else if( strstr( type, "template" ) && !templateLoaded ) {
						want = 1;
						templateLoaded = 1;
					}
					if( !want ) {
						continue;
					}
				}

				trace( "%d. loading %s from mfm...\n", ++count, modelFile.s );
				modelLoad( modelFile.s );
				CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ CMOIL::Main::NumImages - 1 ];
				if( aa ) {
					if( !stricmp( type, "native" ) ) {
						strcpy( native, name );
						aa->properties.putI( "colorIndex", CMOIL::ZC_LtRed );
						for( int c=1; c<aa->numChains; c++ ) {
							aa->DisplayPdbChain( c, 0 );
								// by default only display the first chain
						}
						zMsgQueue( "type=Zmoil_SelectModel model=%ld", long(aa) );
						zMsgQueue( "type=Zmoil_TMAlign __countdown__=2" );
					}
					else if( !stricmp( type, "template" ) ) {
						strcpy( tmpl, name );
						aa->properties.putI( "colorIndex", CMOIL::ZC_LtBlue );
					}
					else if( !strncmp( type, "model", 5 ) ) {
						aa->properties.putS( "native", native );
						aa->properties.putS( "template", tmpl );
						int lp8rank, lp9rank;
						int fields = sscanf( buf, "%*s %*s %f %d %d", &score, &lp8rank, &lp9rank );
						aa->properties.putF( "score", score );
						if( fields == 3 ) {
							aa->properties.putI( "lp8rank", lp8rank );
							aa->properties.putI( "lp9rank", lp9rank );
						}
						int color = CMOIL::ZC_LtGreen;
						if( strstr( type, "_hi1" ) ) {
							// the best model produced by LOOPP
							color = CMOIL::ZC_Yellow;
							strncpy( looppModelTmpl, name, 6 );
							looppModelTmpl[6]=0;
						}
						else if( strstr( type, "_hi2" ) ) {
							// the best model produced by any server, not LOOPP
							color = CMOIL::ZC_Orange;
						}
						else if( strstr( type, "_hi3" ) ) {
							// the LOOPP model produced from the best template
							color = CMOIL::ZC_Yellow;
						}
						else if( strstr( name, looppModelTmpl ) ) {
							aa->properties.putI( "looppBestModel", 1 );
								// indicate this model is also the best model produced by loopp
						}
						aa->properties.putI( "colorIndex", color );
					}
					else {
						trace( "Warning: unknown type %s for model %s in mfm file %s\n", type, name, filespec );
					}
					aa->properties.putS( "type", type );
				}
			}
			else {
				trace( "Error loading %s in mfmLoad of %s\n", modelFile.s, filespec );
			}
		}
		fclose( f );
		operationPending( 4 );
			// indicate an operation has begun that will take 4 frames to complete
	}

	::dataListUpdate();
	syncUIToInputParams();

	zMsgQueue( "type=Zmoil_DisplayMode mode=SKIP val=1 __countdown__=2" );
	zMsgQueue( "type=Zmoil_DisplayMode mode=BACKBONE val=1 __countdown__=2" );
	zMsgQueue( "type=Zmoil_DisplayQuality mode=BACKBONE val=%d __countdown__=2", boneQuality );
	zMsgQueue( "type=Zmoil_SyncUI __countdown__=2" );
	zMsgQueue( "type=ZUISet key=selected val=1 toZUI=AlignButton" );

	return 1;
}


void modelUnload( int index ) {
	assert( index < CMOIL::Main::NumImages );
	assert( CMOIL::Main::AAptrs[ index ] );
	if( CMOIL::Main::AAptrs[ index ] ) {
		delete CMOIL::Main::AAptrs[ index ];
		// @TODO: I have not verified that this does a good job
		// cleaning itself up. (tfb)
		CMOIL::Main::AAptrs[ index ] = 0;
		if( CMOIL::Main::NumImages > 1 ) {
			memcpy( CMOIL::Main::AAptrs+index, CMOIL::Main::AAptrs+index+1, ( CMOIL::Main::NumImages - index - 1 ) * sizeof( CMOIL::Main::AAptrs[0] ) );
		}
		CMOIL::Main::NumImages--;
	}
	
	// Testing looking at first model and last loaded in the case of
	// multiple models
	if( CMOIL::Main::NumImages > 1 ) {
		CMOIL::Main::ImageIndex = CMOIL::Main::NumImages;
		CMOIL::Main::Overlap1 = 0;
		CMOIL::Main::Overlap1 = CMOIL::Main::NumImages-1;
	}
	else {
		CMOIL::Main::ImageIndex = CMOIL::Main::NumImages-1;
	}
}

void modelUnloadAll() {
	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		delete CMOIL::Main::AAptrs[ i ];
		CMOIL::Main::AAptrs[ i ] = 0;
	}
	CMOIL::Main::NumImages = 0;
}

ZTLPVec< char > recentFiles;

void saveRecentFiles() {
	char buffer[1024];
	buffer[0]=0;
	int i;
	for( i=0; i<recentFiles.count; i++ ) {
		if( i>0 ) {
			strcat( buffer, ", " );
		}
		strcat( buffer, recentFiles[i] );
	}
	options.putS( "Zmoil_lastLoaded", buffer );
}

void addRecentFile( char *file ) {
	recentFiles.add( strdup( file ) );
	saveRecentFiles();
}

void removeRecentFile( char *file ) {
	for( int i=0; i<recentFiles.count; i++ ) {
		if( !strcmp( recentFiles[i], file ) ) {
			recentFiles.del( i );
			break;
		}
	}
	saveRecentFiles();
}

void clearRecentFiles() {
	recentFiles.clear();
	saveRecentFiles();
}





ZTLPVec< char > setList;
int setListIndex = -1;
int setListLoad( char *filespec ) {
	// A "set" is a file which contains a list of files in a given order.
	// It is used to allow the user to use "Next" and "previous" UI that show 
	// up only when a set is loaded.  The next and previous case the next/previous
	// file in the list to be loaded (and the previous unloaded).
	//
	// Filenames given in the .set file should be relative to the location of
	// the .set file.
	//
	// This was originally written for use in conjunction with mfm files:
	// I want to present a set of model-families, and allow the user to 
	// easily navigate through them.  I would like them presented in a particular
	// order, which is given by the order in the file.
	//
	//
	
	FILE *f = fopen( filespec, "r" );
	if( !f ) { 
		return 0;
	}
	ZFileSpec fs( filespec );

	setList.clear();
	setListIndex = -1;

	char buf[256];
	while( fgets( buf, 256, f ) ) {
		char *ptr = buf + strlen( buf ) - 1;
		while( ptr >= buf && isspace( *ptr ) ) {
			*ptr-- = 0x0;
		}
		ZTmpStr fname( "%s/%s", fs.getDir(), buf );
		if( zWildcardFileExists( fname.s ) ) {
			setList.add( strdup( fname.s ) );
		}
		else {
			trace( "WARNING: unable to locate %s specified in set file %s\n", fname.s, filespec );
		}
	}
	
	fclose( f );

	if( setList.count > 0 ) { 
		ZUI::zuiFindByName( "setListNext" )->putI( "hidden", 0 );
		ZUI::zuiFindByName( "setListPrev" )->putI( "hidden", 0 );
		ZUI::zuiFindByName( "setListClear" )->putI( "hidden", 0 );
		ZUI::zuiFindByName( "setListIndex" )->putI( "hidden", 0 );
		ZUI::zuiFindByName( "fileOpen" )->putI( "hidden", 1 );

		zMsgQueue( "type=Zmoil_SetListNext" );
	}
	model.viewInfo.putS( "setListName", filespec );
	return 1;
}

ZMSG_HANDLER( Zmoil_SetListPrev ) {
	setListIndex--;
	if( setListIndex < 0 ) {
		setListIndex = setList.count-1;
	}
	ZUI *z = ZUI::zuiFindByName( "setListIndex" );
	z->putS( "text", ZTmpStr( "Index %d of %d", setListIndex+1, setList.count ) );
	zMsgQueue( "type=Zmoil_Open filespec='%s' setListIndex=%d", setList[ setListIndex ], setListIndex );
		// pass setListIndex so it can be saved in recent options if load is successful
}

ZMSG_HANDLER( Zmoil_SetListNext ) {
	setListIndex++;
	if( setListIndex >= setList.count ) {
		setListIndex = 0;
	}
	ZUI *z = ZUI::zuiFindByName( "setListIndex" );
	z->putS( "text", ZTmpStr( "Index %d of %d", setListIndex+1, setList.count ) );
	zMsgQueue( "type=Zmoil_Open filespec='%s' setListIndex=%d", setList[ setListIndex ], setListIndex );
		// pass setListIndex so it can be saved in recent options if load is successful
}

ZMSG_HANDLER( Zmoil_SetListClear ) {
	// clear the list and remove UI for sets
	model.viewInfo.del( "setListName" );
	setList.clear();
	modelUnloadAll();
	clearRecentFiles();
	ZUI::zuiFindByName( "setListNext" )->putI( "hidden", 1 );
	ZUI::zuiFindByName( "setListPrev" )->putI( "hidden", 1 );
	ZUI::zuiFindByName( "setListClear" )->putI( "hidden", 1 );
	ZUI::zuiFindByName( "setListIndex" )->putI( "hidden", 1 );
	ZUI::zuiFindByName( "fileOpen" )->putI( "hidden", 0 );
	::dataListUpdate();
}



ZHashTable recentOpt;
void savePathOptions() {
	recentOpt.clear();
	recentOpt.putS( "Zmoil_openPath", options.getS( "Zmoil_openPath" ) );
	if( options.has( "Zmoil_lastLoaded" ) ) {
		recentOpt.putS( "Zmoil_lastLoaded", options.getS( "Zmoil_lastLoaded" ) );
	}
	recentOpt.putI( "Zmoil_skipTutorial", options.getI( "Zmoil_skipTutorial", 0 ) );
	zConfigSaveFile( "zmoilrecent.opt", recentOpt );
}
void loadPathOptions() {
	zConfigLoadFile( "zmoilrecent.opt", recentOpt );
	options.putS( "Zmoil_openPath", recentOpt.getS( "Zmoil_openPath", "" ) );
	if( recentOpt.has( "Zmoil_lastLoaded" ) ) {
		options.putS( "Zmoil_lastLoaded", recentOpt.getS( "Zmoil_lastLoaded" ) );
	}
	options.putI( "Zmoil_skipTutorial", recentOpt.getI( "Zmoil_skipTutorial", 0 ) );
}
void updateDefaultPaths( char *filespec, char *key ) {
	// CREATE & STORE a path with * wildcard to save user preferred location of loading/saving files.
	if( zWildcardFileExists( filespec ) ) {
		ZFileSpec fs ( filespec );
		char *drive = fs.getDrive();
		char *dir = fs.getDir();
		char *file = fs.getFile( 0 );
		char *ext = fs.getExt();
		ZFileSpec fs2 ( zFileSpecMake(   FS_DRIVE, drive, FS_DIR, dir, FS_FILE, "*", FS_END ) );
			// new: do not store the extension.
		options.putS( key, fs2.path );
	}
}

void setupCmoil() {
	// Build hashtable to map display mode names to enum values
	model.modeNameToValue.putI( "SKIP", CMOIL::SKIP );
	model.modeNameToValue.putI( "SPACEBALL", CMOIL::SPACEBALL );
	model.modeNameToValue.putI( "STICK", CMOIL::STICK );
	model.modeNameToValue.putI( "RIBBON", CMOIL::RIBBON );
	model.modeNameToValue.putI( "BACKBONE", CMOIL::BACKBONE );
	model.modeNameToValue.putI( "RNABACKBONE", CMOIL::RNABACKBONE );
	model.modeNameToValue.putI( "UNITBALL", CMOIL::UNITBALL );
	model.modeNameToValue.putI( "STICKBALL", CMOIL::STICKBALL );
	model.modeNameToValue.putI( "STRUCT2nd", CMOIL::STRUCT2nd );
	model.modeNameToValue.putI( "SURF", CMOIL::SURF );
	model.modeNameToValue.putI( "SUBSURF", CMOIL::SUBSURF );
	model.modeNameToValue.putI( "INSPHERE", CMOIL::INSPHERE );
	model.modeNameToValue.putI( "RAVINORM", CMOIL::RAVINORM );
	model.modeNameToValue.putI( "RAVINORMCOARSE", CMOIL::RAVINORMCOARSE );
	model.modeNameToValue.putI( "RESIDUEPICK", CMOIL::RESIDUEPICK );
	model.modeNameToValue.putI( "REDISPLAY", CMOIL::REDISPLAY );
	model.modeNameToValue.putI( "NONE_STEREO", CMOIL::NONE_STEREO );
	model.modeNameToValue.putI( "HARDWARE_STEREO", CMOIL::HARDWARE_STEREO );
	model.modeNameToValue.putI( "ANAGLYPH_STEREO", CMOIL::ANAGLYPH_STEREO );

	// Build hashtable to map 3 letter amino abbrev to 1 letter codes.
	// I grabbed this data from the tmalign fortran code. (Tfb)
	// todo: PUT THIS IN THE AMINOACID CLASS OF CMOIL
	model.aminoNameToChar.putS( "BCK", "X" );
	model.aminoNameToChar.putS( "GLY", "G" );
	model.aminoNameToChar.putS( "ALA", "A" );
	model.aminoNameToChar.putS( "SER", "S" );
	model.aminoNameToChar.putS( "CYS", "C" );
	model.aminoNameToChar.putS( "VAL", "V" );
	model.aminoNameToChar.putS( "THR", "T" );
	model.aminoNameToChar.putS( "ILE", "I" );
	model.aminoNameToChar.putS( "PRO", "P" );
	model.aminoNameToChar.putS( "MET", "M" );
	model.aminoNameToChar.putS( "ASP", "D" );
	model.aminoNameToChar.putS( "ASN", "N" );
	model.aminoNameToChar.putS( "LEU", "L" );
	model.aminoNameToChar.putS( "LYS", "K" );
	model.aminoNameToChar.putS( "GLU", "E" );
	model.aminoNameToChar.putS( "GLN", "Q" );
	model.aminoNameToChar.putS( "ARG", "R" );
	model.aminoNameToChar.putS( "HIS", "H" );
	model.aminoNameToChar.putS( "PHE", "F" );
	model.aminoNameToChar.putS( "TYR", "Y" );
	model.aminoNameToChar.putS( "TRP", "W" );
	model.aminoNameToChar.putS( "CYX", "C" );
	
	/*
	data aa/ 'BCK','GLY','ALA','SER','CYS','VAL','THR','ILE',
     &     'PRO','MET','ASP','ASN','LEU',
     &     'LYS','GLU','GLN','ARG',
     &     'HIS','PHE','TYR','TRP','CYX'/
      character*1 slc(-1:20)
      data slc/'X','G','A','S','C','V','T','I',
     &     'P','M','D','N','L','K','E','Q','R',
     &     'H','F','Y','W','C'/
	*/

}



	
//-----------------------------------------------------------------------------
// 'Main' Message Handlers

ZMSG_HANDLER( Zmoil_ToggleUI ) {
	// tfb: I've played with a couple options here.  The most obvious thing to do is hide the root UI 
	// node, which I did first.  The downside to this is that keyBindings for UI that is hidden are
	// not processed: so you can't control the program with the UI hidden!  So I tried moving the UI
	// offscreen.  This works, and key bindings work: interestingly, when a keyBinding is executed,
	// the UI comes back onscreen; surely this is because my hacked method of moving it offscreen 
	// is not quite right.  But the bigger downside to this approach is that either the visible
	// performance gain to be had in hiding the UI is not had when just moving it offscreen.  This makes
	// sense: it isn't the triangles that are costing us with the UI, but all the logic.
	//
	assert( g_uiPanel );
	if( g_uiPanel->getI( "hidden" ) ) {
		// MAKE Visible
		zMsgQueue( "type=ZUISet key=hidden val=0 toZUI=%s", g_uiPanel->name );
		zMsgQueue( "type=ZUISet key=layoutManual_x val=%d toZUI=pluginPanel", UI_PANEL_WIDTH );
		zMsgQueue( "type=ZUISet key=layoutManual_w val='W %d -' toZUI=pluginPanel; type=ZUILayout toZUI=pluginPanel", UI_PANEL_WIDTH );
		// SYNC some UI state
		zMsgQueue( "type=ZUISet key=selected val=%d toZUI=movieMode", CMOIL::Main::MovieMode );
		zMsgQueue( "type=ZUISet key=selected val=%d toZUI=pickMode", CMOIL::Main::PickMode );
		zMsgQueue( "type=ZUISet key=selected val=%d toZUI=rockingMode", CMOIL::Main::Rocking );
		tempMessage( "" );
	}
	else {
		// HIDE: 
		zMsgQueue( "type=ZUISet key=hidden val=1 toZUI=%s", g_uiPanel->name );
		zMsgQueue( "type=ZUISet key=layoutManual_x val=0 toZUI=pluginPanel" );
		zMsgQueue( "type=ZUISet key=layoutManual_w val='W 0 -' toZUI=pluginPanel; type=ZUILayout toZUI=pluginPanel");
		tempMessage( "F8 to show UI", 119, 60 );
	}
}

ZMSG_HANDLER( Zmoil_Open ) {
	// opens filespec if passed, else creates filepicker to select file
	//if( CMOIL::Main::NumImages > 9 ) {
	//	zuiMessageBox( "Max Models Loaded", "Please unload one of the models before attempting to load another.", 1 );
	//	return;
	//}
	char cwdbuf[256];
	getcwd( cwdbuf, 256 );
	if( zmsgHas(filespec) ) {
		char filespec[1024];
		strcpy( filespec, zmsgS(filespec) );
		char *fname = strtok( filespec, "," );
		while( fname ) {
			while( isspace( *fname ) ) { fname++; }
			trace( "cwd=%s, Will attempt to load the zmoil file: %s\n", cwdbuf, fname );
			if ( zWildcardFileExists( fname ) ) {
				int loaded=0;
				ZFileSpec fs( fname );
				if( !stricmp( fs.getExt(), "mfm" ) ) {
					modelUnloadAll();
					if( ! zmsgHas( setListIndex ) ) {
						clearRecentFiles();
							// if navigating a set, don't clear the .set file from recent list
					}
					loaded = mfmLoad( fname );
				}
				else if( !stricmp( fs.getExt(), "set" ) ) {
					modelUnloadAll();
					clearRecentFiles();
					loaded = setListLoad( fname );
				}
				else {
					loaded = modelLoad( fname );
				}
				if( loaded ) {
					trace( "%s loaded.\n", fname );
					if( !zmsgHas( setListIndex ) ) {
						// only add to recent files if not navigating a set
						updateDefaultPaths( fname, "Zmoil_openPath" );
						addRecentFile( fname );
					}
					model.viewInfo.putS( "filename", filespec );
				}
				else {
					// modelLoad failed.  Display notice to user & reset UI
					zuiMessageBox( "Load Error", "The file you have selected does not appear to be in a supported format.", 1 );
					break;
				}
			}
			else {
				// bad filename was passed; last model remains open
				trace( "The file passed does not exist.\n" );
			}
			fname = strtok( NULL, "," );
		}
		::dataListUpdate();
		syncUIToInputParams();

		zMsgQueue( "type=Zmoil_ModelDirty dirty=0 nosave=1" );
		zMsgQueue( "type=Zmoil_DisplayControl cmd=reset" );
	}
	else {
		zuiChooseFile( "Open File", options.getS( "Zmoil_openPath" ), 1, "type=Zmoil_Open", 0, true );
	}
}

//-----------------------------------------------------------------------------
// 'Models' Message Handlers

ZMSG_HANDLER( Zmoil_DisplayModels ) {
	setAllListItems( "dataList", zmsgI( all ) );
}

ZMSG_HANDLER( Zmoil_ModelUnload ) {
	CMOIL::Main::PickQ.reset();
	int unloadAll=zmsgI( all );
	if( unloadAll ) {
		manualAlign = 0;
		zMsgQueue( "type=ZUISet toZUI=manualAlign key=selected val=0" );
		modelUnloadAll();
	}
	else {
		CMOIL::AtomArray *aa = modelListGetSelected( "dataList", 0, 0 );
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			if( CMOIL::Main::AAptrs[i] == aa ) {
				if( manualAlign == CMOIL::Main::AAptrs[i] ) {
					manualAlign = 0;
					zMsgQueue( "type=ZUISet toZUI=manualAlign key=selected val=0" );
				}
				removeRecentFile( aa->InParm.fileName );
				modelUnload( i );
			}
		}
	}
	if( unloadAll || CMOIL::Main::NumImages == 0 ) {
		clearRecentFiles();
	}
	dataListUpdate();
}

ZMSG_HANDLER( Zmoil_AutoColor ) {
	int colors = zmsgI( colors );

	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ i ];
		if( colors ) {
			aa->properties.putI( "colorIndex", i % CMOIL::ZC_NumColors );
		}
		else {
			aa->properties.del( "colorIndex" );
		}
	}
	dataListUpdate();
	ZUI::dirtyAll();
}

struct CrdColorDef {
	double origin;
	double binSize;
	CMOIL::AtomArray *crd;
};
ZHashTable crdColorDefinitions;
#define NUM_CRD_COLORRANGE 5


void setColorsFromCrdColorDefinitions() {
	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ i ];
		if( aa->InParm.filetag == CMOIL::fCRD || aa->InParm.filetag == CMOIL::fDCD ) {
			CrdColorDef *ccd = (CrdColorDef*)crdColorDefinitions.bgetp( &aa->numAtom, sizeof( aa->numAtom ) );
			if( ccd && ccd->crd ) {
				for( int i=0; i<ccd->crd->numAtom; i++ ) {
					int bin = 0;
					if( ccd->binSize > 0.0 ) {
						bin = int( ( ccd->crd->atomList[ i ].weight - ccd->origin ) / ccd->binSize );
						bin = min( bin, NUM_CRD_COLORRANGE-1 );
					}
					aa->atomList[i].color = CMOIL::zmoilColors[ bin ];
				}
			}
			else {
				aa->setColor();
					// sets default colors
			}
		}
	}
}

ZTmpStr gBoxColorFile;
ZMSG_HANDLER( Zmoil_BoxColor ) {
	if( zmsgI( newfile ) ) {
		if( CMOIL::AtomArray::readBoxInfo( gBoxColorFile ) ) {
			Zmoil_gridColor = 1;
		}
		else {
			Zmoil_gridColor = 0;
			gBoxColorFile.set( "" );
			zMsgQueue( "type=ZUISet key=selected val=0 toZUI=BoxColor" );
		}
	}
	else {
		Zmoil_gridColor = zmsgI( colors );
		if( Zmoil_gridColor && ! zWildcardFileExists( gBoxColorFile ) ) {
			zMsgQueue( "type=Zmoil_ChooseBoxColorFile" );
			return;
		}
	}

	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ i ];

		// TODO: check that there is some boxInfo data available for this structure.
		// For now, we are just coloring in a range as a test.

		aa->setColor();
	}
}

ZMSG_HANDLER( Zmoil_ChooseBoxColorFile ) {
	zuiChooseFile( "Choose Box Color File", options.getS( "Zmoil_openPath" ), 1, "type=Zmoil_SetBoxColorFile" );
}

ZMSG_HANDLER( Zmoil_SetBoxColorFile ) {
	char *filespec = zmsgS( filespec );
	if( zWildcardFileExists( filespec ) ) {
		gBoxColorFile.set( filespec );
		ZFileSpec fs( filespec );
		zMsgQueue( "type=ZUISet key=text toZUI=boxcolorfile val='BoxColor: %s'", escapeQuotes( fs.getFile() ) );
		zMsgQueue( "type=Zmoil_BoxColor newfile=1" );
	}
}

ZMSG_HANDLER( Zmoil_CrdColor ) {
	// Idea to use the weighting factor column from a CRD file to define 
	// the coloring of atoms.  For now we will take the first structure
	// that is a CRD, and split the range of values seen into 5 subranges.
	// Each subrange will be assigned a unique color.
	//
	// We will walk through the model list, and for any DCD structures which are
	// found to contain the same number of atoms, we will assume they should
	// be colored similarly, allowing a CRD to define the colors for a DCD

	crdColorDefinitions.clear();

	if( zmsgI( colors ) ) {
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ i ];
			if( aa->InParm.filetag == CMOIL::fCRD && !crdColorDefinitions.bgetp( &aa->numAtom, sizeof( aa->numAtom ) ) ) {
				// get range of values, split into NUM_CRD_COLORRANGE
				double minW=1e100, maxW=-1e100;
				for( int i=0; i<aa->numAtom; i++ ) {
					double w = aa->atomList[i].weight;
					minW = min( w, minW );
					maxW = max( w, maxW );
				}
				CrdColorDef *ccd = new CrdColorDef();
				ccd->crd = aa;
					// the atom weights in this crd will be referenced by all crd/dcd structs with same numAtom
				ccd->binSize = ( maxW - minW ) / NUM_CRD_COLORRANGE; 
				ccd->origin = minW;
				crdColorDefinitions.bputP( &aa->numAtom, sizeof( aa->numAtom ), ccd, 1, 1 );
			}
		}
	}
	
	setColorsFromCrdColorDefinitions();
}


ZMSG_HANDLER( Zmoil_TMAlign ) {
	// We'll need some way to specifiy which model is to be aligned to which target.
	// For now we'll align the last loaded model to the first one.  These are the two
	// that are currently displayed.

	CMOIL::AtomArray *aa0,*aa1;
	aa0 = modelListGetSelected( "dataList" ); 
	if( !aa0 ) {
		zuiMessageBox( "Alignment", "Please select the structure you wish to align the other structures to.\n\nThe selected structure will appear highlighted (in green) within the list once it has been selected.", 1  );
		return;
	}

	aa0->clearAlignments();
	aa0->resetCrd();
	aa0->computeCentroid();
	aa0->properties.putI( "tmalignMaster", 1 );

	bool success = true;
	bool badformat = false;
	int trivial = ZUI::zuiFindByName( "tm_trivial" )->getI( "selected" );

	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		aa1 = CMOIL::Main::AAptrs[ i ];
		if( aa1 != aa0  ) {
			aa1->properties.del( "tmalignMaster" );
			aa1->clearAlignments();
			aa1->resetCrd();
			aa1->computeCentroid();
			
			if( ( aa0->InParm.filetag == CMOIL::fCRD || aa0->InParm.filetag == CMOIL::fPDB ) &&
				( aa1->InParm.filetag == CMOIL::fCRD || aa1->InParm.filetag == CMOIL::fPDB ) ) {
				success = success && aa1->setAlignmentTarget( aa0, trivial );
			}
			else {
				badformat = true;
				break;
			}
		}
	}

	if( badformat ) {
		zuiMessageBox( "TM-Align Error", "Alignment may only be made on PDB or CRD models.", 1 );
	}
	else if( !success ) {
		zuiMessageBox( "TM-Align Error", ZTmpStr("Alignment was not successful.  Does %stmalign program exist?  If you are launching zmoil other than through moil.tcl or the zmoil folder, you need to ensure the tmalign program is in your execution path.  If you are aligning CRD files, the program crd2pdb must also be in your execution path.", CMOIL::Main::ExeDir ), 1 );
	}
	if( !ZUI::zuiFindByName( "comNever" )->getI( "selected" ) ) {
		zMsgQueue( "type=Zmoil_UpdateCoM when=now" );
	}
}

ZMSG_HANDLER( Zmoil_ManualAdjust ) {
	if( zmsgI( val ) ) {
		manualAlign = modelListGetSelected( "dataList" );
	}
	else {
		manualAlign = 0;
	}
}

ZMSG_HANDLER( Zmoil_ClearAlignment ) {
	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[i];
		if( aa ) {
			aa->clearAlignments();
		}
	}
	if( !ZUI::zuiFindByName( "comNever" )->getI( "selected" ) ) {
		zMsgQueue( "type=Zmoil_UpdateCoM when=now" );
	}
}

ZMSG_HANDLER( Zmoil_WriteTransformedPDB ) {
	for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[i];
		if( aa && aa->InParm.filetag == CMOIL::fPDB) {
			aa->WriteTransformedPdb();
			// @TODO: there exist writecrd and writepdb we can use that currently
			// write the org coordinates; we can add flags to output the transformed
			// coords instead.
		}
		else {
			zuiMessageBox( "Not Supported", "Some models were not written.  This function is only supported for models loaded from a native PDB file format.", 1 );
		}
	}
}

//-----------------------------------------------------------------------------
// Other Message Handlers

ZMSG_HANDLER( Zmoil_TutorialOK ) {
	ZUI *z = ZUI::zuiFindByName( "skipTutorial" );
	options.putI( "Zmoil_skipTutorial", z->getI("selected") );
	savePathOptions();
}

ZMSG_HANDLER( Zmoil_ProgramMode ) {
	char *mode = zmsgS( mode );
	int   val  = zmsgI( val );

	if( !strcmp( mode, "pick" ) )  {
		// experimental: always in pickmode: here we just reset the queue
		CMOIL::Main::PickMode = 1;
		if( CMOIL::Main::AAptrs[0] ) {
			if( NO_OVERLAP ) {
				CMOIL::Main::PickQ.reset( CMOIL::Main::ImageIndex );
						// erase selected atom for different image
			}
			else {
				CMOIL::Main::PickQ.reset( CMOIL::Main::Overlap1, CMOIL::Main::Overlap2 );
			}
		}
		tempMessage( "Reset picked atoms.  Left-click to pick..." );
		CMOIL::Main::MsgBrd.erase();
	}
	else if( !strcmp( mode, "movie" ) ) {
		CMOIL::Main::MovieMode = val;
		if( val ) {
			CMOIL::Main::Rocking=0;
			zMsgQueue( "type=ZUISet key=selected val=0 toZUI=rockingMode" );
		}
	}
	else if( !strcmp( mode, "rock" ) ) {
		CMOIL::Main::Rocking = val;
		if( val ) {
			CMOIL::Main::MovieMode=0;
			zMsgQueue( "type=ZUISet key=selected val=0 toZUI=movieMode" );
		}
	}
	else if( !strcmp( mode, "allowPPFlip" ) ) {
		Zmoil_PeptidePlaneFlipping = val;
			// allow the plane to flip as dictated by CA-O backbone atoms;
			// this defeats the "normal inversion" that occurs in cmoil_atom.cpp to prevent unsightly twists
			// in the backbone. (search for Zmoil_PeptidePlaneFlipping)
	}
	else if( !strcmp( mode, "atomVelocity" ) ) {
		Zmoil_VelocityDisplay = val;
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ i ];
			if( aa && aa->velfp && !aa->properties.has( "velocityMax" ) ) {
				//aa->getVelocityStats();
				aa->ReadStruct( aa->currentmol - 1 );
					// cause the velocities for the current  frame to be loaded
			}
		}
	}
}

ZMSG_HANDLER( Zmoil_MovieParam ) {
	char *param = zmsgS( param );
	if( !strcmp( param, "fps" ) ) {
		int val  = zmsgI( val );
		if( val > 0 ) {
				// val in frames per sec
			float msPerFrame = 1000.f / val;
			for( int i=0; i<CMOIL::ARRAY_MAXNUM; i++ ) {
				if( CMOIL::Main::AAptrs[i] ) {
					CMOIL::Main::AAptrs[i]->InParm.sleepseconds = msPerFrame;
				}
			}
		}
	}
}

ZMSG_HANDLER( Zmoil_PickInfo ) {
	char *info = zmsgS( info );

	if( !strcmp( info, "angle" ) ) {
		CMOIL::Main::PickQ.print_angle();
	}
	else if( !strcmp( info, "distance") ) {
		CMOIL::Main::PickQ.print_distance();
	}
	else if( !strcmp( info, "torsion" ) ) {
		CMOIL::Main::PickQ.print_torsional();
	}
}


//-----------------------------------------------------------------------------
// Image Capture Message Handlers

ZMSG_HANDLER( Zmoil_ImageType ) {
	model.viewInfo.putS( "captureFile", replaceExtension( model.viewInfo.getS( "captureFile", "zmoil" ), zmsgS(itype) ));
	ZUI::zuiFindByName( "captureFilename" )->putS( "text", model.viewInfo.getS( "captureFile" ) );
}

ZMSG_HANDLER( Zmoil_SetImageCaptureFilename ) {
	bool bLaunchFilePicker = false;
	if( !zmsgI( cancel ) ) {
		if( zmsgHas(filespec) ) {
			// VALIDATE/SET captureFile
			char *fs = zmsgS( filespec );
			ZFileSpec fspec( fs );
			char *ext = fspec.getExt();
			if( !strcmp( ext, "png" ) || !strcmp( ext, "tif" ) || !strcmp( ext, "mpeg" ) || !strcmp( ext, "gif" ) ) {
				model.viewInfo.putS( "captureFile", fs );
				ZUI::zuiFindByName( "captureFilename" )->putS( "text", model.viewInfo.getS( "captureFile" ) );
				ZUI::zuiFindByName( ext )->putI( "selected", 1 );
				zMsgQueue( "type=ZUILayout" );
			}
			else {
				// INVALID format chosen (as spec'd by extension): inform user 
				char *msg = "Zmoil supports png, tif, mpeg and gif extensions for image capture.  Please use a filename with one of these extensions.";
				char *title = "Unsupported Save Format";
				zuiMessageBox( title, msg, 1 );
				bLaunchFilePicker = true;
			}
		}
		else {
			bLaunchFilePicker = true;
		}
	}
	if( bLaunchFilePicker ) {
		// LAUNCH picker
		char *openPath = options.getS( "Zmoil_openPath" );
		ZFileSpec captureFile( model.viewInfo.getS( "captureFile", "zmoil.png" ) );
		char *imagePath = replaceExtension( openPath, captureFile.getExt() );
		zuiChooseFile( "Set Capture Filename (png, tif, mpeg, gif)", imagePath, 0, "type=Zmoil_SetImageCaptureFilename" );
	}
}

ZMSG_HANDLER( Zmoil_SaveImage ) {
	// Saves image to captureFile as spec'd in model.viewInfo
	// Format of saved image is determined by the extension of the filename.

	zGraphFileInit();
		// ensure freeimage is ready to go.

	char *capFile = model.viewInfo.getS( "captureFile" );
	assert( capFile );
	ZFileSpec capFS( capFile );
	char *ext = capFS.getExt();
	assert( ( !strcmp( ext, "png" ) || !strcmp( ext, "tif" ) || !strcmp( ext, "mpeg" ) || !strcmp( ext, "gif" ) ) );

	// @TODO:
	// We may wish to allow the specification of a begin frame and end frame, such that the user
	// may select to save a series of still images, and optionally encode to mpeg.  For now, choosing
	// a still-format means capture the current image, and choosing the mpeg format means
	// save a series of still images (1 for each structure index) and then encode to mpeg.

	if(  ( !strcmp( ext, "png" ) || !strcmp( ext, "tif" ) /*|| !strcmp( ext, "gif" )*/ ) ) {
		// Save current image:
		capturePluginImage( capFile );
	}
	else if ( !strcmp( ext, "mpeg" ) ) {
		// We will capture multiple .png images.  This happens in two instances:
		// 1. a dcd file with multiple frames, for which we will capture 1 png per frame
		// 2. a setlist, for which we will capture 1 png for each file in the setlist

		//
		// Create a the folder that will store the individual frames,
		// and create a basename for the files that will be saved.
		//
		char *capFolder;
		if( setListIndex >= 0 ) {
			capFolder = replaceExtension( model.viewInfo.getS( "setListName", capFile ), "frames" );
		}
		else {
			replaceExtension( capFile, "frames" );
		}
		model.viewInfo.putS( "captureMovieFolder", capFolder );
		char cmdFolder[256];
		strcpy( cmdFolder, capFolder );
		#ifdef WIN32
			char *cmd = "del /F /Q %s\\*";
			char *p = cmdFolder;
			while (*p) {
				if( *p == '/' ) {
					*p = '\\';
				}
				p++;
			}
		#else
			char *cmd = "rm -f %s/*";
		#endif
		
		ZTmpStr zcmd(  cmd, cmdFolder );
		system( zcmd.s );
		#ifdef WIN32
			 mkdir( capFolder );
		#else
			 mkdir( capFolder, S_IRWXO | S_IRWXU | S_IRWXG );
		#endif
		strcat( capFolder, "/" );
		if( setListIndex < 0 ) {
			strcat( capFolder, capFS.getFile( 0 ) );
		}
		model.viewInfo.putS( "captureMovieBasename", capFolder );


		//
		// Go to start Index
		//
		if( setListIndex >= 0 ) {
			// we'll just start capturing fromt the current index
		}
		else {
			CMOIL::AtomArray *aaptr = CMOIL::Main::AAptrs[ CMOIL::Main::ImageIndex ];
			aaptr->ReadStruct( 0 );
		}


		// Temp: since we do not actually encode the mpeg movie for now, let the
		// user know that we will capture a series of frames and that he can use
		// ffmpeg, convert, mpeg2encode etc... to generate movies.
		char *msg = "Zmoil will write a series of .png frames to\n\n  %s\n\nYou may then use any program (such as ffmpeg, convert, mpeg2encode, QuickTime, etc.) to create a movie in the format of your choice.";
		ZTmpStr message( msg, model.viewInfo.getS( "captureMovieFolder" ) );
		char *title = "Movie Creation";
		zuiMessageBox( title, message.s, 0, "type=Zmoil_EnableMultiframeCapture __countdown__=1", "" );
			// note use of non-null cancel message to cause cancel button to appear
			// use __countdown__ to allow time for dialog box to remove from screen
	}
	else if( !strcmp( ext, "gif" ) ) {
		// Multi-frame GIF:

		// Go to start Index
		CMOIL::AtomArray *aaptr = CMOIL::Main::AAptrs[ CMOIL::Main::ImageIndex ];
		aaptr->ReadStruct( 0 );

		// Queue message to enable Multiframe capture
		zMsgQueue( "type=Zmoil_EnableMultiframeCapture" );
	}
}

ZMSG_HANDLER( Zmoil_EnableMultiframeCapture ) {
	// Set Image-Capture Flag
	model.viewInfo.putI( "captureFrames", 1 );
}

//-----------------------------------------------------------------------------
// Display Options Message Handlers

ZMSG_HANDLER( Zmoil_DisplayGeneral ) {
	char *mode = zmsgS( mode );

	if( !strcmp( mode, "WBACKGRND" ) ) {
		int   val  = zmsgI( val );
		if( CMOIL::Main::AAptrs[0] ) {
			CMOIL::Main::AAptrs[0]->InParm.bg.r = val ? 0xFF : 0x0;
			CMOIL::Main::AAptrs[0]->InParm.bg.g = val ? 0xFF : 0x0;
			CMOIL::Main::AAptrs[0]->InParm.bg.b = val ? 0xFF : 0x0;
		}
	}
	else if ( !strcmp( mode, "SHINY" ) ) {
		CMOIL::Main::Shine = zmsgF( val );	
	}
	else if ( !strcmp( mode, "FOG" ) ) {
		gEnableFog = zmsgI( val );
	}
	else if ( !strcmp( mode, "AXES" ) ) {
		gEnableAxes = zmsgI( val );
	}
	else if ( !strcmp( mode, "SINDEX" ) ) {
		gEnableStructInfoDisplay = zmsgI( val );
	}
	else if ( !strcmp( mode, "SEQ" ) ) {
		ZUI *main, *resi;
		main = ZUI::zuiFindByName( "cmoilRender" );
		resi = ZUI::zuiFindByName( "residuePanel" );
		if( main && resi ) {
			int showResi = zmsgI( val );
			//main->putI( "layoutManual_y", showResi ? (int)resi->h : 0 );
			main->putS( "layoutManual_h", showResi ? (char*)ZTmpStr("H %d -", int( resi->h ) ) : (char*)"H" );
			resi->putI( "hidden", !showResi );
			zMsgQueue( "type=ZUILayout" );
		}
	}
	else if ( !strcmp( mode, "WATER" ) ) {
		Zmoil_DisplayWater = zmsgI( val );
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::Main::AAptrs[i]->setSkip();
		}
	}
	else if ( !strcmp( mode, "HBONDS" ) ) {
		Zmoil_DisplayHBonds = zmsgI( val );
		//for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
		//	if( ! CMOIL::Main::AAptrs[i]->properties.getI( "hBondsCalculated" ) ) {
		//		CMOIL::Main::AAptrs[i]->calcHydrogenBonds();
		//	}
		//}
			// the above is unnecessary, since we will compute them on demand
			// for any AtomArray that is rendered with hbonds on; the downside is
			// that if an AA is turned off, the delay will occur when the user turns
			// it on.
	}

}

ZMSG_HANDLER( Zmoil_HBondDistance ) {
	if( !zmsgI( twiddleComplete ) ) {
		Zmoil_HBondDistance = zmsgF( val );
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::Main::AAptrs[ i ]->properties.del( "hBondsCalculated" );
		}
	}
}



ZMSG_HANDLER( Zmoil_DisplayMode ) {
	char *mode = zmsgS( mode );
	int   val  = zmsgI( val );
	int   modeVal = model.modeNameToValue.getI( mode );

	int polymode = 0;
	if( modeVal == CMOIL::RIBBON || modeVal == CMOIL::STRUCT2nd ) {
		polymode = GL_FILL;
	}

	// SETUP clearmask to handle mutually exclusive display modes
	int clearMask = 0;
	if( val && (modeVal==CMOIL::SKIP || modeVal == CMOIL::STICK || modeVal == CMOIL::STICKBALL || modeVal == CMOIL::SPACEBALL) ) {
		clearMask = CMOIL::STICK | CMOIL::STICKBALL | CMOIL::SPACEBALL;
	}
	else if ( val && (modeVal == CMOIL::RIBBON || modeVal == CMOIL::STRUCT2nd) ) {
		clearMask = CMOIL::RIBBON | CMOIL::STRUCT2nd;
	}

	CMOIL::AtomArray *aaptr;
	int ImageIndex = CMOIL::Main::ImageIndex;
	int Overlap1   = CMOIL::Main::Overlap1;
	int NumImages  = CMOIL::Main::NumImages;
	for( int i=0; i<NumImages; i++ ) {
		// read at most two structures
		aaptr = CMOIL::Main::AAptrs[i];
		if ( polymode != 0 ) {
			aaptr->InParm.ribbonpolymode=polymode;
		}
		if( val ) {
			// clear any mutually exclusive bits & turn selected bit on
			aaptr->InParm.dispmode &= ~clearMask;
			aaptr->InParm.dispmode |= modeVal;
		}
		else {
			// turn bit off
			aaptr->InParm.dispmode &= ~modeVal;
		}
		aaptr->ReProcessData();
	}
}

ZMSG_HANDLER( Zmoil_StereoMode ) {
	CMOIL::Main::StereoMode = model.modeNameToValue.getI( zmsgS( mode ) );
	CMOIL::Main::ColorAnaglyph = zmsgI( color ) == 1;
}

ZMSG_HANDLER( Zmoil_DisplayQuality ) {
	char *mode = zmsgS( mode );
	int   val  = zmsgI( val );
	int   modeVal = model.modeNameToValue.getI( mode );

	if( modeVal == CMOIL::STICK ) {
		CMOIL::Main::DrawBondAsLine = val==0;		
	}
	else if ( modeVal == CMOIL::BACKBONE || modeVal == CMOIL::RNABACKBONE ) {
	    CMOIL::Main::DrawBackboneAsCurve = val==0;
	}
}

ZMSG_HANDLER( Zmoil_DisplaySize ) {
	int increase = zmsgI( increase );
	char *element= zmsgS( elem );
	if( !strcmp( element, "ribbon" ) ) {
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::Main::AAptrs[ i ]->InParm.ribbonWidthFactor *= increase ? 1.1f : .9f;
		}
	}
	else if( !strcmp( element, "stick" ) ) {
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::Main::AAptrs[ i ]->InParm.stickWidthFactor += increase ? 1 : -1;
			if( CMOIL::Main::AAptrs[ i ]->InParm.stickWidthFactor < 1.0 ) {
				CMOIL::Main::AAptrs[ i ]->InParm.stickWidthFactor = 1.0;
			}
		}
	}
}


//-----------------------------------------------------------------------------
// Display Controls Message Handlers

void rotate( float xAxis, float yAxis, float zAxis ) {
	CMOIL::AtomArray *aaptr = CMOIL::Main::AAptrs[ NO_OVERLAP?CMOIL::Main::ImageIndex:CMOIL::Main::Overlap1 ];
	float aa[3];
	aa[0]=-xAxis;
	aa[1]=-yAxis;
	aa[2]=-zAxis;
	CMOIL::axis_to_quat(aa, aaptr->InParm.rotateDegree, CMOIL::Main::LastQuat.q );
	if( manualAlign ) {
		CMOIL::add_quats( CMOIL::Main::LastQuat.q, manualAlign->viewQuat.q, manualAlign->viewQuat.q );
	}
	else {
		CMOIL::add_quats( CMOIL::Main::LastQuat.q, CMOIL::Main::CurQuat.q, CMOIL::Main::CurQuat.q );
	}
}

void translate( FVec3 translate ) {
	if( manualAlign ) {
		manualAlign->viewTranslate.add( translate );
	}
	else {
		CMOIL::Main::Translate.add( translate );
	}
}

ZMSG_HANDLER( Zmoil_DisplayControl ) {

	if( CMOIL::Main::NumImages < 1 || g_AltKey ) {
		return;
	}

	char *cmd = zmsgS( cmd );
	int val   = zmsgI( val );

	int shift = g_ShiftKey ? -1 : 1;
		// rotate opposite direction if shift is being held

	if( !strcmp( cmd, "rotX" ) ) {
		rotate( (float)val*shift, 0.f, 0.f );
	}
	else if( !strcmp( cmd, "rotY" ) ) {
		rotate( 0.f, (float)val*shift, 0.f );
	}
	else if( !strcmp( cmd, "rotZ" ) ) {
		rotate( 0.f, 0.f, (float)val*shift );
	}
	if( !strcmp( cmd, "tranX" ) ) {
		double scaleBias = CMOIL::Main::AAptrs[0]->maxDisplacementXYZ * 0.1f;
		translate( FVec3( (float)scaleBias*val*shift, 0.f, 0.f ) );
		//CMOIL::Main::Translate.x +=   scaleBias * val * shift;
	}
	else if( !strcmp( cmd, "tranY" ) ) {
		double scaleBias = CMOIL::Main::AAptrs[0]->maxDisplacementXYZ * 0.1f;
		translate( FVec3( 0.f, (float)scaleBias*val*shift, 0.f ) );
		//CMOIL::Main::Translate.y +=   scaleBias * val * shift;
	}
	else if( !strcmp( cmd, "tranZ" ) ) {
		double scaleBias = CMOIL::Main::AAptrs[0]->maxDisplacementXYZ * 0.1f;
		translate( FVec3( 0.f, 0.f, (float)scaleBias*val*shift ) );
		//CMOIL::Main::Translate.z +=   scaleBias * val * shift;
	}
	else if( !strcmp( cmd, "scale" ) ) {
		CMOIL::Main::ScaleFactor += val * .05f;
		if( CMOIL::Main::ScaleFactor < 0.f ) {
			CMOIL::Main::ScaleFactor = 0.f;
		}
	}
	else if( !strcmp( cmd, "reset" ) ) {
		CMOIL::Main::ScaleFactor = 1.f;
		Zmoil_eyeZ = .9f;
		for( int i=0; i<3; i++ ) {
			float aa[3] = { 0.f, 0.f, 0.f };
			aa[i] = 1.f;
			CMOIL::axis_to_quat( aa, 0, CMOIL::Main::LastQuat.q );
			CMOIL::axis_to_quat( aa, 0, CMOIL::Main::CurQuat.q );
		}
		CMOIL::Main::Translate.origin();
	}
} 


//-----------------------------------------------------------------------------
// Surface Message Handlers

ZMSG_HANDLER( Zmoil_Surface ) {

	char *cmd  = zmsgS( cmd );
	int val    = zmsgI( val );
	float fval = zmsgF( val );

	// The following was taken/modified from CMOIL::Main::keyboard_zcmd:
	CMOIL::AtomArray *aaptr;
	int ImageIndex = CMOIL::Main::ImageIndex;
	int Overlap1   = CMOIL::Main::Overlap1;
	int NumImages  = CMOIL::Main::NumImages;
	/*
	for (int i=(NO_OVERLAP)?ImageIndex:Overlap1; CMOIL::Main::AAptrs[i]; ) {
		// read at most two structures
	*/
	for (int i=0; i<NumImages; i++ ) {
		aaptr = CMOIL::Main::AAptrs[i];
		if( !strcmp( cmd, "Display" ) ) {
			if( val ) {
		        aaptr->InParm.dispmode |= CMOIL::SURF;
				aaptr->InParm.ribbonpolymode = GL_FILL;
					// ? - I took this from cmoil, but doesn't make sense; probably
					// this mode value is used for more than just the ribbon.
			}
			else {
		        aaptr->InParm.dispmode &= ~CMOIL::SURF;
			}
            aaptr->ReProcessData();
		}
		else if( !strcmp( cmd, "UniqueColors" ) ) {
			aaptr->InParm.uniqueSurfaceColors = val != 0;
		}
		else if( !strcmp( cmd, "Transparent" ) ) {
			aaptr->InParm.transparent = 255 - val * 55;
				// val is 0 or 1 (1 means transparent)
		}
		else if( !strcmp( cmd, "Mesh" ) ) {
			aaptr->InParm.showSurfInMesh = val != 0;  
		}
		else if( !strcmp( cmd, "ShowSelected" ) ) {
			aaptr->InParm.showSubSurf = !(aaptr->InParm.showSubSurf);  // highlight the area
		}
		else if( !strcmp( cmd, "CavityOnly" ) ) {
			aaptr->InParm.showCavity = val != 0;  
			if ( val ) {
				aaptr->printCavity();
				aaptr->InParm.cutoffval = 1.e+20f;
			}
		}
		else if( !strcmp( cmd, "RaviNormal" ) ) {
			val ? aaptr->InParm.dispmode |= CMOIL::RAVINORM : aaptr->InParm.dispmode &= ~CMOIL::RAVINORM;
		}
		else if( !strcmp( cmd, "RaviNormalCoarse" ) ) {
			val ? aaptr->InParm.dispmode |= CMOIL::RAVINORMCOARSE : aaptr->InParm.dispmode &= ~CMOIL::RAVINORMCOARSE;
		}
		else if( !strcmp( cmd, "ColorNone" ) ) {
			aaptr->InParm.AAcolorGroups &= ~(CMOIL::SURF_COLOR_ATOM|CMOIL::SURF_COLOR_RES);
		}
		else if( !strcmp( cmd, "ColorAtoms" ) ) {
			aaptr->InParm.AAcolorGroups &= ~CMOIL::SURF_COLOR_RES;
			aaptr->InParm.AAcolorGroups |= CMOIL::SURF_COLOR_ATOM;
		}
		else if( !strcmp( cmd, "ColorResidues" ) ) {
			aaptr->InParm.AAcolorGroups &= ~CMOIL::SURF_COLOR_ATOM;
			aaptr->InParm.AAcolorGroups |= CMOIL::SURF_COLOR_RES;
		}
		else if( !strcmp( cmd, "CutoffClear" ) ) {
			aaptr->InParm.cutoffval = 1.e+20f; 
		}
		else if( !strcmp( cmd, "CutoffPlane" ) ) {
			aaptr->setSurfCutoffPlane(0);
		}
		else if( !strcmp( cmd, "CutoffCycleAxis" ) ) {
			aaptr->setSurfCutoffPlane(1);
		}
		else if( !strcmp( cmd, "CutoffFlip" ) ) {
			aaptr->flipSurfCutoffPlane();
		}
		else if( !strcmp( cmd, "CutoffFarther" ) ) {
			aaptr->InParm.cutoffval -= aaptr->InParm.cutoffstep; 
		}
		else if( !strcmp( cmd, "CutoffNearer" ) ) {
			aaptr->InParm.cutoffval += aaptr->InParm.cutoffstep; 
		}
		else if( !strcmp( cmd, "CutoffJumpDistance" ) ) {
			aaptr->InParm.cutoffstep = fval;
		}
		else {
			printf( "Zmoil_Surfance handler received an unknown command.\n" );
			assert( false && "bad cmd to Zmoil_Surface" );
			return;
		}
/*
		if ( NO_OVERLAP || i == CMOIL::Main::Overlap2 ) {
			break;
		}
		else {
			i = CMOIL::Main::Overlap2;
		}
		*/
    }
}

ZMSG_HANDLER( Zmoil_CutoffJumpDistance ) {
	// We recieve a float value between .5 and 64.  Want to round this
	// to nearest power of 2.

	float vals[] = { .5, 1, 2, 4, 8, 16, 32, 64 };
	const int count = sizeof( vals ) / sizeof( vals[0] );
	float val = zmsgF( val );

	for( int i=0; i<count; i++ ) {
		if( val <= vals[i] ) {
			val = vals[i];
			break;
		}
	}

	// @TODO: there is a way to have the message handler pass us the pointer to ZUI I think.
	ZUI::zuiFindByName( "cutoffJumpDist" )->putF( "val", val );
	zMsgQueue( "type=Zmoil_Surface cmd=CutoffJumpDistance val=%f", val );
}

//-----------------------------------------------------------------------------
// Structure Message Handlers

void updateIndexSlider( bool labelOnly ) {
	// UPDATE the ZUISlider to reflect the state of the model
	static ZUI *slider = 0;
	static ZUI *label = 0;
	if( !slider ) {
		slider = ZUI::zuiFindByName( "structSlider" );
		label  = ZUI::zuiFindByName( "structIndex" );
	}
	if( slider ) {
		CMOIL::AtomArray *aa = modelListGetSelected( "dataList" );
		if( !aa ) {
			aa = CMOIL::Main::AAptrs[ 0 ];
		}
		if( aa ) {
			if( !labelOnly ) {
				float val = (float)aa->currentmol / aa->InParm.structno;
				slider->putF( "val", val );
				slider->putF( "lastVal", val );
					// prevent generating message from ZUISlider
			}
			label->putS( "text", ZTmpStr("Structure Index:  %d / %d", aa->currentmol, ( aa->InParm.structno > 0 ? aa->InParm.structno : 1 ) ) );
		}
		slider->dirty();
		label->dirty();
	}
}

ZMSG_HANDLER( Zmoil_StructSlider ) {
	CMOIL::AtomArray *aa = modelListGetSelected( "dataList" );
	if( !aa ) {
		aa = CMOIL::Main::AAptrs[ 0 ];
	}
	if( aa ) {
		int frame = (int)min( zmsgF(val) * aa->InParm.structno + 1, aa->InParm.structno );
		int step = frame - aa->currentmol;
		if( step ) {
			zMsgQueue( "type=Zmoil_StructNext step=%d fromSlider=1", step );
		}
	}
}

ZMSG_HANDLER( Zmoil_StructNext ) {
	CMOIL::AtomArray *aa = modelListGetSelected( "dataList" );
	if( !aa ) {
		aa = CMOIL::Main::AAptrs[ 0 ];
	}

	int step  = zmsgI( step );

	if ( aa->InParm.filetag != CMOIL::fPDB ) {
		// use random access version: this only works for some types of
		// data formats.
		int currentmol = aa->currentmol - 1 + step;
			// note: aaptr->currentmol is post-incremented after each read thus appears
			// as 1-based, so we convert to 0 based by -1.
		if ( currentmol < 0 ) {
			currentmol = aa->InParm.structno - 1;
		}
		if ( currentmol >= aa->InParm.structno ) {
			currentmol = 0;
		}
		aa->ReadStruct( currentmol );
	} 

	CMOIL::Main::PickQ.print_pick_info();

	// UPDATE the slider & index label
	updateIndexSlider( zmsgI(fromSlider)!=0 );
}
/*
ZMSG_HANDLER( Zmoil_StructNext ) {
	CMOIL::AtomArray *aaptr;
	int ImageIndex = CMOIL::Main::ImageIndex;
	int Overlap1   = CMOIL::Main::Overlap1;
	int Overlap2   = CMOIL::Main::Overlap2;
	int NumImages  = CMOIL::Main::NumImages;

	int step  = zmsgI( step );

	// CMOIL::Main::MsgBrd.erase();

	int i;
	for (i=(NO_OVERLAP)?ImageIndex:Overlap1; CMOIL::Main::AAptrs[i]; ) {
		aaptr = CMOIL::Main::AAptrs[i];
		if ( aaptr->InParm.filetag != CMOIL::fPDB ) {
			// use random access version: this only works for some types of
			// data formats.
			int currentmol = aaptr->currentmol - 1 + step;
				// note: aaptr->currentmol is post-incremented after each read thus appears
				// as 1-based, so we convert to 0 based by -1.
			if ( currentmol < 0 ) {
				currentmol = aaptr->InParm.structno - 1;
			}
			if ( currentmol >= aaptr->InParm.structno ) {
				currentmol = 0;
			}
			aaptr->ReadStruct( currentmol );
		} 
		if ( NO_OVERLAP || i==Overlap2 ) {
			break;
		}
		else {
			i = Overlap2;
		}
	}

	// cause pick information to be displayed for new structure index
	// NOPICKMODE
	if( true || CMOIL::Main::PickMode ) {
		CMOIL::Main::PickQ.print_pick_info();
	}
	// UPDATE the slider & index label
	updateIndexSlider( zmsgI(fromSlider)!=0 );
}
*/
ZMSG_HANDLER( Zmoil_ModelDraw ) {
	CMOIL::AtomArray *aa = (CMOIL::AtomArray*)STRING_TO_PTR( msg, model );
	if( aa ) {
		aa->bDraw = zmsgI( draw ) == 1;
	}
}

ZMSG_HANDLER( Zmoil_ModelBlink ) {
	CMOIL::AtomArray *aa = (CMOIL::AtomArray*)STRING_TO_PTR( msg, model );
	if( aa ) {
		aa->bBlink = zmsgI( blink ) == 1;
	}
}

ZMSG_HANDLER( Zmoil_ModelCycle ) {
	CMOIL::AtomArray *aa = (CMOIL::AtomArray*)STRING_TO_PTR( msg, model );
	if( aa && aa->bCycle != (zmsgI( cycle )==1) ) {
		aa->bCycle = zmsgI( cycle ) == 1;
		if( !zmsgI( cycle ) ) {
			// revert to display of the models that were previously selected...
			for( int i=0; i<aa->numModels; i++ ) {
				ZUI *z = ZUI::zuiFindByName( ZTmpStr( "pdbmodel%d", i ) );
				if( z ) {
					if( z->getI( "selected" ) ) {
						aa->models[ i ].display = 1;
					}
					else {
						aa->models[ i ].display = 0;
						trace( "  ********** display  for model %d was turned off!\n", i );
					}
				}
			}
			aa->setSkip();
			aa->setModelChainDisplay();
		}
	}
}

ZMSG_HANDLER( Zmoil_ModelShow ) {
	CMOIL::AtomArray *aa = (CMOIL::AtomArray*)STRING_TO_PTR( msg, model );
	if( aa ) {
		int show  = zmsgI( show );
		int modelIdx;
		if( zmsgHas( modelIdx ) ) {
			modelIdx = zmsgI( modelIdx );
		}
		else {
			int modelName = zmsgI( modelName );
				// this 'modelName' is the numeric identifier used in the actual PDB
			for( modelIdx=0; modelIdx<aa->numModels; modelIdx++ ) {
				if( aa->models[ modelIdx ].modelName == modelName ) {
					break;
				}
			}
		}
		if( modelIdx < aa->numModels && aa->models[ modelIdx ].display != show ) {
			aa->DisplayPdbModel( modelIdx, show == 1 );
		}
		// toggling chain status automatically updates the centroid
		// and applies the translation to center the view: if we don't
		// want this, revert it:
		if( ZUI::zuiFindByName( "comNever" )->getI( "selected" ) ) {
			aa->resetCrd();
		}
	}
}

ZMSG_HANDLER( Zmoil_ChainShow ) {
	CMOIL::AtomArray *aa = (CMOIL::AtomArray*)STRING_TO_PTR( msg, model );
	if( aa ) {
		int show  = zmsgI( show );
		int chain;
		if( zmsgHas( chain ) ) {
			chain = zmsgI( chain );
		}
		else {
			int chainId = zmsgI( chainId );
			for( chain=0; chain<aa->numChains; chain++ ) {
				if( (int) aa->chains[ chain ].id == chainId ) {
					break;
				}
			}
		}
		
		if( chain < aa->numChains && aa->chains[ chain ].display != show ) {
			aa->DisplayPdbChain( chain, show == 1 );
		}
		// toggling chain status automatically updates the centroid
		// and applies the translation to center the view: if we don't
		// want this, revert it:
		if( ZUI::zuiFindByName( "comNever" )->getI( "selected" ) ) {
			aa->resetCrd();
		}
	}
}

ZMSG_HANDLER( Zmoil_UpdateCoM ) {
	CMOIL::AtomArray *aaptr;

	// Get picked coordinate first in case we want to center on it;
	// if we have multiple structs, they'll all want to translate
	// by this amount to preserve any alignment, and this may 
	// alter the pickedAtom coordinate!
	DVec3 pickedAtom;
	CMOIL::Coordinate pickedCrds[3];
	int nPicked = CMOIL::Main::PickQ.getPickedCrds(pickedCrds);
	if( nPicked ) {
		pickedAtom.x = pickedCrds[0].x;
		pickedAtom.y = pickedCrds[0].y;
		pickedAtom.z = pickedCrds[0].z;
	}

	// Two passes, such that if a model is aligned to a target,
	// we ensure updating the target centroid first, which is 
	// needed by the second.
	int i,pass;
	for( pass=0; pass<2; pass++ ) {
		for (i=0; i<CMOIL::Main::NumImages; i++ ) {
			aaptr = CMOIL::Main::AAptrs[i];

			if( (!aaptr->alignToModel && pass==0) || (aaptr->alignToModel && pass==1) ) {

				if( zmsgI( clear ) ) {
					aaptr->resetCrd();
				}
				else if( zmsgIs( when, never ) ) {
					aaptr->updateCoM = CMOIL::AtomArray::CoM_Never;
				}
				else if ( zmsgIs( when, first ) ) {
					aaptr->updateCoM = CMOIL::AtomArray::CoM_First;
				}
				else if ( zmsgIs( when, always ) ) {
					aaptr->updateCoM = CMOIL::AtomArray::CoM_Always;
				}
				else if ( zmsgIs( when, now ) ) {
					aaptr->resetCrd();
					aaptr->computeCentroid();
					if( aaptr->alignToModel ) {
						aaptr->applyTransform( aaptr->alignMatrix );
						aaptr->applyTranslation( aaptr->alignToModel->centroid, true );
					}
					else {
						aaptr->applyTranslation( aaptr->centroid, true );
					}
				}
				else if ( zmsgIs( which, lastpickedatom ) ) {
					aaptr->applyTranslation( pickedAtom, true );
				}
			}
		}
	}
}

//-----------------------------------------------------------------------------
// Utility (program generated) Message Handlers

ZMSG_HANDLER( Zmoil_KillZUI ) {
	ZUI *z = ZUI::zuiFindByName( zmsgS( zui ) );
	if( z ) {
		z->detachFrom();
		z->killChildren();
		z->die();
		zMsgQueue( "type=ZUILayout" );
	}
}


// General-purpose message handlers:



//-----------------------------------------------------------------------------
// plugin

ZPLUGIN_BEGIN( zmoil );

void startup() {
	trace( "Running %s\n", getVersionString() );
	//emailMsgTo( "Starting zmoil!",  "blomcode@gmail.com" );
	
	// BUILD handy lookup tables for mapping ZLAB info to CMOIL info
	setupCmoil();
	
	// LOAD the control panel, get ptr to main plugin window
	char buf[255];
	getcwd( buf, 255 );
	char *ppath = pluginPath( "_zmoil.zui" );
	ZUI::zuiExecuteFile( pluginPath( "_zmoil.zui" ) );
	zMsgQueue( "type=ZUILayout" );
	g_pluginPanel = ZUI::zuiFindByName( "pluginPanel" );
	g_uiPanel     = ZUI::zuiFindByName( "cmoilZUI" );
	
	// RESIZE the plugin panel (main drawing area) based on the UI size:
	ZTmpStr uiWidth( "%d", UI_PANEL_WIDTH );
	zMsgQueue( "type=ZUISet key=layoutManual_x val='%s' toZUI=pluginPanel", uiWidth.s );
	zMsgQueue( "type=ZUISet key=layoutManual_w val='W %s -' toZUI=pluginPanel; type=ZUILayout toZUI=pluginPanel", uiWidth.s );
	
	
	// LOAD fonts
	// zglFontLoad( "verdana11", "verdana.ttf", 11 );
	zglFontLoad( "courier10",  zlabCorePath( "cour.ttf" ),  10 );
	
	
	// LOAD user defaults
	loadPathOptions();
	if( options.getI("Zmoil_skipTutorial") ) {
		zMsgQueue( "type=ZUIHide toZUI=zmoilTutorial" );
	}
	
	// Initialize openGL matrices, Lighting
	CMOIL::trackball( CMOIL::Main::CurQuat.q, 0.0, 0.0, 0.0, 0.0);
	CMOIL::Main::InitLightSource();
	glEnable(GL_LIGHT0);
		// can we enable a given light and leave it enabled, and control whether it
		// actually gets used just be enable/disabling lighting?  (with no perf hit?)
	glEnable( GL_NORMALIZE );

	// On Darwin(osx), the open command is used to launch zmoil, and cmdline
	// args aren't passed. In this case, the file to open is written into
	// our app bundle as zmoil.arg:
#ifdef __APPLE__
	char *argsFile = "./zlab.app/zmoil.args";
	if( zWildcardFileExists( argsFile ) ) {
		// Get the name of the input file, set this as nonOptionArg2 (used below),
		// and delete the file (it is written freshly my MOIL each time,
		// and we don't want it read when starting standalone)
		ZStr *fileLines = zStrReadFile( argsFile );
		if( fileLines ) {
			options.putS( "nonOptionArg2", fileLines->getS( 0 ) );
			strcpy( EXEFolder, fileLines->getS( 1 ) ); 
			trace( "zmoil.args was successfully read.\n" );
			trace( "   fileToLoad = %s\n", options.getS( "nonOptionArg2" ) );
			char buf[256];
			getcwd( buf, 256 );
			trace( "   working from folder %s\n", buf );
		}
		zStrDelete( fileLines );
		//unlink( argsFile );
	}
#endif
	// If we were launched from the TCL gui for MOIL, the input file
	// to open will be found at argv[3]; otherwise look for the last file
	// that was loaded. In either case, load with a delay to allow UI init
	// to have completed.
	char *fileToLoad = options.getS( "nonOptionArg2", 0 );
	if( !fileToLoad ) {
    	fileToLoad = options.getS( "Zmoil_lastLoaded", 0 );
    }
	if( fileToLoad /*&& zWildcardFileExists( fileToLoad )*/ ) {
		trace( "At startup, fileToLoad is %s\n", fileToLoad );
		zMsgQueue( "type=Zmoil_Open filespec='%s' __countdown__=4", fileToLoad );
	}
	strcpy( EXEFolder, options.getS( "nonOptionArg3", "" ) );
	if( strlen( EXEFolder ) ) {
		strcat( EXEFolder, "/" );
	}
	CMOIL::Main::ExeDir = EXEFolder;
	if( !EXEFolder[1] ) {
		strcpy( EXEFolder, exeName );
		int slashCount=0;
		char *p = EXEFolder + strlen( EXEFolder ) - 1;
		while( p > EXEFolder ) {
			if( *p == '\\' || *p == '/' ) {
				if( ++slashCount == 2 ) {
					p++;
					*p=0;
					break;
				}
			}
			p--;
		}
	}
	
	trace( "EXEFolder has been set to '%s'\n", EXEFolder );
	trace( "gTimeout is %.4lf\n", gTimeout );
}

void shutdown() {
	// SAVE user defaults
	savePathOptions();
}
ZMSG_HANDLER( Zmoil_Shutdown ) {
	shutdown();
	exit( 0 );
}

void update() {
	
	operationPending( -1 );
		// decrement the frame counter for operations pending

	useDirtyRects = Zmoil_UIdirtyRects;
	drawDirtyRects = Zmoil_UIdirtyRectsDebug;
	// this is toggling vars in zui.cpp
	
#define SLEEP_TIMEOUT_SECS 3.0
#define SLEEP_LENGTH_MSECS 200
	ZTime now = zTimeNow();
	static int mouseLastX, mouseLastY, keyLast;
	if( zMouseMsgX!=mouseLastX || zMouseMsgY!=mouseLastY || zuiKeyEvent!=keyLast ) {
		mouseLastX = zMouseMsgX;
		mouseLastY = zMouseMsgY;
		keyLast    = zuiKeyEvent;
		gLastActivityTime = now;
	}
	if( now - gLastActivityTime > SLEEP_TIMEOUT_SECS  ) {
		if( !(CMOIL::Main::MovieMode || CMOIL::Main::Spinning) ) {
			zTimeSleepMils( SLEEP_LENGTH_MSECS );
		}
	}
	// Update realtime UI state
	if( CMOIL::Main::MovieMode ) {
		updateIndexSlider();
	}
	
	// check status of debug vars
	static int peptidePlaneAtoms=0;
	if( Zmoil_PeptidePlaneAtoms != peptidePlaneAtoms ) {
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::Main::AAptrs[i]->setSkip();
		}
		peptidePlaneAtoms = Zmoil_PeptidePlaneAtoms;
	}
	static int HBondsWater=0;
	if( Zmoil_DisplayHBondsWater != HBondsWater ) {
		for( int i=0; i<CMOIL::Main::NumImages; i++ ) {
			CMOIL::Main::AAptrs[i]->properties.del( "hBondsCalculated" );
		}
		HBondsWater = Zmoil_DisplayHBondsWater; 
	}
}

void handleMsg( ZMsg *msg ) {
	gLastActivityTime = zTimeNow();
	
	if( zmsgIs( type, KeyDown ) ) {
		if( zmsgIs( which, f8 ) )  {
			// toggle visibility of the main UI
			zMsgQueue( "type=Zmoil_ToggleUI" );
			zMsgUsed();
		}
		// @TODO: wanted to make ctrl-versions of these key bindings to rotate the opposite way.
		// But appears that the normal keybinding still gets triggered in addition to this 
		// handler when ctrl-key is pressed.  (This seems different than alt_z for fps?)
		// Another options is to look at the async keystate in the handler to see if the shift
		// or ctrl key is pressed, and rotate opposite direction if so.
		/*
		 else if( zmsgIs( which, ctrl_x) ) {
		 zMsgQueue( "type=Zmoil_DisplayControl cmd=rotX val=1" );
		 }
		 else if( zmsgIs( which, ctrl_y) ) {
		 zMsgQueue( "type=Zmoil_DisplayControl cmd=rotY val=1" );
		 }
		 else if( zmsgIs( which, ctrl_z) ) {
		 zMsgQueue( "type=Zmoil_DisplayControl cmd=rotZ val=1" );
		 }
		 */
		
		else if( zmsgIs( which, lshift ) || zmsgIs( which, rshift )) {
			g_ShiftKey = 1;
		}
		if( zmsgIs( which, lalt ) || zmsgIs( which, ralt )) {
			g_AltKey = 1;
		}
		
		
		else if( g_uiPanel->getI( "hidden" ) ) {
			// do keys we want to support while UI is hidden:
			if( !zmsgIs( which, lshift ) ) {
				int debug=1;
			}
			if( zmsgIs( which, s ) ) {
				zMsgQueue( "type=Zmoil_SaveImage" );
			}
			else if( zmsgIs( which, m ) ) {
				zMsgQueue( "type=Zmoil_ProgramMode mode=movie val=%d", !CMOIL::Main::MovieMode );
			} 
			else if( zmsgIs( which, p ) ) {
				zMsgQueue( "type=Zmoil_ProgramMode mode=pick val=%d", !CMOIL::Main::PickMode );
			} 
			else if( zmsgIs( which, l ) ) {
				zMsgQueue( "type=Zmoil_UpdateCoM which=lastpickedatom" );
			} 
			else if( zmsgIs( which, a ) ) {
				if( 1 || CMOIL::Main::PickMode ) {
					zMsgQueue( "type=Zmoil_PickInfo info=angle" );
				}
			} 
			else if( zmsgIs( which, d ) ) {
				if( 1 || CMOIL::Main::PickMode ) {
					zMsgQueue( "type=Zmoil_PickInfo info=distance" );
				}
			} 
			else if( zmsgIs( which, t ) ) {
				if( 1 || CMOIL::Main::PickMode ) {
					zMsgQueue( "type=Zmoil_PickInfo info=torsion" );
				}
			} 
			else if( zmsgIs( which, r ) ) {
				zMsgQueue( "type=Zmoil_ProgramMode mode=rock val=%d", !CMOIL::Main::Rocking );
			} 
			else if( zmsgIs( which, n ) ) {
				zMsgQueue( "type=Zmoil_StructNext step=1" );
			} 
			else if( zmsgIs( which, b ) ) {
				zMsgQueue( "type=Zmoil_StructNext step=-1" );
			} 
			else if( zmsgIs( which, c ) ) {
				// zMsgQueue( "type=Zmoil_ChainNext" );
			} 
			else if( zmsgIs( which, x ) ) {
				zMsgQueue( "type=Zmoil_DisplayControl cmd=rotX val=-1" );
			}
			else if( zmsgIs( which, y ) ) {
				zMsgQueue( "type=Zmoil_DisplayControl cmd=rotY val=-1" );
			}
			else if( zmsgIs( which, z ) ) {
				zMsgQueue( "type=Zmoil_DisplayControl cmd=rotZ val=-1" );
			}
			else if( zmsgIs( which, = ) ) {
				zMsgQueue( "type=Zmoil_DisplayControl cmd=scale val=1" );
			}
			else if( zmsgIs( which, - ) ) {
				zMsgQueue( "type=Zmoil_DisplayControl cmd=scale val=-1" );
			}
			else if( zmsgIs( which, 1 ) ) {
				zMsgQueue( "type=Zmoil_DisplayControl cmd=reset" );
			}
		}
	}
	else if( zmsgIs( type, KeyUp ) ) {
		if( zmsgIs( which, lshift ) || zmsgIs( which, rshift )) {
			g_ShiftKey = 0;
		}
		if( zmsgIs( which, lalt ) || zmsgIs( which, ralt )) {
			g_AltKey = 0;
		}
	}
}


ZPLUGIN_EXPORT_PROPERTY( shadowGardenInterface, "1" );
ZPLUGIN_EXPORT_SYMBOL( startup );
ZPLUGIN_EXPORT_SYMBOL( shutdown );
// ZPLUGIN_EXPORT_SYMBOL( render );
	// not exporting this causes an early out in zuiPluginView::render() that saves us a lot
	// of opengl state management, and helps our framerate.

ZPLUGIN_EXPORT_SYMBOL( update );
ZPLUGIN_EXPORT_SYMBOL( handleMsg );

ZPLUGIN_END;

//-----------------------------------------------------------------------------
// Some utility fns that might get factored out of this plugin later...

char * replaceExtension( char *filepath, char *newExt ) {
	// SET or REPLACE extension in filepath. This will replace multitple
	// extenstions in cases like myfile.ext.ext
	// RETURNS pointer to this static: (be careful!)
	static char replaceExtensionResults[512];
	assert( newExt && strlen( filepath ) + strlen( newExt ) < 512 );

	ZFileSpec fs( filepath );
	char *ext = fs.getExt();
	if( !*ext ) {
		sprintf( replaceExtensionResults, "%s.%s", fs.get(), newExt );
	}
	else if( !strcmp( ext, newExt ) ) {
		strcpy( replaceExtensionResults, filepath );
	}
	else {
		char extWithPeriodRegEx[16];
		char newExtWithPeriod[16];
		sprintf( extWithPeriodRegEx, "\\.%s", ext );
		sprintf( newExtWithPeriod, ".%s", newExt );
		ZStr fname( fs.get() );
		zStrReplace( &fname, extWithPeriodRegEx, newExtWithPeriod );	
		strcpy( replaceExtensionResults, fname.getS() );
	}
	return replaceExtensionResults;
}

void glPrint( int x, int y, char *text,  char *font ) {
	glPushAttrib( GL_ALL_ATTRIB_BITS );
	glDisable( GL_LIGHTING );
	glDisable( GL_COLOR_MATERIAL );
	glDisable( GL_CULL_FACE );
	glDisable( GL_DEPTH_TEST );
	glDisable( GL_BLEND );
	glDisable( GL_TEXTURE_2D );
	glDisable( GL_FOG );
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL );

	glMatrixMode( GL_PROJECTION );
	glPushMatrix();
	glLoadIdentity();
	glMatrixMode( GL_MODELVIEW );
	glPushMatrix();
	zglPixelMatrixFirstQuadrant();

	char *_font = font ? font : (char*)"controls";
	zglFontPrint( text, (float)x, (float)y, _font );

	glMatrixMode( GL_PROJECTION );
	glPopMatrix();
	glMatrixMode( GL_MODELVIEW );
	glPopMatrix();

	glPopAttrib();
}

void tempMessage( char *msg, char color, ZTime secondsToDisplay ) {
	switch( color ) {
		case 'r':
			tmpMessageColor[0] = 255;
			tmpMessageColor[1] = 0;
			tmpMessageColor[2] = 0;
		break;

		default:
			tmpMessageColor[0] = tmpMessageColor[1] = tmpMessageColor[2] = 255;
	} 
	strcpy( tmpMessage, msg );
	tmpMessageDisplayUntil = zTimeNow() + secondsToDisplay;
}

int imageWidth;
int imageHeight;
#define BYTES_PER_PIXEL 4
unsigned char * capturePluginImage( char *filename, bool notify ) {
	assert( g_pluginPanel );

	imageWidth = (int)g_pluginPanel->w & ~1;
		// force even - required if image will be converted to mpeg movie
	imageHeight = (int) g_pluginPanel->h;

	// ALLOC pixel buffer to fill with screen data
	static unsigned char * pixels = 0;
	static long pixelsize = 0;
	long sizeRequest = imageWidth * imageHeight * BYTES_PER_PIXEL;
	if ( sizeRequest > pixelsize ) {
		if( pixels ) {
			delete [] pixels;
		}
		pixels = new unsigned char[ sizeRequest ];
		pixelsize = sizeRequest;
	}

	// READ screen to pixel buffer
	glReadPixels( (int)g_pluginPanel->x, (int)g_pluginPanel->y, imageWidth, imageHeight, GL_RGBA,  GL_UNSIGNED_BYTE, pixels);

	// OPTIONALLY WRITE TO FILE
	if( filename ) {
		// CONVERT to ZBitmapDesc format & invert lines.
		ZBitmapDesc bmd;
		bmd.init( imageWidth, imageHeight, BYTES_PER_PIXEL, (char *)pixels );
		zBitmapDescInvertLines( bmd );

		#ifdef __BIG_ENDIAN__
		zBitmapDescInvertBGR( bmd );
			// Flips the byte order of the bits in a bitmap of depth 3 or 4
		#endif

		// WRITE pixels to file
		int result = zGraphFileSave( filename, &bmd );
		if( notify ) {
			if( result ) {
				tempMessage( ZTmpStr( "Saved %s", filename ) );
			}
			else {
				tempMessage( ZTmpStr( "Error: Can't save %s", filename ), 'r' );
			}
		}
	}

	// RETURN raw pixels
	return pixels;
}

void updateImageCapture() {

	if( operationPending() ) {
		// some multi-frame opertaion (such as loading an mfm and preparing
		// for display) have not completed, so wait to capture image.
		return;
	}

	// called by render for multiframe captures
	static FIMULTIBITMAP *multi = NULL;
		// for animated GIF capture

	char *capFile = model.viewInfo.getS( "captureFile" );
	assert( capFile );
	ZFileSpec capFS( capFile );
	char *ext = capFS.getExt();
	assert( ( !strcmp( ext, "png" ) || !strcmp( ext, "tif" ) || !strcmp( ext, "mpeg" ) || !strcmp( ext, "gif" ) ) );

	CMOIL::AtomArray *aa = CMOIL::Main::AAptrs[ CMOIL::Main::ImageIndex ];

	//ZTime time = zTimeNow();
	//trace( "Writing image frame %d at time %g ... ", aa->currentmol, time );

	if( !strcmp( ext, "mpeg" ) || !strcmp( ext, "png" ) ) {
		// CONSTRUCT filename for current frame
		char filename[256];
		strcpy( filename, model.viewInfo.getS( "captureMovieBasename" ) );
		if( setListIndex >= 0 ) {
			// this is a quick hack to capture an image for each file in a set, but
			// it precludes the possibility of capturoing multiple structs from a single
			// dcd file, which could have been loaded via the set mechanism
			char *setFile = model.viewInfo.getS( "filename", "" );
			ZFileSpec fs( setFile );
			strcat( filename, ZTmpStr( "%04d_%s.png", setListIndex+1, fs.getFile() ) );

		}
		else {
			strcat( filename, ZTmpStr( "%04d.png", aa->currentmol ) );
		}

		// CAPTURE
		capturePluginImage( filename, false );
	}
	else if( !strcmp( ext, "gif" ) ) {	// writing to a multi-frame gif
		unsigned char *pixels = capturePluginImage( NULL, false );

		// I don't understand why I need to invert to bgr, but it seems I do.
		ZBitmapDesc bmd;
		bmd.init( imageWidth, imageHeight, BYTES_PER_PIXEL, (char *)pixels );
		zBitmapDescInvertBGR( bmd );
	
		//#define FI_RGBA_RED_MASK		0x00FF0000
		//#define FI_RGBA_GREEN_MASK	0x0000FF00
		//#define FI_RGBA_BLUE_MASK		0x000000FF
			// these are defined by freeimage, but messing with them doesn't change 
			// the inverted red and blue, thus the call to InvertBGR above.
		FIBITMAP *dib = FreeImage_ConvertFromRawBits( pixels, imageWidth, imageHeight, BYTES_PER_PIXEL*imageWidth, BYTES_PER_PIXEL*8,		
													   FI_RGBA_RED_MASK, FI_RGBA_GREEN_MASK, FI_RGBA_BLUE_MASK, 1 );

		FIBITMAP *dib24= FreeImage_ConvertTo24Bits( dib );
		FIBITMAP *dib8 = FreeImage_ColorQuantize( dib24, FIQ_WUQUANT );
			// convert from to 24bit then to 8bit palettized (required for GIF
			
		// STORE frame with animation metadata:
		if( !multi ) {
			multi = FreeImage_OpenMultiBitmap( FIF_GIF, capFile, true, false, true );
			assert( multi );
		}
		if( multi ) {
			long frameTime = aa ? (long)aa->InParm.sleepseconds : 200;
			FreeImage_SetMetadata( FIMD_ANIMATION, dib8, NULL, NULL );
			FITAG *tag = FreeImage_CreateTag();
			if( tag ) {
				FreeImage_SetTagKey( tag, "FrameTime" );
				FreeImage_SetTagType( tag, FIDT_LONG );
				FreeImage_SetTagCount( tag, 1 );
				FreeImage_SetTagLength( tag, 4 );
				FreeImage_SetTagValue( tag, &frameTime );
				FreeImage_SetMetadata( FIMD_ANIMATION, dib8, FreeImage_GetTagKey( tag ), tag );
				FreeImage_DeleteTag( tag );
			}
			FreeImage_AppendPage( multi, dib8 );
		}

		FreeImage_Unload( dib );
		FreeImage_Unload( dib24 );
		FreeImage_Unload( dib8 );
	}
	else {
		assert( false && "bad captureFile spec in updateImageCapture()" );
	}

	//ZTime time2 = zTimeNow();
	//trace( "done (took %g seconds)\n", time2 - time );

	// INCREMENT to next struct, clear capture flag if done.
	if( setListIndex >= 0 ) {
		zMsgQueue( "type=Zmoil_SetListNext" );
		if( setListIndex == setList.count - 1 ) {
			model.viewInfo.putI( "captureFrames", 0 );
			tempMessage( ZTmpStr( "Saved %d frames to %s", setList.count, !strcmp( ext, "mpeg" ) ?
							 model.viewInfo.getS( "captureMovieFolder" ) : capFile ) );
		}
	}
	else if( aa ) {
		aa->ReadNextStruct();

		if( aa->currentmol == 1 ) {
			// 1 because currentmol is post-incremented after read
			model.viewInfo.putI( "captureFrames", 0 );
			tempMessage( ZTmpStr( "Saved %d frames to %s", aa->InParm.structno, !strcmp( ext, "mpeg" ) ?
							 model.viewInfo.getS( "captureMovieFolder" ) : capFile ) );

			// CLOSE the multiframe GIF if open:
			if( multi ) {
				FreeImage_CloseMultiBitmap( multi );

				//ZTime time3 = zTimeNow();
				//trace( "CloseMultiBitmap completed after %g seconds.\n", time3 - time2 );

				multi = NULL;
			}
		}
	}
}
//-------------------------------------------------------------------------

int operationFramesLeft=0;
int operationPending( int inc ) {
	operationFramesLeft += inc;
	if( operationFramesLeft < 0 ) {
		operationFramesLeft = 0;
	}
	return operationFramesLeft;
}











