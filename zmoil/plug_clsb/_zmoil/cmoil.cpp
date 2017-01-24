//*************************************************************************************************
//*  Filename:   cmoil.cpp
//*
//*  Description: 
//*    Main loop of graphic display in cmoil program.
//*
//*  Modification History:
//*  
//*  Date           Developer   Description
//*  ------------   ----------- ---------------------------------------------------------------------
//*  Aug. ?? 2000   ??      Initial Development
//*  Oct. 09 2000   Baohua Wang Seperated header file cmoil.hpp
//*  Mar.    2001   Baohua Wang Add text display
//*  Jul. 18 2001   Baohua Wang Reshape to support mutiple structure display and overlap
//*  Jul. 28 2005   B. Wang     Linux support
//*************************************************************************************************

// module 
#include "_zmoil.h"
#include "cmoil.h"
#include "cmoil_globals.h"
#include "cmoil_pick.h"
#include "cmoil_camera.h"
#include "cmoil_bspline.h"

// zbslib
#include "ztime.h"
#include "zhashtable.h"
#include "zui.h"
#include "ztmpstr.h"
#include "zwildcard.h"
#include "zgltools.h"
#include "tbutil.h"


#pragma warning( once : 4244 )
	// warning: conversion double to float

extern int gEnableAxes;
extern float Zmoil_eyeZ;
extern int gEnableStructInfoDisplay;
extern int Zmoil_PeptidePlaneVectors;
extern int Zmoil_PthDisplayStructs;
extern int Zmoil_PthDisplayAlphaLinear;
extern int Zmoil_VelocityDisplay;

extern double Zmoil_raviTranslateX;
extern double Zmoil_raviTranslateY;
extern double Zmoil_raviTranslateZ;

extern void trace( char *fmt, ... );
extern ZTime globalTime;
extern int globalTimerOn;


extern CMOIL::AtomArray *manualAlign;

namespace CMOIL {

tRGBA red   ( 255,0,0  );
tRGBA green ( 0,255,0  );
tRGBA blue  ( 0,0,255  );

tRGBA yellow( 255,255,0 );
tRGBA orange( 255,127,0 );
tRGBA cyan  ( 0,255,255  );
tRGBA purple( 255,0,255 );

tRGBA ltRed ( 255,128,128  );
tRGBA ltGreen ( 128,255,128  );
tRGBA ltBlue ( 128,128,255  );
tRGBA ltYellow( 128,128,0 );


tRGBA white ( 255,255,255  );

tRGBA zmoilColors[12] = { red, orange, yellow, green, blue, purple, cyan, ltRed, ltGreen, ltBlue, ltYellow, white };




/**************************************************************************************************/

// static variables of class State
AtomArray *State::AAptrs[ARRAY_MAXNUM];	// All input Structures
PickedAtoms State::PickQ;	// most recently selected 4 atoms
MessageBoard State::MsgBrd;
MessageBoard State::MsgBrd2;
Camera State::camera;
char *State::ExeDir = NULL;	// set this up with calls to getenv()
int State::ImageIndex = 0;	// index of current displayed/processed AtomArray, [0, NumImages), & NumImages(overlap)
int State::NumImages = 0;	// total number of images for display
int State::Overlap1 = 0;	// remember 2 recent access structure indexes for overlap purpose
int State::Overlap2 = 0;
int State::PickMode = 1;	// global Pick Mode (0/1)
// towards eliminating this mode: you can always pick
	State::PickInfoType State::PickModeInfo = piNone;
int State::StereoMode = NONE_STEREO;
GLfloat State::Shine = 0.5f;	//0-1.0
bool State::ColorAnaglyph = false;
bool State::HasHardwareStereoSupport = false;

// draw backbone/loop in either line or tube mode: start,stop,draw
bool State::DrawBondAsLine = true;
bool State::DrawBackboneAsCurve = true;

// static variables of class Main
int Main::MovieMode = 0;
int Main::Spinning = 0;
int Main::Moving = 0;
int Main::Rocking = 0;
int Main::Translating = 0;
int Main::BeginX, Main::BeginY;
int Main::Scaling;
FQuat Main::CurQuat;
FQuat Main::LastQuat;
float Main::ScaleFactor = 1.0;
FVec3 Main::Translate;
DVec3 Main::BoundsMax;
DVec3 Main::BoundsMin;
int Main::MainWindow = -1;

//-----------------------------------------------------------------------------
int Main::main( int argc, char **argv ) {
	AtomArray *aaptr = NULL;
	char filename[FILENAME_MAXLEN], *filetype;

	NumImages = 0;
	for( int inFile = 1; inFile < argc; ) {
		if( argv[inFile][0] == '-' )
			filetype = argv[inFile++];
		else
			filetype = "";
		strncpy( filename, argv[inFile++], FILENAME_MAXLEN );
		SetValidPath( filename );
		ReadFileParams( filename, filetype, argv[0] );
	}
	ImageIndex = ( NumImages > 0 ) ? 0 : -1;	// display first image by default 
	trackball( CurQuat.q, 0.0, 0.0, 0.0, 0.0 );

	// Window init
	InitLightSource(  );

	//sphere list 
	for( int j = 0; j < NumImages; j++ )
		AddSphereToList(  );
	glFinish(  );

	return 0;
}

//-----------------------------------------------------------------------------

void Main::recalcmodelView( AtomArray *model ) {

	//TODO: the axes are being drawn in here for some reason, but they should be 
	// drawn elsewhere.  I removed the lighting disable/enable calls wrapping
	// those draws, so the axes may look weird until this is fixed.

	// view transformation already applied; model transform is now appended
	glMatrixMode( GL_MODELVIEW );

	GLfloat m[4][4];
	build_rotmatrix( m, CurQuat.q );
	glMultMatrixf( &m[0][0] );
		// rotate the scene (all models)

	if( camera.eye == camera.LEFT && gEnableAxes ) {
		double lw;
		glGetDoublev( GL_LINE_WIDTH, &lw );
		glLineWidth( 5.0 );
		float scale = BoundsMax.z ? BoundsMax.z / 2.f : 10.f;
		zglDrawStandardAxis( scale );
		glLineWidth( lw );
	}

	glTranslatef( Translate.x, Translate.y, Translate.z );
		// translate the scene (all models)

	if( model ) {
		glTranslatef( model->viewTranslate.x, model->viewTranslate.y, model->viewTranslate.z );
			// translate the model


		// translate to model origin: if CoM is already updated, it means that
		// translate and centroid will cancel (centroid translation has already been applied)
		// otherwise will translate to the centroid to rotate about the center of the model.
		DVec3 toCenter = model->centroid;
		toCenter.add( model->translate );

		glTranslated( toCenter.x, toCenter.y, toCenter.z );
	

		build_rotmatrix( m, model->viewQuat.q );
		glMultMatrixf( &m[0][0] );
			// rotate the model

		if( camera.eye == camera.LEFT && gEnableAxes ) {
				double lw;
				glGetDoublev( GL_LINE_WIDTH, &lw );
				glLineWidth( 2.0 );
				float scale = BoundsMax.z ? BoundsMax.z / 2.f : 10.f;
				zglDrawStandardAxis( scale );
				glLineWidth( lw );
		}

		glTranslated( -toCenter.x, -toCenter.y, -toCenter.z );
			// translate to model origin

	}

	glScalef( ScaleFactor, ScaleFactor, ScaleFactor );
		// scale everything
}

//-----------------------------------------------------------------------------

void Main::animate( void ) {
	AtomArray *aaptr = EFFECT_AAptr;
	
	if( NumImages <= 0 || !aaptr)
		return;

	if( Rocking > 0 ) {
		ax[0] = 0;
		ax[1] = ( Rocking == 1 ) ? 1 : -1;
		ax[2] = 0;
		angle = aaptr->InParm.rotateDegree;
		axis_to_quat( ax, angle, LastQuat.q );
		add_quats( LastQuat.q, CurQuat.q, CurQuat.q );
		Rocking = Rocking % 2 + 1;
	}
	else {
		if( MovieMode ) {
			ZTime elapsed = zTime - aaptr->lastStructAnimate;
			if( elapsed >= aaptr->InParm.sleepseconds / 1000.0 ) {
				// (sleepseconds, despite its name, appears to hold values in milliseconds)
				aaptr->lastStructAnimate = zTime;
				for( int k = ( NO_OVERLAP ) ? ImageIndex : Overlap1;; ) {
					aaptr = AAptrs[k];
					if( aaptr->InParm.filetag == fPTH || aaptr->InParm.filetag == fDCD ) {
						aaptr->ReadNextStruct(  );	// only path and dcd support movie
						aaptr->lastStructAnimate = zTime;
					}
					if( NO_OVERLAP )
						break;
					else if( k != Overlap2 )
						k = Overlap2;
					else
						break;
				}
			}
		}
		if( Spinning == 1 ) {

			if( !manualAlign ) {
				add_quats( LastQuat.q, CurQuat.q, CurQuat.q );
			}
			else {
				add_quats( LastQuat.q, manualAlign->viewQuat.q, manualAlign->viewQuat.q );
			}
		}
	}

	if( ( MovieMode || Rocking ) && ( PickMode || 1 ) ) {
		PickQ.print_pick_info(  );
	}
}

//-----------------------------------------------------------------------------

void Main::motion( int x, int y, float W, float H ) {
	// This is called by mouse message handlers in the _cmoil.cpp plugin file.
	// See ZUIRenderCmoil::handleMsg()

	if( NumImages <= 0 )
		return;

	// experimental: allow the motion to modify the transform of a particular 
	// model if desired:
	FVec3 *pTrans  = &Translate;
	if( manualAlign ) {
		pTrans = &manualAlign->viewTranslate;
	}

	if( Moving ) {
		// should be called "Rotating"
		W *= 0.5;
		H *= 0.5;
		if( !manualAlign ) {
			//trace( "screen: p1 = (%d,%d) p2 = (%d,%d) delta=(%d,%d)\n", BeginX, BeginY, x, y, x-BeginX, y-BeginY );
			trackball( LastQuat.q, ( BeginX - W ) / W, ( H - BeginY ) / H, ( x - W ) / W, ( H - y ) / H );
		}
		else {
			// the LastQuat above is in the fixed global coordinate system.  If we are modifiy the rotation
			// of an individual model, we need to transform the axis into the local system.

			FMat4 rotLocal   = manualAlign->viewQuat.mat();
			FMat4 rotGlobal  = CurQuat.mat();
			trackball( LastQuat.q, ( BeginX - W ) / W, ( H - BeginY ) / H, ( x - W ) / W, ( H - y ) / H,
				&rotLocal, &manualAlign->viewTranslate, &rotGlobal, &Translate );
		}
		BeginX = x;
		BeginY = y;
		Spinning = 1;
	}
	if( Translating ) {
		float scaleBias =(float) AAptrs[0]->maxDisplacementXYZ;
		FVec3 trans(  scaleBias * float( x - BeginX ) / float ( W ),
					  -scaleBias * float( y - BeginY ) / float ( H ),
					  0.f );

		// transform this displacement according to current rotation
		FMat4 transMat;
		build_rotmatrix( transMat.m, CurQuat.q );
		transMat.transpose();
		trans = transMat.mul( trans );

		//trace( "--> %g, %g, %g\n", trans.x, trans.y, trans.z );

		pTrans->x += trans.x;
		pTrans->y += trans.y;
		pTrans->z += trans.z;
		//trace( "translate by %6.2f %6.2f %62.f\n", pTrans->x, pTrans->y, pTrans->z );

		BeginX = x;			// save last location
		BeginY = y;			// save last location
	}
}


//-----------------------------------------------------------------------------

float Main::ax[3];
float Main::angle = 0;
void Main::animate_rot(  ) {
	// unused currently?
	// tfb: feb 2009

	if( MainWindow < 0 || ImageIndex < 0 )
		return;
	axis_to_quat( ax, angle, LastQuat.q );
	add_quats( LastQuat.q, CurQuat.q, CurQuat.q );
}

//-----------------------------------------------------------------------------
// (x,y,z)Axis= 0, 1, -1
void Main::Rotates( float xAxis, float yAxis, float zAxis )	// deltaX=0, int deltaY=0, int deltaZ=0)
{
	// unused currently?
	// tfb: feb 2009

	if( MainWindow < 0 || ImageIndex < 0 )
		return;

	AtomArray *aaptr = AAptrs[NO_OVERLAP ? ImageIndex : Overlap1];
	ax[0] = -xAxis;
	ax[1] = -yAxis;
	ax[2] = -zAxis;
	if( ax[0] != 0 || ax[1] != 0 || ax[2] != 0 ) {
		angle = aaptr->InParm.rotateDegree;
		axis_to_quat( ax, angle, LastQuat.q );
		add_quats( LastQuat.q, CurQuat.q, CurQuat.q );
	}
}

//-----------------------------------------------------------------------------
int Main::ZclipFront = 0, Main::ZclipBack = 0;
void Main::animate_zclip(  ) {
	if( MainWindow < 0 || ImageIndex < 0 )
		return;
	camera.Clip( ZclipFront, ZclipBack );
}

//-----------------------------------------------------------------------------
// front, back = -1,0,1,   -2:reset,  2: on|off
void Main::ZClip( int front, int back ) {
	if( front == back ) {
		if( front == ClipOnOff ) {
			if( ZclipFront != NoClip || ZclipBack != NoClip ) {
				ZclipFront = NoClip;
				ZclipBack = NoClip;
			}
			else {
				ZclipFront = ClipOnOff;
				ZclipBack = ClipOnOff;
			}
		}
		else if( front == ResetClip ) {
			camera.Clip( ResetClip, ResetClip );
			ZclipFront = NoClip;
			ZclipBack = NoClip;
		}
	}
	else {
		ZclipFront = front;
		ZclipBack = back;
		camera.Clip( front, back );
	}
	if( MainWindow < 0 || ImageIndex < 0 )
		return;
}

//-----------------------------------------------------------------------------
void Main::display_clip(  ) {
	// clip front and back of the screen picture
	if( ZclipFront != 0 || ZclipBack != 0 ) {
		glMatrixMode( GL_MODELVIEW );
		GLdouble m[16];
		glGetDoublev( GL_MODELVIEW_MATRIX, m );
		double z = BoundsMax.z * m[10];	// convert maxz to eye coordinate
		double stepDelta = z * 0.02;
		if( camera.clipBack != 0 )	// back cut
		{
			GLdouble eqn[4] = { 0, 0, 1, z - camera.clipBack * stepDelta };
			glClipPlane( GL_CLIP_PLANE1, eqn );
			glEnable( GL_CLIP_PLANE1 );
		}
		if( camera.clipFront > 0 )	// front cut
		{
			GLdouble eqn[4] = { 0, 0, -1, z - camera.clipFront * stepDelta };
			glClipPlane( GL_CLIP_PLANE0, eqn );
			glEnable( GL_CLIP_PLANE0 );
		}
	}
}

//-----------------------------------------------------------------------------
unsigned char gGlobalAlpha;
void Main::display_draw( int windowSize ) {
	INIT_TIMERS( 5 );
	if( ( nTimed < 5 || globalTimerOn ) )
	trace( " *** %9.8lf elapsed to display_draw\n", zTimeNow() - globalTime );


	START_TIMER( display_draw );
	START_TIMER( display_draw_START );
	gGlobalAlpha=255;
	// this is set per imageindex to allow us to easily control global alpha setting
	// per image, to aid in viewing multiple images at once (allow user to twiddle)

	AtomArray *aaptr;
	int colorFilter =
		( StereoMode !=
		ANAGLYPH_STEREO ) ? ColorFilterNone : ( ( camera.eye == camera.LEFT ) ? ColorFilterRed : ColorFilterCyan );

	display_clip(  );

	SphereShiness(  );

	int ii = 0;
	int displayedStructs[64];
	int displayedStructCount = 0;


	END_TIMER( display_draw_START );
	ZTime timeNow = zTimeNow();

#define BLINK_PERIOD 1.0

	tRGBA *colors[5] = {&orange, &purple, &green, &blue, &red };
		// quick hack to dispaly things in different colors

	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
	glEnable( GL_BLEND );
	glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
		// this is used for the alpha scaling of different structures, but is 
		// this causing a serious perf hit?

	int displayCount=0;
	for( int k=0; k<NumImages; k++ ) {
		aaptr = AAptrs[k];
		if( aaptr && aaptr->bDraw ) {
			displayCount++;
//			int colorIndex = aaptr->properties.getI( "colorIndex", -1 );
			tRGBA *aaColor = displayCount == 1 ? &orange : &purple;
//			if( colorIndex >= 0 ) {
//				aaColor = &zmoilColors[ colorIndex ];
//			}
			START_TIMER( display_draw_PRE_DisplayStruct );
			glPushMatrix();
			recalcmodelView( aaptr );

			// NEW: loopp over number of struct we want to draw from this PTH if 
			// desired
			int structIndex = aaptr->currentmol-1;
			double skipStructs = 999999.0;
			if( Zmoil_PthDisplayStructs != 1 && aaptr->InParm.structno > 1 ) {
				Zmoil_PthDisplayStructs = min( Zmoil_PthDisplayStructs, aaptr->InParm.structno / 2 );
				skipStructs = (double)(aaptr->InParm.structno-1) / double(Zmoil_PthDisplayStructs-1);
				structIndex = 0;
					// we always start with 1st frame when displaying multiple
			}
			int structCount = 0;
			int currentMolSaved = aaptr->currentmol;
			for( ; structIndex < aaptr->InParm.structno; structIndex = int((double)structCount * skipStructs) ) {  
				structCount++;
				if( structIndex != aaptr->currentmol-1 ) {
					aaptr->ReadStruct( structIndex );
				}
				if( Zmoil_PthDisplayStructs > 1 ) {
					if( Zmoil_PthDisplayAlphaLinear ) {
						gGlobalAlpha = int((float)structCount/(float)Zmoil_PthDisplayStructs * 255.f);
							// linear
					}
					else {
						gGlobalAlpha = float(structCount * structCount) / (float)( Zmoil_PthDisplayStructs * Zmoil_PthDisplayStructs ) * 255.f;
							// quadratic
					}
				}

				ZTime elapsed = 0;
				if( aaptr->bCycle || aaptr->bBlink ) {
						// note bCycle has been tacked on to bBlink functionality; if cycle gets turned on
						// then blinking is disabled & we cycle through alternate models instead
					elapsed = timeNow - aaptr->lastBlinkOn;
					if( elapsed > BLINK_PERIOD ) {
						aaptr->lastBlinkOn += elapsed;
						elapsed -= BLINK_PERIOD;
						if( aaptr->bCycle ) {
							aaptr->DisplayNextPdbModel();
						}
					}
					// aaptr->alpha = fabs( elapsed - BLINK_PERIOD/2.0 ) / (BLINK_PERIOD/2.0) * 255;
					// experiment using alpha fades for blink
				}
				END_TIMER( display_draw_PRE_DisplayStruct );
				if( aaptr->bCycle || elapsed < BLINK_PERIOD / 2.0 ) {
					if( Zmoil_PthDisplayStructs == 1 ) {
						gGlobalAlpha = aaptr->alpha;
					}

					START_TIMER( display_draw_callDisplayStruct );
					if( ( nTimed < 5 || globalTimerOn ) )
					trace( " *** %9.8lf elapsed to before call to DisplayStructure\n", zTimeNow() - globalTime );

					//
					// DisplayStructure will render the struct as ball/stick/spacefill etc
					//
					aaptr->DisplayStructure( ii, colorFilter, windowSize );
					displayedStructs[displayedStructCount++] = aaptr->currentmol;

					if( ( nTimed < 5 || globalTimerOn ) )
					trace( " *** %9.8lf elapsed after call to DisplayStructure\n", zTimeNow() - globalTime );

					END_TIMER( display_draw_callDisplayStruct );


					START_TIMER( display_draw_POST_DisplayStruct );

					ZTime bbone;
					int time2nd = ( aaptr->dispmodeComb & (RNABACKBONE|BACKBONE|STRUCT2nd|RIBBON) ); 
					if( time2nd && ( nTimed < 5 || globalTimerOn ) ) {
						bbone = zTimeNow();
					}

					//
					// backbone
					//
					if( ( aaptr->dispmodeComb & BACKBONE ) != 0 ) {
						DrawBackboneAsCurve ? glDisable( GL_LIGHTING ) : glEnable( GL_LIGHTING );
						aaptr->backbone( windowSize, colorFilter, BACKBONE );
					}
					if( ( aaptr->dispmodeComb & RNABACKBONE ) != 0 ) {
						DrawBackboneAsCurve ? glDisable( GL_LIGHTING ) : glEnable( GL_LIGHTING );
						aaptr->backbone(windowSize, colorFilter, RNABACKBONE );
					}
					
					// 
					// 2nd-struct OR ribbon, never both
					//
					if( ( aaptr->dispmodeComb & STRUCT2nd ) != 0 ) {
						glEnable( GL_LIGHTING );
						aaptr->struct2nd( windowSize, colorFilter );
					}
					else if( ( aaptr->dispmodeComb & RIBBON ) != 0 ) {
						glEnable( GL_LIGHTING );
						aaptr->ribbon( windowSize, colorFilter );
					}

					//
					// velocity vectors
					//
					if( Zmoil_VelocityDisplay && aaptr->properties.has( "velocityMax" ) ) {
						glDisable( GL_LIGHTING );
						glLineWidth( 0.1f );
						aaptr->velocityVectors();
						glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
					}

					if( time2nd && ( nTimed < 5 || globalTimerOn ) ) {
						bbone = zTimeNow() - bbone;
						trace( "    [%9.8lf secs] -->  display_draw, backbone, 2nd, etc\n", bbone );
					}

					if( Zmoil_PeptidePlaneVectors ) {
						glDisable( GL_LIGHTING );
						aaptr->peptidePlaneVectors();
					}

					if( ( aaptr->dispmodeComb & SURF ) != 0 ) {
						if( ( nTimed < 5 || globalTimerOn ) ) {
							bbone = zTimeNow();
						}

						// draw ravi normals
						if( aaptr->InParm.dispmode & RAVINORM ) {
							// glEnable( GL_LIGHTING );
							// aaptr->raviAtomDraw( colorFilter );
							// glDisable( GL_LIGHTING );
							glDisable( GL_LIGHTING );
							aaptr->raviNormalDraw( colorFilter );
						}
						if( aaptr->InParm.dispmode & RAVINORMCOARSE ) {
							glDisable( GL_LIGHTING );
							aaptr->raviNormalDraw( colorFilter, true );
						}
						if( aaptr->surface ) {
							tRGBA *surfColor = colors[k%5];
							glEnable( GL_LIGHTING );
							//glEnable( GL_BLEND );
							//glBlendFunc( GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA );
							aaptr->surf( colorFilter, aaptr->InParm.uniqueSurfaceColors ? surfColor : 0 );
							//glDisable( GL_BLEND );
						}
						if( ( nTimed < 5 || globalTimerOn ) ) {
							bbone = zTimeNow() - bbone;
							trace( "    [%9.8lf secs] -->  display_draw, surface + pick sphere\n", bbone );
						}
					}


					// draw sphere net
					//if ( (aaptr->pickSphereCenter)>=0 && aaptr->InParm.showSphereNet )
					//{
					//  aaptr->SphereNet(colorFilter);
					//}
					// the above is mostly useful for debugging: the sphere radius is in Angstroms, which helps with scale

					glEnable( GL_LIGHTING );
					aaptr->DisplayPickedResidue( *aaColor, colorFilter );

					END_TIMER( display_draw_POST_DisplayStruct );
				}
			}
			if( aaptr->currentmol != currentMolSaved ) {
				aaptr->ReadStruct( currentMolSaved - 1 );
			}
			glPopMatrix();
		}
	}

	glDisable( GL_COLOR_MATERIAL );
	glDisable( GL_BLEND );
	glDisable( GL_DEPTH_TEST );


	START_TIMER( display_draw_Messages );
	gGlobalAlpha = 255;
	PickQ.display_pick( windowSize, colorFilter );
	glDisable( GL_CLIP_PLANE1 );	// if any clip on, 
	glDisable( GL_CLIP_PLANE0 );
	if( gEnableStructInfoDisplay ) {
		MsgBrd.print( 1 );		// 1:  bg = black, do not apply clip
	}

	// Second MessageBoard added for lower-righth-hand display: currently always displays index of displayed struct
	if( gEnableStructInfoDisplay ) {
		char buffer[512];
		strcpy( buffer, "Index: " );
		for( int n = 0; n < displayedStructCount; n++ ) {
			char buf[8];
			sprintf( buf, "%d", displayedStructs[n] );
			strcat( buffer, buf );
			if( n + 1 < displayedStructCount ) {
				strcat( buffer, ", " );
			}
		}
		MsgBrd2.erase(  );
		if( displayedStructCount ) {
			MsgBrd2.set( buffer );
			MsgBrd2.print( 1, 1 );
			// print lower right corner
		}
	}
	END_TIMER( display_draw_Messages );
	if( nTimed<5 ) {
		trace( " *** %9.8lf elapsed end of display_draw\n", zTimeNow() - globalTime );
	}
	END_TIMER( display_draw );
	TRACE_TIMERS();
}							// display_draw()


//-----------------------------------------------------------------------------
void Main::display( int wWidth, int wHeight ) {

	INIT_TIMERS( 5 );
	if( ( nTimed < 5 || globalTimerOn ) )
	trace( " *** %9.8lf elapsed to MainDisplay\n", zTimeNow() - globalTime );

	START_TIMER( MainDisplayStartup );

	int windowSize;
	glDrawBuffer( GL_BACK_LEFT );
	if( NumImages <= 0 ) {
		glClearColor( 0, 0, 0, 1 );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		return;
	}

	glClearColor( AAptrs[0]->InParm.bg.r / 255., AAptrs[0]->InParm.bg.g / 255., AAptrs[0]->InParm.bg.b / 255., 1 );
	glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );

	if( HasHardwareStereoSupport ) {
		glDrawBuffer( GL_BACK_RIGHT );
		glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
	}
	if( StereoMode == ANAGLYPH_STEREO ) {
		glClearAccum( 0.0, 0.0, 0.0, 0.0 );
		glClear( GL_ACCUM_BUFFER_BIT );
	}


	// classify window to small(1), middle(2), large(3)
	if( wWidth < 350 || wHeight < 350 )
		windowSize = 1;
	else if( wWidth < 600 || wHeight < 600 )
		windowSize = 2;
	else
		windowSize = 3;
	END_TIMER( MainDisplayStartup );

	START_TIMER( MainDisplayCamera );
	camera.SetSizes( wWidth, wHeight, BoundsMax.z, BoundsMin.z );
	camera.LookAt( camera.LEFT, StereoMode );
	END_TIMER( MainDisplayCamera );

	if( nTimed<5 ) {
		trace( " *** %9.8lf elapsed about to call display_draw\n", zTimeNow() - globalTime );
	}
	display_draw( windowSize );
	if( nTimed<5 ) {
		trace( " *** %9.8lf elapsed back from display_draw\n", zTimeNow() - globalTime );
	}


	START_TIMER( MainDisplayStereo );
	if( StereoMode != NONE_STEREO ) {
		if( StereoMode == ANAGLYPH_STEREO ) {
			glAccum( GL_LOAD, 1.0 );
			glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
		}
		camera.LookAt( camera.RIGHT, StereoMode );
		display_draw( windowSize );

		if( StereoMode == ANAGLYPH_STEREO ) {
			glAccum( GL_ACCUM, 1.0 );
			glAccum( GL_RETURN, 1.0 );
		}
	}
	if( StereoMode != HARDWARE_STEREO && HasHardwareStereoSupport ) {
		glAccum( GL_LOAD, 1.0 );
		glDrawBuffer( GL_BACK_RIGHT );
		glAccum( GL_RETURN, 1.0 );
	}
	END_TIMER( MainDisplayStereo );
	if( nTimed<5 ) {
		trace( " *** %9.8lf elapsed end of Main::Display\n", zTimeNow() - globalTime );
	}
	TRACE_TIMERS();
}							// display()

//-----------------------------------------------------------------------------

GLUquadricObj *Main::qobj = NULL;
/* Add to sphere List */
void Main::AddSphereToList(  ) {
	if( NumImages <= 0 )	/* check last image to add missing sphere to list */
		return;

	AtomArray *aaptr = AAptrs[NumImages - 1];
	Atom *listi = aaptr->atomList;
	int i, l;
	Coordinate crd;

	if( !glIsList( 1 ) )	// add unit sphere, redius differ from Stick and backbone Tube
	{
		glNewList( 1, GL_COMPILE );
		gluSphere( qobj, /*Radius */ Stick::StickRadius * 3., /*Slices */ 8, /*Stacks */ 8 );
			// unit sphere, typically used for ball/stick display, or picked atoms is 
			// less polys to help frame rate
		glEndList(  );
	}


	for( i = 0; i < aaptr->numAtom; i++, listi++ ) {
		l = int ( listi->radius * 100 + 0.5 );
		if( !glIsList( l ) ) {
			//printf( "building glList for radius %g\n", listi->radius );
			glNewList( l, GL_COMPILE );	// Create sphere display list.
			gluSphere( qobj, listi->radius, 20, 20 );
			glEndList(  );
		}
	}
	Tube::AddToQuadricList( qobj );
}

//-----------------------------------------------------------------------------
void Main::SphereShiness(  ) {
	// specular color
	GLfloat sp[4] = { Shine, Shine, Shine, 1.0 };
	glMaterialfv( GL_FRONT_AND_BACK, GL_SPECULAR, sp );
	glMaterialf( GL_FRONT_AND_BACK, GL_SHININESS, 100. );
}

//-----------------------------------------------------------------------------
void Main::InitLightSource(  ) {
	qobj = gluNewQuadric(  );
	gluQuadricDrawStyle( qobj, GLU_FILL );

	// diffuse color 0.75
	GLfloat light_diffuse[4] = { 0.75, 0.75, 0.75, 1.0 };
	glLightfv( GL_LIGHT0, GL_DIFFUSE, light_diffuse );

	GLfloat light_position[4] = { 1.0, 1.0, 1.0, 0.0 };
	//glMatrixMode( GL_MODELVIEW );
	//glPushMatrix(  );
	//glLoadIdentity(  );
	glLightfv( GL_LIGHT0, GL_POSITION, light_position );
	//glPopMatrix(  );

	glLightf( GL_LIGHT0, GL_SPOT_EXPONENT, 5.0 /* *ScaleFactor */  );

	GLfloat sp[4] = { 1.0, 1.0, 1.0, 1.0 };
	glLightfv( GL_LIGHT0, GL_SPECULAR, sp );

	// ambient color
	GLfloat tv[4] = { 0.5, 0.5, 0.5, 1.0 };	// grey color 
	glLightModelfv( GL_LIGHT_MODEL_AMBIENT, tv );
	glLightModelf( GL_LIGHT_MODEL_TWO_SIDE, 0.0 );

	SphereShiness(  );

	GLfloat clr[] = { 0.1f, 0.5f, 0.8f, 1.0f };
	glMaterialfv( GL_BACK, GL_AMBIENT_AND_DIFFUSE, clr );
}


//-----------------------------------------------------------------------------
// set globals: ImageIndex  NumImages AAptrs[ImageIndex] 
int Main::ReadFileParams( char *fileName, char *fileType, char *cmoiexe ) {
	AtomArray *aaptr = NULL;
	int currCrd = 0;

	if( fileName == NULL )
		return 0;

	aaptr = new AtomArray;
	if( aaptr == NULL ) {
		fprintf( stderr, "Error: Unable to allocate memory allocation at ReadFileParams().\n" );
		exit( -1 );
	}
	strcpy( aaptr->InParm.fileName, fileName );
	aaptr->InParm.pdbName[0] = 0;
	if( strcmp( fileType, "-pdb" ) == 0 ) {
		strcpy( aaptr->InParm.crdName, fileName );	// no instruction file
		aaptr->InParm.filetag = fPDB;
	}
	else if( strcmp( fileType, "-xyz" ) == 0 ) {
		strcpy( aaptr->InParm.crdName, fileName );	// no instruction file
		aaptr->InParm.filetag = fXYZ;
	}
	else if( strcmp( fileType, "-crd" ) == 0 ) {
		strcpy( aaptr->InParm.crdName, fileName );	// no instruction file;
		// we will try to use a similarly named wcon file, if it exists.
		char buf[256];
		strcpy( buf, replaceExtension( fileName, "wcon" ) );
		if( zWildcardFileExists( buf ) ) {
			strcpy( aaptr->InParm.connName, buf );
			trace( "Found wcon file %s and will try to use it.\n", buf );
		}
		else {

		}
		aaptr->InParm.filetag = fCRD;
	}
	else if( strcmp( fileType, "-dcd" ) == 0 ) {
		// we need to be able to find a .wcon to load this natively, 
		int haveWCON = 0;
		char buf[256];
		strcpy( buf, fileName );
		char *p = strstr( buf, "_dyn.dcd" );
		if( p ) {
			*p = 0;
			strcat( buf, ".wcon" );
		}
		else {
			strcpy( buf, replaceExtension( buf, "wcon" ) );
		}
		if( zWildcardFileExists( buf ) ) {
			strcpy( aaptr->InParm.connName, buf );
			strcpy( aaptr->InParm.dcdName, fileName );
			aaptr->InParm.swapmode = 1;
			aaptr->InParm.filetag = fDCD;
			aaptr->InParm.structno = 1000; // default to load all structs
		}
		else {
			ZTmpStr message( "If you want to load a DCD file natively, you should ensure that %s exists with the DCD file.  Otherwise, use the MOIL graphical interface to View Structure, from which you can specify the DCD and WCON to use.", buf );
			zuiMessageBox( "Connectiving File", message.s, 1 );
			delete aaptr;
			return 0;
		}
	}
	else
		currCrd = aaptr->GetInputParameters( fileName, currCrd, cmoiexe );

	if( currCrd >= 0 ) {
		ImageIndex = NumImages;
		AAptrs[ImageIndex] = aaptr;
		NumImages++;
		aaptr->ProcessData(  );
		return 1;
	}
	else
		delete aaptr;
	return 0;
}

}
