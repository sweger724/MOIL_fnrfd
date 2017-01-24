//*************************************************************************************************
//*  Filename:   cmoil_pick.cpp
//*
//*  Description: 
//*    Methods of PickedAtoms class for disply information of selected atoms
//*
//*  Modification History:
//*  
//*  Date           Developer   Description
//*  ------------   ----------- ---------------------------------------------------------------------
//*  Aug. 2000  ??      Initial Development
//*  Mar. 2001    Baohua Wang Seperated from cmoil.cpp and modified for message display
//*  July 2001  Baohua Wang Reshape to support mutiple structure display
//*  Sept 2007      Thomas Blom     Update for Linux/OSX
//*************************************************************************************************
#ifdef __APPLE__
#include "sys/uio.h"			// for ssize_t
#include "pthread.h"
#endif


#include <stdlib.h>
#include "cmoil_pick.h"
#include "cmoil.h"
#include "cmoil_globals.h"


namespace CMOIL {
#define IsOverlap    (image2 >= 0)

PickedAtoms::PickedAtoms(  ) {
	currImage = 0;
	reset(  );
} 

void PickedAtoms::reset( int imageIndex1, int imageIndex2 )	// should pass in two overlaped image index also
{
	PickModeInfo = piNone;
	image1 = imageIndex1;
	image2 = imageIndex2;

	// Reset the display flags for any atoms that were picked, and clear numentries:
	for( int i = 0; i < numentries; i++ ) {
		aptr[i]->dispmode &= ~UNITBALL;
		aaptr[i] = 0;
	}
	numentries = 0;
}

int PickedAtoms::isPicked( Atom * atom ) {
	for( int i = 0; i < numentries; i++ ) {
		if( aptr[i] == atom )
			return 1;
	}
	return 0;
}

// set at most 3 recent picked Atom Crds, return number of picked crds
int PickedAtoms::getPickedCrds( Coordinate * crds ) {
	int numPicked = ( numentries > 3 ) ? 3 : numentries;
	for( int i = 0; i < numPicked; i++ )
		crds[i] = *aptr[i];
	return numPicked;
}

// mouse clicked at (x,y), return the select atom number.
// Since GL_SELECT render mode is not supported by hardware acceleration,
// implement the select in GL_RENDER mode.
//
int PickedAtoms::get_atomnum( int x, int y, int & imgIndex ) {

	// Note: for picking we assume all of the loaded models are using the same
	// transform -- that is, NOT using the facility to manually align models.
	if( Main::AAptrs[0] ) {
		glPushMatrix();
		Main::recalcmodelView( Main::AAptrs[0] );
	}

	GLint viewport[4];
	GLdouble pM[16], mM[16];

	glGetIntegerv( GL_VIEWPORT, viewport );
	glGetDoublev( GL_PROJECTION_MATRIX, pM );
	glGetDoublev( GL_MODELVIEW_MATRIX, mM );

	// debug: look at matrices etc.
	/*
	   double m[16];
	   FILE *f = fopen( "modelview.txt", "w" );
	   glGetDoublev( GL_MODELVIEW_MATRIX, m );
	   fprintf( f, "ModeView:\n" );
	   for( int i=0; i<4; i++ ) {
	   for( int j=0; j<4; j++ ) {
	   fprintf( f, "\t%g", m[i*4+j] );
	   }
	   fprintf( f, "\n" );
	   }
	   glGetDoublev( GL_PROJECTION_MATRIX, m );
	   fprintf( f, "\nProjection:\n" );
	   for( i=0; i<4; i++ ) {
	   for( int j=0; j<4; j++ ) {
	   fprintf( f, "\t%g", m[i*4+j] );
	   }
	   fprintf( f, "\n" );
	   }
	   fprintf( f, "\nViewport:\n" );
	   for( i=0; i<4; i++ ) {
	   fprintf( f, "\t%d", viewport[i] );     
	   }
	   fprintf( f, "\n\n" );
	   fclose(f);
	 */


	y =  ( viewport[3] + viewport[1] ) - ( GLint ) y - 1;
	x += viewport[0] ;
	double minz = 0., ox, oy, oz, wx, wy, wz;
	double minR2 = 900, R2;

	Atom *a = NULL;			// , *lastA=NULL;
	imgIndex = -1;
	int atomIndex = -1, structN, atomCnt = 0;
	for( structN = 0; structN <NumImages; structN++ ) {
		a = AAptrs[structN]->atomList;
		atomCnt = AAptrs[structN]->numAtom;

		for( int i = 0; i < atomCnt; i++, a++ ) {
			if( a->skip == 0 && gluProject( a->x, a->y, a->z, mM, pM, viewport, &wx, &wy, &wz ) ) {
				if( wz >= 0 ) {
					wx += sqrt( ( wx - x ) * ( wx - x ) + ( wy - y ) * ( wy - y ) );	// get distance in window space
					gluUnProject( wx, wy, wz, mM, pM, viewport, &ox, &oy, &oz );	// map to object space
					R2 = sqrt( ( a->x - ox ) * ( a->x - ox ) + ( a->y - oy ) * ( a->y - oy ) + ( a->z - oz ) * ( a->z - oz ) ) / a->radius;	// compare tdistance to radius
					if( R2 < minR2 || ( R2 == minR2 && wz < minz ) ) {
						minR2 = R2;
						atomIndex = i;
						minz = wz;
						imgIndex = structN;
					}
				}
			}
		}
	}
	if( Main::AAptrs[0] ) {
		glPopMatrix();
	}
	return atomIndex;
}

// put selected atom in to picked list
void PickedAtoms::put( int x, int y ) {
	int atomIndex, imgIndex;
	atomIndex = get_atomnum( x, y, imgIndex );

	if( atomIndex >= 0 ) {
		// UNSELECT oldest in queue if queue is full 
		if( numentries == ATOM_MAX_PICK ) {
			aptr[ATOM_MAX_PICK-1]->dispmode &= ( ~UNITBALL );
		}

		// SHIFT queue
		for( int j = numentries; j >= 1; j-- ) {	
			q[j]     = q[j - 1];
			aaptr[j] = aaptr[j - 1];
			aptr[j]  = aptr[j - 1];
		}

		// SET new entry
		q[0]     = atomIndex;
		aaptr[0] = AAptrs[imgIndex];
		aptr[0]  = aaptr[0]->atomList + atomIndex;
		aptr[0]->dispmode |= UNITBALL;

		if( numentries < ATOM_MAX_PICK )
			numentries++;

		// print coordinates only for XYZ type crd file
		//
		print_atom_info(  );
	}

}							// PickedAtoms::put()

// get full atom name for display, aaptr[atomIndex] -- error?
char *PickedAtoms::atom_name( int atomIndex, char fileTag ) {
	static char name[160];

	if( atomIndex > 3 || atomIndex < 0 )	// support 0-3 index at same time
	{
		name[0] = '\0';
		return name;
	}

	int atomNum = aptr[atomIndex] - aaptr[atomIndex]->atomList + 1;
	char *myname = name + atomIndex * 40;
	if( fileTag == fPDB ) {
		Model *model = aaptr[atomIndex]->models + aaptr[atomIndex]->chains[ aptr[atomIndex]->chainIndex ].modelIdx;
		sprintf( myname, "%s:%.0d:%c:%d:%s:%d:%s", aaptr[atomIndex]->moleculeName, model->modelName,	// display model, chain for PDB file
				 ( aaptr[atomIndex]->chains + aptr[atomIndex]->chainIndex )->id,
				 aaptr[atomIndex]->resNum( aptr[atomIndex] ),
				 aaptr[atomIndex]->resName[aptr[atomIndex]->rn], atomNum, aptr[atomIndex]->atomName );
	}
	else if( fileTag == fPTH || fileTag == fDCD )	// binary crd to display structure #
	{
		sprintf( myname, "%s:%0d:%d:%s:%d:%s", aaptr[atomIndex]->moleculeName, aaptr[atomIndex]->currentmol,	// display currentmol for pth or dcd file 
				 aaptr[atomIndex]->resNum( aptr[atomIndex] ),
				 aaptr[atomIndex]->resName[aptr[atomIndex]->rn], atomNum, aptr[atomIndex]->atomName );
	}
	else {
		sprintf( myname, "%s:%d:%s:%d:%s",
				 aaptr[atomIndex]->moleculeName,
				 aaptr[atomIndex]->resNum( aptr[atomIndex] ),
				 aaptr[atomIndex]->resName[aptr[atomIndex]->rn], atomNum, aptr[atomIndex]->atomName );
	}
	return myname;
}
// atom_name()

// print info for last selected atom
void PickedAtoms::print_atom_info(  ) {
	char text[LINE_MAXLEN];	// temp for Footer display
	if( aaptr[0]->InParm.filetag == fXYZ )
		sprintf( text, "(%.5g, %.5g, %.5g)", aptr[0]->orgX, aptr[0]->orgY, aptr[0]->orgZ );
	else
		sprintf( text, "%s (%.5g, %.5g, %.5g) Radius=%.5gA",
				 atom_name( 0, ( aaptr[0]->InParm.filetag ) ),
				 aptr[0]->orgX, aptr[0]->orgY, aptr[0]->orgZ, aptr[0]->radius );

	MsgBrd.set( text );
}

// print the distance btwn first two pick entries 
void PickedAtoms::print_distance(  ) {
	Coordinate t, t1, t2;
	double d;
	char text[LINE_MAXLEN];

	if( numentries < 2 ) {
		MsgBrd.set( "Please select two atoms for a distance\n" );
		return;
	}

//	t1 == *aptr[0];
//	t2 == *aptr[1];

	t1 = *aptr[0];
	t2 = *aptr[1];
		// tfb: the first (original) method uses operator== (!) to copy the original
		// coordinates; but this means aligned structures will not report correct
		// inter-atomic distances; so I changed to dynamic coords; this shouldn't matter
		// anyway, since all atoms are transformed equally under rotations etc...

	t = t1 - t2;
	d = t.len(  );

	char hasPDB = aaptr[0]->InParm.filetag;
	if( aaptr[0]->InParm.filetag == fPDB || aaptr[1]->InParm.filetag == fPDB )
		hasPDB = 'b';
	else
		hasPDB = aaptr[0]->InParm.filetag;	// set to the crd filetag
	sprintf( text, "Distance of (%s, %s): %.5g A", atom_name( 0, hasPDB ), atom_name( 1, hasPDB ), d );

	MsgBrd.set( text );

	PickModeInfo = piDistance;

}							// print_distance()

void PickedAtoms::print_angle(  ) {
	Coordinate t1, t2;
	Coordinate c0, c1, c2;
	double ac;
	char text[LINE_MAXLEN];

	if( numentries < 3 ) {
		MsgBrd.set( "Please select three atoms for an angle" );
		return;
	}

//	c0 == *aptr[0];
//	c1 == *aptr[1];
//	c2 == *aptr[2];

	c0 = *aptr[0];
	c1 = *aptr[1];
	c2 = *aptr[2];
		// see comments above in ::print_distance (tfb)

	t1 = ( c0 - c1 ).norm(  );
	t2 = ( c2 - c1 ).norm(  );
	ac = acos( t1 ^ t2 ) * 180 / PI;	//dot

	char hasPDB;
	if( aaptr[0]->InParm.filetag == fPDB || aaptr[1]->InParm.filetag == fPDB || aaptr[2]->InParm.filetag == fPDB )
		hasPDB = fPDB;
	else
		hasPDB = aaptr[0]->InParm.filetag;

	sprintf( text, "Angle of (%s, %s, %s): %.5g Deg.",
			 atom_name( 0, hasPDB ), atom_name( 1, hasPDB ), atom_name( 2, hasPDB ), ac );

	MsgBrd.set( text );

	PickModeInfo = piAngle;

}

void PickedAtoms::print_torsional(  ) {
	char text[LINE_MAXLEN];

	if( numentries < 4 ) {
		MsgBrd.set( "Please select four atoms for a torsinal angle" );
		return;
	}

	char hasPDB;
	if( aaptr[0]->InParm.filetag == fPDB || aaptr[1]->InParm.filetag == fPDB ||
		aaptr[2]->InParm.filetag == fPDB || aaptr[3]->InParm.filetag == fPDB )
		hasPDB = fPDB;
	else
		hasPDB = aaptr[0]->InParm.filetag;

	sprintf( text, "Tors.Angle of (%s, %s, %s, %s): %.5g Deg.",
			 atom_name( 0, hasPDB ), atom_name( 1, hasPDB ),
			 atom_name( 2, hasPDB ), atom_name( 3, hasPDB ), angleT( *aptr[0], *aptr[1], *aptr[2], *aptr[3] ) );

	MsgBrd.set( text );

	PickModeInfo = piTorsion;
}

void PickedAtoms::print_pick_info(  ) {
	switch ( PickModeInfo ) {
	case piNone:
		if( numentries > 0 ) {
			print_atom_info(  );
		}
		break;
	case piDistance:
		print_distance(  );
		break;
	case piAngle:
		print_angle(  );
		break;
	case piTorsion:
		print_torsional(  );
	}
}

void PickedAtoms::highlight( char rgbtype ) {
	Atom *aptr = aaptr[0]->atomList + q[0];

	switch ( rgbtype ) {
	case 'r':
		aptr->color.r += DELTA_COLOR;
		break;
	case 'g':
		aptr->color.g += DELTA_COLOR;
		break;
	case 'b':
		aptr->color.b += DELTA_COLOR;
		break;
	default:
		aptr->color.r += DELTA_COLOR;
		aptr->color.g += DELTA_COLOR;
		aptr->color.b += DELTA_COLOR;
	}
	aptr->color.setGrey(  );
	// TOZLAB glutPostRedisplay();       
}

void PickedAtoms::darken( char rgbtype ) {
	Atom *aptr = aaptr[0]->atomList + q[0];

	switch ( rgbtype ) {
	case 'r':
		aptr->color.r -= DELTA_COLOR;
		break;
	case 'g':
		aptr->color.g -= DELTA_COLOR;
		break;
	case 'b':
		aptr->color.b -= DELTA_COLOR;
		break;
	default:
		aptr->color.r -= DELTA_COLOR;
		aptr->color.g -= DELTA_COLOR;
		aptr->color.b -= DELTA_COLOR;
	}
	aptr->color.setGrey(  );
	// TOZLAB glutPostRedisplay();       
}

void PickedAtoms::reverse(  ) {
	Atom *aptr = aaptr[0]->atomList + q[0];

	aptr->color.r ^= 0xFF;
	aptr->color.g ^= 0xFF;
	aptr->color.b ^= 0xFF;
	aptr->color.setGrey(  );
	//  TOZLAB glutPostRedisplay();      
}

void PickedAtoms::footer_print( int pickmode ) {
	MsgBrd.print( pickmode );
}

void PickedAtoms::footer_erase(  ) {
	MsgBrd.erase(  );
}

void PickedAtoms::display_pick( int windowSize, int colorFilter, GLubyte ballR, GLubyte ballG, GLubyte ballB ) {
	Atom *listi;
	int kk;

	// draw picked unitball
	glEnable( GL_LIGHTING );

	tRGBA c( ballR, ballG, ballB );
	glColor3ub( Color2R( colorFilter, &c ), Color2G( colorFilter, &c ), Color2B( colorFilter, &c ) );	// yellow ball

	if( aaptr[0] ) {
		glPushMatrix();
		Main::recalcmodelView( aaptr[0] );
	}

	for( kk = 0, listi = aptr[0]; ( kk < numentries ) && ( listi != NULL ); listi = aptr[++kk] ) {
		glTranslatef( listi->x, listi->y, listi->z );	// Move to the proper location
		glCallList( 1 );	// Render sphere display list
		glTranslatef( -listi->x, -listi->y, -listi->z );	// Move back
	}

	// draw pick lines
	glEnd(  );
	glDisable( GL_LIGHTING );
	glLineStipple( 1, 0x5555 );
	glEnable( GL_LINE_STIPPLE );
	glLineWidth( ( float ) windowSize / 2.0 );	//glLineWidth(3.0);
	glLoadName( NoSelectLoadName );	// no select for dash lines
	glBegin( GL_LINE_STRIP );
	for( kk = 0, listi = aptr[0]; kk < numentries && listi != NULL; listi = aptr[++kk] ) {
		glVertex3d( listi->x, listi->y, listi->z );
	}
	glEnd(  );
	glDisable( GL_LINE_STIPPLE );

	if( aaptr[0] ) {
		glPopMatrix();
	}
}
}
