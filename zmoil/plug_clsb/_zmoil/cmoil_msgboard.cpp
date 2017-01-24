//*************************************************************************************************
//*  Filename:   cmoil_msgboard.cpp
//*
//*  Description: 
//*    Display text messages at the bottom of the window
//*
//*  Modification History:
//*  
//*  Date           Developer   Description
//*  ------------   ----------- ---------------------------------------------------------------------
//*  Mar. 2001  Baohua Wang Initial Release
//*
//*************************************************************************************************


#include "zgltools.h"
#include "zglfont.h"
#include "zui.h"

#include "cmoil.h"
#include "cmoil_msgboard.h"

namespace CMOIL {
/* Our pick queue stores the last 4 atoms picked (0 is the most recent) */
MessageBoard::MessageBoard(  ) {
	text[0] = '\0';
	displayonce = 0;
} MessageBoard::MessageBoard( char *mytext ) {
	strncpy( text, mytext, LINE_MAXLEN );
	// TOZLAB glutPostRedisplay();        /* reflesh text buff */
}

void MessageBoard::set( char *mytext ) {
	if( displayonce == 0 ) {
		strcpy( text, mytext );
		//TOZLAB glutPostRedisplay();       /* reflesh text buff */
	}
}

void MessageBoard::append( char *mytext ) {
	strcat( text, mytext );
	// TOZLAB glutPostRedisplay();        /* reflesh text buff */
}

void MessageBoard::erase(  ) {
	text[0] = '\0';
	// TOZLAB glutPostRedisplay();     /* reflesh text buff */
}


// print text as footer
//
void MessageBoard::print( int inBlackBg, int lowerRight ) {
	if( text[0] == '\0' ) {
		return;				// nothing to print
	}

	glPushAttrib( GL_ALL_ATTRIB_BITS );
	glDisable( GL_LIGHTING );
	glDisable( GL_COLOR_MATERIAL );
	glDisable( GL_CULL_FACE );
	glDisable( GL_DEPTH_TEST );
	glDisable( GL_BLEND );
	glDisable( GL_TEXTURE_2D );
	glDisable( GL_FOG );
	glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

	glMatrixMode( GL_PROJECTION );
	glPushMatrix(  );
	glLoadIdentity(  );
	glMatrixMode( GL_MODELVIEW );
	glPushMatrix(  );
	::zglPixelMatrixFirstQuadrant(  );

	ZUI *z = ZUI::zuiFindByName( "pluginPanel" );
	int x = 5, y = 1;
	if( lowerRight ) {
		x = ( int ) ( z->w - 100 );
	}
	if( inBlackBg ) {
		// black quad to print against
		glBegin( GL_QUADS );
		glColor3f( 0.0, 0.0, 0.0 );
		glVertex3d( x - 5, 0.0, 0.0 );
		glVertex3d( lowerRight ? z->w : z->w - 100, 0.0, 0.0 );
		glVertex3d( lowerRight ? z->w : z->w - 100, 15, 0.0 );
		glVertex3d( x - 5, 15, 0.0 );
		glEnd(  );
	}

	glColor3f( 1.f, 1.f, 0.f );	// text color
	::zglFontPrint( text, (float)x, (float)y, "controls" );

	glMatrixMode( GL_PROJECTION );
	glPopMatrix(  );
	glMatrixMode( GL_MODELVIEW );
	glPopMatrix(  );

	glPopAttrib(  );
	return;
}

void MessageBoard::once(  ) {
	displayonce = 1;
}

}
