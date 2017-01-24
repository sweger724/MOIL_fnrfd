//*************************************************************************************************
//*  Filename:   cmoil_msgboard.cpp
//*
//*  Description: 
//*    Display text messages at the bottom of the window
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Mar. 2001	Baohua Wang	Initial Release
//*
//*************************************************************************************************
#include "cmoil.h"
#include "cmoil_msgboard.h"

namespace CMOIL {
/* Our pick queue stores the last 4 atoms picked (0 is the most recent) */
MessageBoard::MessageBoard() 
{
   text[0] = '\0';
   displayonce=0;
}

MessageBoard::MessageBoard(char *mytext) 
{
  strncpy(text, mytext, LINE_MAXLEN);
  glutPostRedisplay();		/* reflesh text buff */
}

void MessageBoard::set(char *mytext) 
{
  if (displayonce==0)
  {
    strcpy(text, mytext);
    glutPostRedisplay();		/* reflesh text buff */
  }
}

void MessageBoard::append(char *mytext) 
{
  strcat(text, mytext);
  glutPostRedisplay();		/* reflesh text buff */
}

void MessageBoard::erase() 
{
  text[0]='\0';
  glutPostRedisplay();     /* reflesh text buff */
}

// print text as footer
//
void MessageBoard::print(int inBlackBg, int lowerRight ) 
{   
   if ( text[0] == '\0' )
      return;			// nothing to print

   /* Push current matrix mode and viewport attributes */
   glPushAttrib(GL_TRANSFORM_BIT | GL_VIEWPORT_BIT);
 
   glMatrixMode( GL_MODELVIEW );
   glPushMatrix();
   glLoadIdentity();

   glMatrixMode( GL_PROJECTION );
   glPushMatrix();
   glLoadIdentity();
   gluOrtho2D( 0.0, 1.0, 0.0, 1.0 );
 
   if (inBlackBg)		// patch a black area in the bottom of the window
   {
     glBegin(GL_QUADS);
     glColor3f ( 0.0, 0.0, 0.0);
     glVertex3d( 0.0, 0.0, 0.0);
     glVertex3d( 1.0, 0.0, 0.0);
     glVertex3d( 1.0, 0.035, 0.0);
     glVertex3d( 0.0, 0.035, 0.0);
     glEnd();
   }

   glColor3f ( 2.0, 2.0, 0.0 );		      // text color
   float xpos = lowerRight ? .80 : .01;
   glRasterPos4f( xpos, 0.01, 1.0, 1.0);
   for (unsigned int i=0; i < strlen(text); i++ )
     glutBitmapCharacter( GLUT_BITMAP_HELVETICA_12, text[i]);

   glPopMatrix();
   glMatrixMode(GL_MODELVIEW); 
   glPopMatrix(); 

   glPopAttrib();
}

void MessageBoard::once() 
{
  displayonce=1;
}

}
