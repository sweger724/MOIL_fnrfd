//*************************************************************************************************
//*  Filename:   cmoil_pick.cpp
//*
//*  Description: 
//*    Methods of PickedAtoms class for disply information of selected atoms
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. 2000	??		Initial Development
//*  Mar. 2001    Baohua Wang Seperated from cmoil.cpp and modified for message display
//*  July 2001	Baohua Wang	Reshape to support mutiple structure display
//*  Sept 2007		Thomas Blom		Update for Linux/OSX
//*************************************************************************************************
#ifdef __APPLE__
	#include "sys/uio.h"  // for ssize_t
	#include "pthread.h"
#endif


#include <stdlib.h>
#include "cmoil_pick.h"
#include "cmoil.h"
#include "cmoil_globals.h"

 
namespace CMOIL {
#define IsOverlap    (image2 >= 0)
const int UnitballLoadname1st=1000000;

PickedAtoms::PickedAtoms()
{
  currImage = 0;
  reset();
}

void PickedAtoms::reset(int imageIndex1, int imageIndex2)	// should pass in two overlaped image index also
{
  char msg[80];
  image1 = imageIndex1;
  image2 = imageIndex2;

  strcpy (msg , "In PickMode:  please pick atoms");

  if (imageIndex1 >= 0 || imageIndex2 >= 0  )
  {
    if ( AAptrs[imageIndex1]->InParm.filetag == fPTH ||
         AAptrs[imageIndex1]->InParm.filetag == fDCD )
    {
      sprintf ( msg+strlen(msg), " struct #%d", AAptrs[imageIndex1]->currentmol);
    }
    MsgBrd.set(msg);
  }

  if ( image1 != imageIndex1 || image2 != imageIndex2 )
  {
    numentries=0;					// reset for changed image only
    image1 = imageIndex1;
    image2 = imageIndex2;
  }
}

int PickedAtoms::isPicked(Atom *atom)
{
  for (int i=0; i<numentries; i++)
  {
    if (aptr[i]==atom)
         return 1;
  }
  return 0;
}

// set at most 3 recent picked Atom Crds, return number of picked crds
int PickedAtoms::getPickedCrds(Coordinate *crds)
{
  int numPicked = (numentries>3)?3:numentries;
  for (int i=0; i<numPicked ; i++ )
    crds[i] = *aptr[i];
  return numPicked;
}

// mouse clicked at (x,y), return the select atom number.
// Since GL_SELECT render mode is not supported by hardware acceleration,
// implement the select in GL_RENDER mode.
//
int PickedAtoms::get_atomnum(int x, int y) 
{
  GLint viewport[4];
  GLdouble pM[16], mM[16];

  glGetIntegerv(GL_VIEWPORT, viewport);
  glGetDoublev(GL_PROJECTION_MATRIX, pM);
  glGetDoublev(GL_MODELVIEW_MATRIX, mM);
  
  Atom *a=NULL;  // , *lastA=NULL;
  int atomIndex=-1, structN=1, atomCnt=0;
  for (structN=1; structN<=2; structN++)
  {
    if ( structN==1 && image1 >=0  )
    {
        a=AAptrs[image1]->atomList;
        atomCnt= AAptrs[image1]->numAtom;
    } 
    else if ( IsOverlap)
    {
        a=AAptrs[image2]->atomList;
        atomCnt= AAptrs[image2]->numAtom;
    } 
    else
    {
      atomCnt=0;
    }
 
    y=viewport[3]-(GLint)y -1;
    double minz=0., ox, oy, oz, wx, wy, wz ;
    double minR2=900, R2;

    for (int i=0; i<atomCnt; i++, a++)
    {
       if ( a->skip==0 && gluProject(a->x, a->y, a->z,
                  mM, pM, viewport, &wx, &wy, &wz) )
       {
         if (wz>=0)
         {
           wx += sqrt((wx-x)*(wx-x)+(wy-y)*(wy-y));  // get distance in window space
           gluUnProject(wx, wy, wz, mM, pM, viewport, &ox, &oy, &oz);  // map to object space
           R2 = sqrt((a->x-ox)*(a->x-ox)+(a->y-oy)*(a->y-oy)+(a->z-oz)*(a->z-oz))/a->radius; // compare tdistance to radius
           if ( R2 < minR2 || (R2==minR2 && wz<minz) )
           {
              minR2 =R2;
              atomIndex=i;
              minz=wz;
           }   
         }         
       }
    }
  }

  if (structN == 2 && atomIndex >= 0 )
    atomIndex += AAptrs[image1]->numAtom; 

  return atomIndex;
}

// put selected atom in to picked list
void PickedAtoms::put(int x, int y) {
  char text[LINE_MAXLEN];	// temp for Footer display
  int i, a;	

  a = get_atomnum(x,y);		// get selected sequece number from buffer 
  
  //printf("pick #=%d, %d, %d\n", a, AAptrs[image1]->numAtom, AAptrs[image2]->numAtom);

  if ( a < 0 ) // || a >= UnitballLoadname1st+ATOM_MAX_PICK-1 )
    return;				// nothing selected
                 
  for(int j=numentries; j>=1; j--)       // shift circular array
  {
    q[j] = q[j-1];
    aaptr[j] = aaptr[j-1];
    aptr[j] = aptr[j-1];
  }

  if ( a < UnitballLoadname1st )          // not already selected unit ball
  {
    i = image1;
    if ( a >= AAptrs[i]->numAtom )
    {
      a -= AAptrs[i]->numAtom;           // overlap case
      i = image2;
    }
    if (i<0 || a>=AAptrs[i]->numAtom )
      return;

    if ( numentries>=ATOM_MAX_PICK ) 
      aptr[ATOM_MAX_PICK]->dispmode &= (~UNITBALL);      // unset the flag for oldPicked;

    q[0] = a;
    aaptr[0] = AAptrs[i];
    aptr[0]  = aaptr[0]->atomList+a;
    aptr[0]->dispmode |= UNITBALL; 
  }
  else
  { 
    a=a-UnitballLoadname1st;
    if (a >0 && a<ATOM_MAX_PICK)         
    {
      a++;
      q[0]    = q[a];         // copy over selected
      aaptr[0]= aaptr[a];
      aptr[0] = aptr[a];
    }
  }

  if (numentries<ATOM_MAX_PICK )
      numentries++;
 
  // print coordinates only for XYZ type crd file
  //
  
  if (aaptr[0]->InParm.filetag == fXYZ )
    sprintf(text, "(%.5g, %.5g, %.5g)", aptr[0]->orgX, aptr[0]->orgY, aptr[0]->orgZ );
  else
    sprintf(text, "%s (%.5g, %.5g, %.5g) Radius=%.5gA", 
             atom_name(0, (aaptr[0]->InParm.filetag)) ,       
             aptr[0]->orgX, aptr[0]->orgY, aptr[0]->orgZ, aptr[0]->radius);

  MsgBrd.set(text);
  //display_pick();

}	// PickedAtoms::put()

// get full atom name for display, aaptr[atomIndex] -- error?
char  *PickedAtoms::atom_name ( int atomIndex,  char fileTag ) 
{
   static char name[160];
  
   if ( atomIndex >3 || atomIndex < 0 )			// support 0-3 index at same time
   {
     name[0] ='\0';
     return name;
   }

   char   *myname= name + atomIndex*40;
   if ( fileTag == fPDB )
   {
      sprintf(myname, "%s:%.0d:%c:%d:%s:%s",
          IsOverlap?aaptr[atomIndex]->moleculeName:"",  
          (aaptr[atomIndex]->chains+aptr[atomIndex]->chainIndex)->model,  // display model, chain for PDB file
          (aaptr[atomIndex]->chains+aptr[atomIndex]->chainIndex)->id,  
          aaptr[atomIndex]->resNum(aptr[atomIndex]),  
          aaptr[atomIndex]->resName[aptr[atomIndex]->rn], 
          aptr[atomIndex]->atomName  );
   }
   else if ( fileTag == fPTH ||fileTag == fDCD )   // binary crd to display structure #
   {
      sprintf(myname, "%s:%.0d:%d:%s:%s",
          IsOverlap?aaptr[atomIndex]->moleculeName:"",  
          aaptr[atomIndex]->currentmol,         // display currentmol for pth or dcd file 
          aaptr[atomIndex]->resNum(aptr[atomIndex]),  
          aaptr[atomIndex]->resName[aptr[atomIndex]->rn], 
          aptr[atomIndex]->atomName  );
   }
   else
   {
      sprintf(myname, "%s:%d:%s:%s",
          IsOverlap?aaptr[atomIndex]->moleculeName:"", 
          aaptr[atomIndex]->resNum(aptr[atomIndex]), 
          aaptr[atomIndex]->resName[aptr[atomIndex]->rn], 
          aptr[atomIndex]->atomName);
   }
   return myname;
}   // atom_name()

// print the distance btwn first two pick entries 
void PickedAtoms::print_distance() 
{
  Coordinate t, t1, t2;
  double d;
  char text[LINE_MAXLEN];

  if (numentries <2) {
    MsgBrd.set("Please select two atoms for a distance\n");
    return;
  } 
  
  t1 == *aptr[0];
  t2 == *aptr[1]; 
  t = t1 -  t2;
  d = t.len();

  char hasPDB = aaptr[0]->InParm.filetag;
  if ( aaptr[0]->InParm.filetag == fPDB || aaptr[1]->InParm.filetag == fPDB )
    hasPDB = 'b';
  else
    hasPDB = aaptr[0]->InParm.filetag ;                   // set to the crd filetag
  sprintf(text, "Distance of (%s, %s): %.5g A",  atom_name(0,hasPDB), atom_name(1,hasPDB), d);

  MsgBrd.set(text);

}	// print_distance()

void PickedAtoms::print_angle() 
{
  Coordinate t1,t2;
  Coordinate c0,c1,c2;
  double     ac;
  char       text[LINE_MAXLEN];
 
  if (numentries <3) 
  {
    MsgBrd.set("Please select three atoms for an angle");
    return;
  } 

  c0 == *aptr[0];   
  c1 == *aptr[1];   
  c2 == *aptr[2];   
  t1 = (c0-c1).norm();
  t2 = (c2-c1).norm();
  ac = acos(t1^t2)*180/PI;     //dot

  char hasPDB ;
  if ( aaptr[0]->InParm.filetag == fPDB ||  aaptr[1]->InParm.filetag == fPDB || 
                 aaptr[2]->InParm.filetag == fPDB )
    hasPDB = fPDB; 
  else
    hasPDB = aaptr[0]->InParm.filetag ;
  
  sprintf(text, "Angle of (%s, %s, %s): %.5g Deg.", 
      atom_name(0,hasPDB), atom_name(1,hasPDB), atom_name(2,hasPDB),  ac );

  MsgBrd.set(text);		
}

void PickedAtoms::print_torsional() 
{
  char       text[LINE_MAXLEN];

  if (numentries <4) 
  {
    MsgBrd.set("Please select four atoms for a torsinal angle");
    return;
  } 

  char hasPDB ;
  if ( aaptr[0]->InParm.filetag == fPDB ||  aaptr[1]->InParm.filetag == fPDB || 
                 aaptr[2]->InParm.filetag == fPDB ||  aaptr[3]->InParm.filetag == fPDB )
    hasPDB = fPDB;
  else 
    hasPDB = aaptr[0]->InParm.filetag ;
  
  sprintf(text, "Tors.Angle of (%s, %s, %s, %s): %.5g Deg.", 
         atom_name(0,hasPDB), atom_name(1,hasPDB), 
         atom_name(2,hasPDB), atom_name(3,hasPDB), 
         angleT(*aptr[0], *aptr[1], *aptr[2], *aptr[3])  );

  MsgBrd.set(text);
}

void PickedAtoms::highlight(char rgbtype) 
{
   Atom *aptr = aaptr[0]->atomList+q[0];

   switch (rgbtype)
   {
      case 'r':
        aptr->color.r += DELTA_COLOR ;
        break;
      case 'g':
        aptr->color.g += DELTA_COLOR ;
        break;
      case 'b':
        aptr->color.b += DELTA_COLOR ;
        break;
      default:  
        aptr->color.r += DELTA_COLOR ;
        aptr->color.g += DELTA_COLOR;
        aptr->color.b += DELTA_COLOR;
   }
   aptr->color.setGrey();
   glutPostRedisplay();		
}

void PickedAtoms::darken(char rgbtype) 
{
   Atom *aptr = aaptr[0]->atomList+q[0];

   switch (rgbtype)
   {
      case 'r':
        aptr->color.r -= DELTA_COLOR ;
        break;
      case 'g':
        aptr->color.g -= DELTA_COLOR ;
        break;
      case 'b':
        aptr->color.b -= DELTA_COLOR ;
        break;
      default:  
        aptr->color.r -= DELTA_COLOR ;
        aptr->color.g -= DELTA_COLOR;
        aptr->color.b -= DELTA_COLOR;
   }
   aptr->color.setGrey();
   glutPostRedisplay();		
}

void PickedAtoms::reverse() 
{
   Atom *aptr = aaptr[0]->atomList+q[0];

   aptr->color.r ^= 0xFF;
   aptr->color.g ^= 0xFF;
   aptr->color.b ^= 0xFF;
   aptr->color.setGrey();
   glutPostRedisplay();		
}

void PickedAtoms::footer_print(int pickmode) 
{
  MsgBrd.print(pickmode);
}

void PickedAtoms::footer_erase() 
{
  MsgBrd.erase();
}

void PickedAtoms::display_pick(int windowSize, int colorFilter, GLubyte ballR, GLubyte ballG, GLubyte ballB)
{
    Atom *listi;
    int kk;

    // draw picked unitball
    glEnable(GL_LIGHTING);

    tRGBA c(ballR, ballG, ballB);
    glColor3ub(Color2R(colorFilter,&c),Color2G(colorFilter,&c),Color2B(colorFilter,&c));    // yellow ball

    for (kk=0, listi=aptr[0]; (kk<numentries)&&(listi!=NULL); listi=aptr[++kk]) 
    {
      glTranslatef(listi->x, listi->y, listi->z);	// Move to the proper location
      glCallList(1); 					      // Render sphere display list
      glTranslatef(-listi->x, -listi->y, -listi->z);  // Move back
    }

    // draw pick lines
    glEnd();
    glDisable(GL_LIGHTING);
    glLineStipple(1, 0x5555);
    glEnable(GL_LINE_STIPPLE);
    glLineWidth((float)windowSize/2.0);                //glLineWidth(3.0);
    glLoadName(NoSelectLoadName);    // no select for dash lines
    glBegin(GL_LINE_STRIP);
    for (kk=0, listi=aptr[0]; kk<numentries && listi!=NULL ; listi=aptr[++kk]) 
      glVertex3d(listi->x, listi->y, listi->z);
    glEnd();
    glDisable(GL_LINE_STIPPLE);
}
}

