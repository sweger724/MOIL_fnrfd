//*************************************************************************************************
//*  Filename:   cmoil_util.cpp
//*
//*  Description: 
//*    common used utility functions
//*
//*  Modification History:
//*  
//*  Date           Developer   Description
//*  ------------   ----------- ---------------------------------------------------------------------
//*  Aug. ?? 2000   ??      Initial Development
//*  July 2001  Baohua Wang Extracted from cmoil.cpp
//*  Sept 2007  Thomas Blom Update for Linux/OSX
//*************************************************************************************************

#include "wingl.h"
#include "GL/gl.h"
#include "GL/glu.h"

#ifdef __APPLE__
#include "sys/wait.h"			// for id_t
#include "i386/types.h"			// for u_int32_t
#include "sys/stat.h"			// for dev_t
#include "sys/uio.h"			// for ssize_t
#endif

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "cmoil_crd.h"
#include "cmoil_atom.h"
#include "cmoil_camera.h"

extern float Zmoil_eyeZ;

 // debug: let's look at the modelview matrix
void printModelViewMatrix( char *msg )
{
	double m[16];
	int i;
	FILE *f = fopen( "modelview.txt", "w" );
	glGetDoublev( GL_MODELVIEW_MATRIX, m );
	if( msg ) {
		fprintf( f, "%s\n", msg );
	}
	fprintf( f, "ModeView:\n" );
	for( i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			fprintf( f, "\t%g", m[i * 4 + j] );
		}
		fprintf( f, "\n" );
	}
	glGetDoublev( GL_PROJECTION_MATRIX, m );
	fprintf( f, "\nProjection:\n" );
	for( i = 0; i < 4; i++ ) {
		for( int j = 0; j < 4; j++ ) {
			fprintf( f, "\t%g", m[i * 4 + j] );
		}
		fprintf( f, "\n" );
	}
	fclose( f );
}


namespace CMOIL {
void Coordinate::set( Atom & a ) {
	x = a.x;
	y = a.y;
	z = a.z;
};

const Coordinate & Coordinate::norm(  ) {
	double l = sqrt( x * x + y * y + z * z );
	if( l != 0. ) {
		l = 1 / l;
		x *= l;
		y *= l;
		z *= l;
	}
	return *this;
}

// copy dynamic coordinates for graphic display
const Coordinate & Coordinate::operator=( const Atom & right ) {
	x = right.x;
	y = right.y;
	z = right.z;
	return *this;
}

// copy original coordinate
const Coordinate & Coordinate::operator==( const Atom & right ) {
	x = right.orgX;
	y = right.orgY;
	z = right.orgZ;
	return *this;
}

Coordinate operator-( const Coordinate & a, const Coordinate & b ) {
	Coordinate answer;
	answer.x = a.x - b.x;
	answer.y = a.y - b.y;
	answer.z = a.z - b.z;
	return answer;
}

Coordinate operator+( const Coordinate & a, const Coordinate & b ) {
	Coordinate answer;
	answer.x = a.x + b.x;
	answer.y = a.y + b.y;
	answer.z = a.z + b.z;
	return answer;
}


// dot
double operator^( const Coordinate & a, const Coordinate & b ) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}

//cross
Coordinate operator*( const Coordinate & a, const Coordinate & b ) {
	Coordinate answer;
	answer.x = a.y * b.z - a.z * b.y;
	answer.y = a.z * b.x - a.x * b.z;
	answer.z = a.x * b.y - a.y * b.x;
	return answer;
}

//scale
Coordinate operator*( const double &l, const Coordinate & a ) {
	Coordinate answer;
	answer.x = l * a.x;
	answer.y = l * a.y;
	answer.z = l * a.z;
	return answer;
}

Coordinate avg( Coordinate & a, Coordinate & b ) {
	Coordinate answer;
	answer.x = ( a.x + b.x ) / 2.0;
	answer.y = ( a.y + b.y ) / 2.0;
	answer.z = ( a.z + b.z ) / 2.0;
	return answer;
}

// compute torsional angle of (c0, c1, c2, c3)
double angleT( Coordinate c0, Coordinate c1, Coordinate c2, Coordinate c3 ) {
	Coordinate e1 = c0 - c1, e2 = c1 - c2, e3 = c2 - c3;
	double ang;
	ang = acos( ( ( e1 * e2 ).norm(  ) ) ^ ( ( e2 * e3 ).norm(  ) ) ) * 180 / PI;
	if( ( ( e2 * e1 ) ^ ( c0 - c3 ) ) < 0 )
		ang = -ang;
	return ang;
}

// compute torsional angle of (c0, c1, c2, c3)
double angleT( Atom & a0, Atom & a1, Atom & a2, Atom & a3 ) {
	Coordinate c0( a0 ), c1( a1 ), c2( a2 ), c3( a3 );
	return angleT( c0, c1, c2, c3 );
}

#define M2(a1, b1, a2, b2)        ( (a1)*(b2)-(a2)*(b1) )
#define M3( a1, b1, c1, a2, b2, c2, a3, b3, c3) \
( (a1)*M2(b2,c2,b3,c3)-(b1)*M2(a2,c2,a3,c3)+(c1)*M2(a2,b2,a3,b3) )
//
// Ax+By+Cz=D, return [A,B,C,D] as Coordinate type
// A(x-x0)+B(y-y0)+C(z-z0)=0,  D=x0A+y0B+z0C
//
// one/two point:    A=0, B=0
// three points: general form
// notes:  normal's always at +Z.
//
float *pointsPlane( int nPoints, Coordinate * c, char cutoffAxe ) {
	static float coeff[4];
	Coordinate n( 0., 0., 0. );	// normall

	for( int i = 0; i < 4; i++ )
		coeff[i] = 0.0;		// A=B=C=D=0          init  

	switch ( nPoints ) {
	case 1:
	case 2:
		if( cutoffAxe == 'X' )
			n.x = 1.;
		else if( cutoffAxe == 'Y' )
			n.y = 1.;
		else
			n.z = 1.;
		break;
	default:
		n = normal( c[0], c[1], c[2] );
		if( n.z < 0 )
			n = ( -1 ) * n;	// normal points to +Z
	}
	coeff[0] = ( float ) n.x;
	coeff[1] = ( float ) n.y;
	coeff[2] = ( float ) n.z;
	coeff[3] = ( float ) ( n ^ c[0] );

	return coeff;
}

// m[16] :     transformation matrix
// plane[4]:   plane equation (0,0,1,-Z)
// return:     newplane[4]
double *ZplaneTransform( double *m, double z, double *newplane ) {
	if( newplane == NULL )
		return NULL;

	// input plane normal  (0,0,1), cross point (0,0,Z)
	Coordinate n( m[2] * z, m[6] * z, m[10] * z );
	n.norm(  );
	Coordinate p( m[2] * z + m[3], m[6] * z + m[7], m[10] * z + m[11] );
	double d = -( n ^ p );
	*newplane = n.x;
	*( newplane + 1 ) = n.y;
	*( newplane + 2 ) = n.z;
	*( newplane + 3 ) = d;
	return newplane;
}

// for cutoff edge display
//

// get a line definition from input two points
// return A1x+B1y+C1z=D1 & A2x+B2y+C2z=D2 as 8 floats
float *pointsLine( Coordinate & pt1, Coordinate & pt2 ) {
	static float line[8];
	line[0] = ( float ) ( pt2.y - pt1.y );	// A1
	line[1] = ( float ) ( pt1.x - pt2.x );	// B1
	line[2] = 0.;			// C1
	line[3] = ( float ) ( line[0] * pt1.x + line[1] * pt1.y );	// D1

	line[4] = 0.;			// A2
	line[5] = ( float ) ( pt2.z - pt1.z );	// B2
	line[6] = ( float ) ( pt1.y - pt2.y );	// B
	line[7] = ( float ) ( line[5] * pt1.y + line[6] * pt1.z );	// D
	return line;
}

// Use Gaussian elimination to calculate the solution to a 3D linear system, AX=B.
// argMatrix = |eq1(a11 a12 a13)|
//             |eq2(a21 a22 a23)| == [A|B]
//             |eq3(a31 a32 a33)|     
// return true if a sole solution exists.  Otherwise return false.
// X is the solution [x y z] when return value is true
// eqI = { aI1, aI2, aI3, bI};
bool lequation3( float *eq1, float *eq2, float *eq3, float *X ) {
	const int N = 3, M = 4;	// 3D only
	int i, j, k, m, ii;
	float A[N * M], tmp[M], a;
	int rowbytes = sizeof( float ) * 4;

	memcpy( A, eq1, rowbytes );
	memcpy( &A[M], eq2, rowbytes );
	memcpy( &A[M + M], eq3, rowbytes );

	try {
		// get upright matrix
		for( k = 0; k < N; k++ )	// loop on all rows
		{
			// k :  should also be the first non-zero element of the ith row
			// i :  loop after k element
			// m :  memory position of ith element
			//
			ii = -1;		// the row with the kth element != 0
			for( i = k; i < N; i++ )	// loop start at kth row.  m: ith row kth element 
			{
				m = i * M + k;
				if( A[m] != 0 )	// set kth' coef=1
				{
					a = 1.f / A[m];
					for( j = 1; j <= N - k; j++ )
						A[m + j] *= a;
					A[m] = 1.;

					if( ii < 0 )
						ii = i;	// the row with the kth' coef != 0
				}
			}
			if( ii == -1 )
				return false;	// no solution
			else if( ii != k )	// swap
			{
				memcpy( tmp, A + ii * M, rowbytes );
				memcpy( A + ii * M, A + k * M, rowbytes );
				memcpy( A + k * M, tmp, rowbytes );
			}
			//now, k row's kth element is not zero. do elimination to get upright matrix
			ii = k * ( N + 2 );	// kth element of kth row , 
			for( i = k + 1; i < N; i++ )	// loop from k+1 row
			{
				m = i * M + k;
				if( A[m] != 0 )	// A[m]: kth element in ith row
				{
					for( j = 0; j <= N - k; j++ )
						A[m + j] -= A[ii + j];	// kth coeff are all 1 now
				}
			}
		}

		// get [1] matrix
		for( k = N - 1; k >= 0; k-- )	// row indexed as [0,N-1],  A[ii]=1 for now
		{
			ii = k * ( N + 2 );
			for( i = 0; i < k; i++ ) {
				m = i * M + k;	// kth element in ith row
				a = A[m];
				if( A[m] != 0 ) {
					for( j = 0; j <= N - k; j++ ) {
						A[m + j] -= A[ii + j] * a;
					}
				}
			}
		}
		// result in last column
		for( i = 0; i < N; i++ )
			*( X + i ) = A[N + i * M];
	}
	catch( char *str )		// (char *str) 
	{
		str = NULL;
		return false;		// any overflowe indicates a unsolvable case
	}
	return true;
}							// lequation3()



// AX=B  equations for plane
// plane = { a1, a2, a3, b};
bool linePlaneCrossPoint( float *plane, Coordinate & pt1, Coordinate & pt2, float *crossPt ) {
	float *line;
	if( ( line = pointsLine( pt1, pt2 ) ) != NULL )
		return lequation3( plane, line, line + 4, crossPt );
	else
		return false;
}

//**********************************************
//*   Camera clall
//**********************************************
Camera::Camera(  ) {
	focallength = 0;
	clipFront = 0;
	clipBack = 0;
}

void Camera::SetSizes( int wWidth, int wHeight, double maxZ, double minZ ) {
	ratio = wWidth / ( float ) wHeight;

	double aperture = 45;	/* Camera aperture : 45 degree  */
	radians = aperture * AngleToRadius * 0.5;	// half of aperture for each eye from middle

	maxz = maxZ;
	minz = minZ;
	nearz = 30.;
	farz = nearz + fabs( maxZ - minZ ) + 3;	// 10000;
	SetVp( 0, 0, Zmoil_eyeZ );
}

// input is view/camera position
void Camera::SetVp( double x, double y, double z )	// focal point is origin
{
	if( z == 0 )			// compute camera position
		z = fabs( maxz ) + nearz;
	focallength = z;		// 70;
	eyesep = 0.5 * focallength / 20;	// from middle  = eye_separation/2

	wd2 = nearz * tan( radians );
	ndfl = nearz / focallength;

	pr.set( 0, 0, 0 );		// focus to origin
	vu.set( 0, 1, 0 );		// up is Y axis

	vp.set( x, y, z );
	vd = ( -1 ) * vp;
	r = ( vd * vu ).norm(  );
	r = eyesep * r;
}

//
//********** end of Camera ************************************


// simulate ModelViews for left or right eye  ( -1|1|0 )
//
void Camera::LookAt( int leftRight, int stereomode ) {
	//SetVp();
	eye = leftRight;		// left, middle, right

	glMatrixMode( GL_MODELVIEW );
	//glPopMatrix(  );
	//glLoadIdentity(  );

	// right buffer used only with HARDWARE_STEREO
	if( eye == RIGHT && stereomode == HARDWARE_STEREO )
		glDrawBuffer( GL_BACK_RIGHT );
	else
		glDrawBuffer( GL_BACK_LEFT );

	Coordinate p = vp + eye * r;	// eye position 
	Coordinate q = p + vd;	// center
	gluLookAt( p.x, p.y, p.z, q.x, q.y, q.z, vu.x, vu.y, vu.z );	// up
	//glPushMatrix(  );
}

void Camera::Clip( int front, int back ) {
	if( front == ResetClip || back == ResetClip ) {
		clipFront = 0;
		clipBack = 0;
	}
	else if( front >= MinusClipDelta && front <= PlusClipDelta ) {
		clipFront += ( front - SameClipDelta );
		clipBack += ( back - SameClipDelta );
	}
}

// check/chang path for valid dir separator, limited to FILENAME_MAXLEN
void SetValidPath( char *path ) {
	char *cptr = path;
	while( ( cptr = strchr( cptr, DIRDELIMalt ) ) != NULL && ( cptr - path ) < FILENAME_MAXLEN ) {
		*cptr = DIRDELIM;
	}
}


}
