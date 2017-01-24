#ifndef _CMOIL_CRD_H
#define _CMOIL_CRD_H
//*************************************************************************************************
//*  Filename:   cmoil_crd.h
//*
//*  Description: head file for Coordinate class
//*
//*  Modification History:
//*  
//*  Date        	Developer	Description
//*  ------------	-----------	---------------------------------------------------------------------
//*  Aug. 02 2005	Baohua Wang	Reshape
//*
//*************************************************************************************************

#include <math.h>
namespace CMOIL { 

class Atom;

class Coordinate
{
public:
  double x;
  double y;
  double z;

  Coordinate() {x=0; y=0; z=0;};
  Coordinate(double u, double v, double w) {x=u; y=v; z=w;};
  Coordinate(float *u) {x=*(u++); y=*(u++); z=*u;};
  Coordinate(Atom &a)  {set(a);}

  void set(double u, double v, double w) {x=u; y=v; z=w;};
  void set(float *u) {x=*(u++); y=*(u++); z=*u;};
  void set(Atom &a);

  const Coordinate &operator= (const Atom &b) ;
  const Coordinate &operator==(const Atom &b) ;
  double len() { return sqrt(x*x+y*y+z*z);  } ;    // length
  double mag2() { return x*x+y*y+z*z; }
  const Coordinate &norm();
  bool  inSphere(Coordinate &center, double radius) {return ((*this-center).len()<=1.005*radius);};

  friend Coordinate operator+ (const Coordinate &a, const Coordinate &b);
  friend Coordinate operator- (const Coordinate &a, const Coordinate &b);
  friend double     operator^ (const Coordinate &a, const Coordinate &b);  // dot
  friend Coordinate operator* (const Coordinate &a, const Coordinate &b);  // cross
  friend Coordinate operator* (const double &l, const Coordinate &a);      // scale
};

const Coordinate CrdZero;
const Coordinate ZAxis(0.,0.,1.);

inline double area(Coordinate &a, Coordinate &b, Coordinate &c)    {return 0.5*((a-c)*(b-c)).len(); }
inline double areaM2(Coordinate &a, Coordinate &b, Coordinate &c)  {return ((a-c)*(b-c)).len(); }
inline Coordinate normal(Coordinate oldCA, Coordinate oldCB, Coordinate newCA){return ((newCA-oldCA)*(newCA-oldCB)).norm();}
double     angleT(Coordinate c0, Coordinate c1, Coordinate c2, Coordinate c3);
Coordinate avg(Coordinate &a, Coordinate &b);

}

#endif
