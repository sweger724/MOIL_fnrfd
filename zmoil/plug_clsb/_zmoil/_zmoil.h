#ifndef _ZMOIL_H
#define _ZMOIL_H


#include "ztime.h"
#include "tbutil.h"

extern int gEnableFog;
extern float Zmoil_BondLenMax;
extern int Zmoil_VelocityDisplay;
extern int Zmoil_PeptidePlaneAtoms;
extern int Zmoil_DisplayHBondsWater;
extern int Zmoil_DisplayHBonds;
extern int Zmoil_HBondThickness;
extern float Zmoil_HBondDistance;
extern int Zmoil_DisplayWater;
extern float Zmoil_SpaceFillingScale;
extern float Zmoil_CavityMinSurfArea;
extern int Zmoil_gridColor;
extern double Zmoil_gridColorMax;
extern double Zmoil_gridColorMin;


extern void tempMessage( char *msg, char color='w', ZTime displaySeconds=3.0 );
extern void updateIndexSlider( bool labelOnly );


extern ZTime globalTime;
extern int globalTimerOn;


#endif
