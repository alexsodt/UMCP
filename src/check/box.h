#ifndef __boxh__
#define __boxh__

#define BOX_SPHERE	0
#define BOX_CYLINDER	1
#define BOX_PLANE	2
#define BOX_3D		3
#define BOX_HEX		4

double normalize( double *dr );
void initializeBox(  int type, int nparticles, double R, double Lx, double Ly, double Lz, double cutoff );
void updateBox( int particle, double *r );
void removeBox( int particle );
int getNear( int particle, int *neighbors );
int getNear( double *r, int *neighbors );
void place( int particle, double *r );
void report( int p1, int p2 );
int getNear( double *r, int *neighbors, int *del );

#endif


