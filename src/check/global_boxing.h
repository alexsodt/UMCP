#ifndef __global_boxingh__
#define __global_boxingh__
#include <stdlib.h>

struct box
{
	int npspace;
	int np;
	int *plist;
};
typedef struct global_boxing
{
	double PBC[3][3];
	int nx,ny,nz;
	box *boxes; // number of particles in this box.

	int addp( double *r, int pind );
	void clearBoxing( int *list= NULL, int nclear=0 );
	int getNearPts( double *r, int *plist, double rad_search );
} global_boxing;

void setup_global_boxing( double target_box_width, double PBC[3][3] );

#ifndef __global_boxingc__
extern int global_boxing_init;
extern global_boxing *boxing;
#endif

#endif
