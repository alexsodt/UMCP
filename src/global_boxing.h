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

	void setPBC( double PBC_vec[3][3], double *alphas );
	int addp( double *r, int pind );
	void clearBoxing( int *list= NULL, int nclear=0 );
	int getNearPts( double *r, int *plist, double rad_search );
	void setup_boxing( double target_box_width, double PBC[3][3] );
} global_boxing;



#endif
