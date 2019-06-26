#ifndef __rd_kaylah__
#define __rd_kaylah__

#include "interp.h"
#include "input.h"


struct RD
{
	int id;
        int ntracked;
        int *tracked_id;
	
	double k_on;
	double k_off;
	double binding_radius;
	

	void get_tracked(surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, RD **tracked);
};

#endif
