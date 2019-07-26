#ifndef __rd_kaylah__
#define __rd_kaylah__

#include "interp.h"
#include "input.h"
#include "fpr_subroutines/2D.h"

typedef struct
{
	int id;
	int sid1; // sub id 1 (of mine)
	int sid2; // sub id 2 (of the complex with int id; above)
	double prev_norm;
	double prev_sep;
	double curr_sep;
	double curr_norm;
	int info;
} RD_tracked_info;

typedef struct
{
        int ntracked;
        int ntracked_space;
        RD_tracked_info *tracked_info;
	RD_tracked_info *tracked_new;
} RD_tracked;

struct RD
{
	RD_tracked **tracked;	
	
	double dt;
	double Dtot;
	double k_on;
	double k_off;
	double binding_radius;
	double Rmax;
	
	void init(Simulation *theSimulation, double time_step);
	void get_tracked( Simulation *theSimulation );
	void do_rd( Simulation *theSimulation );
};

#endif
