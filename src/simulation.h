#ifndef __simulationh__
#define __simulationh__

typedef struct surface_record
{
	int id; // for locating bound particles
	surface *theSurface;
	FILE *tFile; // for writing trajectory
	struct surface_record *next;
	double *r;
	double *g;
	double *pp;
	double *qdot;
	double *qdot0;

	// for harmonic modes
	int NQ;
	double *gen_transform;
	double *scaling_factor;
	int do_gen_q;
	int doing_spherical_harmonics;
	int doing_planar_harmonics;
	double *output_qvals;
	

	SparseMatrix *EFFM;
	SparseMatrix *MMat;
	force_set *theForceSet;
} surface_record;

typedef struct
{
	surface_record *allSurfaces;

	pcomplex **allComplexes;	
} Simulation;


#endif
