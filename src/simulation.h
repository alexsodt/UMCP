#ifndef __simulationh__
#define __simulationh__

#ifndef __interph__
struct surface;
#endif

#ifndef __sparsematrixh__
struct SparseMatrix;
#endif

#ifndef __lagrangeh__
struct force_set;
#endif

#ifndef __pcomplexh__
struct pcomplex;
#endif

#ifndef __inputh__
struct parameterBlock;
#endif

#include <stdio.h>

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
	int nc;

	// for minimization
	int temp_min_offset;
	double *temp_r;
	double *temp_g;

	// for harmonic modes
	int NQ;
	double *gen_transform;
	double *scaling_factor;
	int do_gen_q;
	int doing_spherical_harmonics;
	int doing_planar_harmonics;
	double *output_qvals;

	double *next_pp;
	double *del_pp;
	int n_real_q;
	double *nav_Q;
	double *av_Q;
	double *av_Q2;
	double *av_Q_T;
	double *nav_Q_T;
	double *QV;
	double *Qdot;
	double *Qdot0;
	double *Qdot0_trial;

	double *qdot_temp;

	SparseMatrix *EFFM;
	SparseMatrix *MMat;
	force_set *theForceSet;
} surface_record;

typedef struct Simulation
{
	surface_record *allSurfaces;

	pcomplex **allComplexes;	
	int ncomplex;	
	double alpha[3];
	double PBC_vec[3][3];

// BEGIN TRIVIAL UTILITY FUNCTIONS

	void wrapPBC( double *dr, double *alphas );

	void fetch( int sid, surface **theSurface, double **r )
	{
		*theSurface = NULL;
		*r = NULL;
		for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		{
			if( sRec->id == sid )
			{
				*theSurface = sRec->theSurface;
				*r = sRec->r;
				return;
			}
		}
	}
	struct surface_record * fetch( int sid )
	{
		for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		{
			if( sRec->id == sid )
			{
				return sRec;
			}
		}
		return NULL;
	}

// MEMBER FUNCTIONS

	void sample_B_hist( double *B_hist, double *A2dz2_sampled,
			int sample_type, 
			int nsamples, double maxr, int nbins, int shape_correction );
	void area_MC_move( double beta, int move_type, parameterBlock *block );
	void writeLimitingSurfacePSF(FILE *theFile );
	void writeLimitingSurface( FILE *tpsf );
	void minimize( int freeze_membrane  );
	void saveRestart( FILE *theFile, int seed ); 
	void saveRestart( char **buf, int seed);
	void loadRestart( FILE *loadFile, int *seed );
} Simulation;


#endif
