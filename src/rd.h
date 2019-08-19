#ifndef __rd_kaylah__
#define __rd_kaylah__

#define SITE_TYPE_ANY -1

#include "interp.h"
#include "input.h"
#include "fpr_subroutines/2D.h"

// forward declarations..
struct parameterBlock;
struct Simulation;
// ..done


typedef struct
{
	double k_on;		// on rate vol / s depending on the equilibrium.
	double k_off;		// off rate, per second
	double Dtot;		// sum of diffusion constants.
	double Rmax;	       // max distance to look.
	double binding_radius; // distance between sites
	char *pcomplex_name1;
	char *pcomplex_name2;
	char *productName;
	int site1;
	int site2;
} ReactionInformation;

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
	int pid; // index of pcomplex
	int sid; // index of site.
        int ntracked; // number of site-site distances tracked at the moment by this site.
        int ntracked_space;
        RD_tracked_info *tracked_info;
	RD_tracked_info *tracked_new;
} RD_tracked;

typedef struct
{
	int stype;
	char *complex_name;
} ReactantType;

struct RD
{
	int nsites_tracked;
	int nsites_tracked_space;
	RD_tracked **tracked;	
	
	double dt;
	double Rmax_test;
	
	ReactantType *reactants;
	int nreactants;
	int nreactantsSpace;

	ReactionInformation *allReactions;
	int nreactions;
	int nreactionsSpace;	
	int *rxn_lookup_table;

	void parseRDFile( const char *fileName );
	void registerReaction( int site_type1, const char *pcomplex_name1,
			       int site_type2, const char *pcomplex_name2,
					double k_on, double k_off, double binding_radius, const char *productName );
	void registerSite( pcomplex *theComplex, int pid, int sid );
	
	// returns -1 for no reaction, otherwise, site in the allReactions
	int rxnLookup( int code1, int code2 );
	int siteLookup( const char *pname, int site_type );
	
	int getPComplexSiteID( const char *pcomplex_name, int site_type );
	void init(Simulation *theSimulation, double time_step, parameterBlock *pblock );
	void get_tracked( Simulation *theSimulation );
	void do_rd( Simulation *theSimulation );
};

#endif
