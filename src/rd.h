#ifndef __rdh__
#define __rdh__

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
	double prevnorm;
	double prevsep;
	double curr_sep;
	double currnorm;
	double ps_prev;
	int info;
} RD_tracked_info;

typedef struct
{
	int reacted;
	int pid; // index of pcomplex
	int sid; // index of site.
        int ntracked; // number of site-site distances tracked at the moment by this site.
        int ntracked_space;
        RD_tracked_info *tracked_info;
	RD_tracked_info *tracked_new;
	double Rmax_for_site;
} RD_tracked;

typedef struct
{
	int stype;
	char *complex_name;
} ReactantType;

struct RD
{
	parameterBlock *params;

	int nsites_tracked;
	int nsites_tracked_space;
	RD_tracked **tracked;	
	
	double dt;
	double Rmax_test;

	double max_binding_radius;
	double max_Rmax;	

	ReactantType *reactants;
	int nreactants;
	int nreactantsSpace;

	ReactionInformation *allReactions;
	int nreactions;
	int nreactionsSpace;	
	int *rxn_lookup_table;
	
	// BEGIN temporary storage for boxed particle tracking
	int ntotp;
	int *to_clear;
	int *complex_for_id;
        int *subp_for_id;
	int *nearlist;
	// END temporary storage for boxed particle tracking
	
	void do_rd( Simulation *theSimulation );
	void init(Simulation *theSimulation, double time_step, parameterBlock *pblock );
	
	// RD Setup:
	void parseRDFile( const char *fileName );
	void registerReaction( int site_type1, const char *pcomplex_name1,
			       int site_type2, const char *pcomplex_name2,
					double k_on, double k_off, double binding_radius, const char *productName );
	void registerSite( pcomplex *theComplex, int pid, int sid );
	

	// Tracking reactions/collisions
	void get_tracked( Simulation *theSimulation );
	void box_reactants( Simulation *theSimulation  );
	void unbox_reactants( void ); // cleans up boxing.
	int check_RD_blocked( Simulation *theSimulation, int p, int s);
	
	// Utility
	int rxnLookup( int code1, int code2 ); // returns -1 for no reaction, otherwise, site in the allReactions
	int siteLookup( const char *pname, int site_type );
	int getPComplexSiteID( const char *pcomplex_name, int site_type );
};

#endif
