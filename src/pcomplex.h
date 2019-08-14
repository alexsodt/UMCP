#ifndef __pcomplexh__

#include "simulation.h"
#include "interp.h"
#include "input.h"

#define __pcomplexh__

#define DEBUG_OFF	0
#define DEBUG_NO_V	1
#define DEBUG_NO_T	2

struct pcomplex
{
	int disabled;
	int nwatchers;

	int is_inside;
	int nsites;
	int nattach;
	int bound;
	int debug;
	int do_bd;

	char *complex_name;
	
	// diffusion constants.
	double *DC;

	// all masses.
	double *mass;

	// all WCA/LJ radii
	double *sigma;
	double *att_eps;
	double *att_sigma;

	// all coordinates in 3D
	double *rall;
	// attach coordinates in f/uv; can exceed face bounds
	int *sid; // id of surface we are on.
	int *stype; // if doing rxn/diff
	int *fs;	
	double *puv;

	// in the final face bounds.
	int *grad_fs;	
	double *grad_puv;
	
	double *last_p; // Hamiltonian conjugate momenta
	double *p; // Hamiltonian conjugate momenta
	double *qdot; // time derivatives of generalized coordinates.
	double *save_grad;
	
	double *PBC_ext; // the current PBC vector of each site.
	double *last_pos;

	double *p_area;
	double *p_c0;

	double *coord_transform;

// This is added by Kayla for keeping track of when to "delete" a complex
	int alive; //1: alive, 0:dead -- remove from memory during garbage cleanup  
// End of what Kayla added

	virtual void alloc( void );
	virtual void pfree( void );

	virtual int getNBonds( void ) { return 0; };
	virtual void putBonds( int *bond_list );

	virtual void init( surface *theSurface, double *, int f, double u, double v ); 
	virtual void init( double *r ); 

	virtual void bind( int f, double u, double v );
	virtual void unbind( void );

	virtual void applyParamsToComplex( double *p );
	virtual void getParamsFromComplex( double *p );
	virtual int nparams( void );
	virtual void orient( surface *theSurface, double *rsurf ) { }

	virtual void refresh( Simulation *theSimulation  );
	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation *theSimulation, double *surface_g, double *particle_g ); 
	virtual void fd_grad_debug(surface *theSurface, double *rsurf ); 
	virtual void loadParams( parameterBlock *block );

	virtual void move_inside( void );
	virtual void move_outside( void );

	virtual int packLenF( void );
	virtual int packLenI( void );
	virtual void packageF( double * );
	virtual void unpackF( double * );
	virtual void packageI( int *);
	virtual void unpackI( int *); 
	
	virtual int isElastic(void) { return 0; }	

	void activateBrownianDynamics(void);

	void printType( char **type );
	void cacheVelocities( void );
	void prepareForGradient( void );
	void loadCoords( surface *theSurface, double *rsurf, double *r, double *n );
	void setrall( Simulation *theSimulation  );
	void base_init( void );
	double T(Simulation *theSimulation, int subp=-1 );
	double update_dH_dq( Simulation *theSimulation, double time_step=-1, double timestep_total=0 );
	void propagate_p( Simulation *theSimulation, double dt );
	//void compute_qdot( surface *theSurface, double *rsurf, double *mesh_qdot0, double *mesh_qdot, double frac_mult=1.0 );
	void compute_qdot( Simulation *theSimulation,  double frac_mult=1.0 );
	void propagate_surface_q( Simulation *theSimulation, double dt );
	void copyParentParameters( pcomplex *parent );

	void debug_dPinv_dq( surface * theSurface, double *rsurf  );
	void getMeshQxdot( surface *theSurface, double *rsurf, double *Minv, double *mesh_p, double *mesh_qdot, double *mesh_qdot0, double *mesh_der_qdot );
	void applyLangevinFriction( Simulation *theSimulation, double dt, double gamma );
	void applyLangevinNoise( Simulation *theSimulation, double dt, double gamma, double temperature );
	int saveComplex( char *buffer, int *buflen, int bufLenMax );
	void saveComplex( FILE * theFile );
	void loadComplex( FILE * theFile, Simulation *theSimulation, int load_to );
	double AttachV( Simulation *theSimulation );
	double AttachG( Simulation *theSimulation, double *pg );
	void evaluate_momentum( surface *theSurface, double *rsurf, double *pout );
	double local_curvature( Simulation *theSimulation);
	void print_type( char **outp );


	void disable( void) { disabled = 1; }
	void destroy( void);
	void watch( void );
	void forget( void );
};

struct actin : pcomplex
{
	int *interfaces; // array of avaiable interfaces by sub ID (for no branching - "0" and "nattach-1")
	double *r_int; // position of the available interfaces
	double k_on;
	double k_off;
};

struct simpleParticle : pcomplex
{
};

struct simpleLipid : pcomplex
{
	double c0_val;
	void loadParams( parameterBlock *block );
	virtual void init( surface *theSurface, double *, int f, double u, double v ); 
};

struct simpleDimer : simpleLipid
{
	double c0_val;
};

struct NBAR : pcomplex
{
	double bond_length;
	double k_phi;
	double k_theta;
	double bond_k;
	double theta_0;
	double phi_0;


	void init( surface *, double *rsurf, int f, double u, double v ); 
	void init( double *r );
	void bind( int f, double u, double v);
	void unbind( void );
	void loadParams( parameterBlock *block );
	void orient( surface *theSurface, double *rsurf );
	
	int getNBonds( void );
	void putBonds( int *bond_list );

	double V( Simulation *theSimulation );
	double grad( Simulation *theSimulation, double *surface_g, double *particle_g );

	void move_inside(void);
	void move_outside(void);	
};

struct dimer : pcomplex
{
	double bond_length;
	double bond_k;
	void init( surface *, double *rsurf, int f, double u, double v ); 
	void init( double *r );
	void loadParams( parameterBlock *block );
	
	int getNBonds( void );
	void putBonds( int *bond_list );
	
	 void bind( int f, double u, double v );
	 void unbind( void );

	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation*theSimulation, double *surface_g, double *particle_g );
};

struct MAB : dimer
{
	// dtheta_0 is the optimal angle between the normals, projected along the vector between the particles
	double dtheta_0;
	// k_theta is the force constant for dtheta_0.
	double k_theta;

	void move_inside( void );
	void move_outside( void );

	void loadParams( parameterBlock * block );
	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation *theSimulation, double *surface_g, double *particle_g );
};

struct elasticCrowder : pcomplex
{
	double hard_r;
	double d; // the distance from the surface to the sphere center.
	double bond_k;
	double att_eps_val;
	double att_sigma_val;
	double att_sigma_inner;
	double crowder_mass;

	void init( surface *, double *rsurf, int f, double u, double v ); 
	void loadParams( parameterBlock * block );
	int isElastic(void) { return 1; }	
	void move_inside( void );
	void move_outside( void );
	int getNBonds( void ) { return 1; }
	void putBonds( int *bond_list );
	virtual double V( Simulation *theSimulation );
	virtual double grad( Simulation *theSimulation, double *surface_g, double *particle_g );
};

void propagateSolutionParticles( Simulation *theSimulation, double dt );
pcomplex *loadComplex( const char *name );

#endif
