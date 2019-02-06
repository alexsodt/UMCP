#ifndef __srdh__
#define __srdh__

#include "lagrange.h"

#define POINT_INSIDE -1
#define POINT_UNKNOWN 0
#define POINT_OUTSIDE 1


typedef struct 
{
	int *list;
	int np;
	int npSpace;
} bin;

typedef struct srd_integrator
{
	int np;
	double *rp;
	double *rp_prev; // previous position for tracking collisions
	double *rp_dist; // last position for which we calculated a membrane-particle distance
	double *vp;
	double *distance; // min distance to the membrane.
	double a; // box length approx
	double temp;
	double eta;
	double dt_c;
	double mass;
	double membrane_point_mass;
	double srdCollisionEstimator;

	double membrane_min[3];
	double membrane_max[3];
	
	double dKE_in_collision;
	double dKE_out_collision;
	
	double dKE_in_collision_sweep;
	double dKE_out_collision_sweep;

	double dKE_in_thermostat;
	double dKE_out_thermostat;

	int planar;
	int debug_mode;
	double *topological_tag;
	double *alt_topological_tag;
	double *last_known_tag;
	double *alt_last_known_tag;
	int *did_collide;
	double rotation_angle;

	int hard_z_boundary;
	int *inside_outside; // particle is inside or outside topologically speaking.
	int grain_x;
	int grain_y;
	int grain_z;
	double PBC_vec[3][3];	
	bin *bins;
	

//	void init( double mesh_spacing, double PBC_vec[3][3], double temp_in, double eta_in, double dt_si, double dt, int doPlanarTopology, double binw=50.0 );
	void init( double mesh_spacing, double PBC_vec[3][3], double temp_in, double eta_in, double dt_si, double dt, int doPlanarTopology, double mass, double M, int hard_z_boundary );
	void clear(void);

	void setMembraneMinMax(double *r, surface *theSurface);
	void clearDKE( void ) { dKE_in_thermostat=0; dKE_out_thermostat=0; dKE_in_collision_sweep=0; dKE_out_collision_sweep=0; dKE_in_collision=0; dKE_out_collision=0; }
	double KE(void);
	void activateDebugMode( void) { debug_mode = 1; }
	void writeXYZ( FILE *theFile );
	void stream( double dt );
	void collide( void );
	int stream_and_collide( double *r, double *g, double *drdt, double *vmem, double *effm, surface *theSurface, double **M  , int mlow, int mhigh, double delta_hull, force_set *theForceSet , double running_time, double integrate_time, double dt_col );
	int collide_with_membrane( double *r, double *g, surface *theSurface, double **M  , int mlow, int mhigh, double delta_hull, force_set *theForceSet , double time_step, double force_factor, double *vmem, double *effm );
	int collide_with_membrane_inside_outside( double *r, double *g, surface *theSurface, double **M  , int mlow, int mhigh, double delta_hull, force_set *theForceSet , double time_step, double force_factor, double *vmem, double *effm );
	int collide_with_membrane_planar( double *r, double *g, surface *theSurface, double **M  , int mlow, int mhigh, double delta_hull, force_set *theForceSet , double time_step, double force_factor, double *vmem, double *effm );
	int collide_with_membrane_planar_multi_collision( double *r, double *g, surface *theSurface, double **M  , int mlow, int mhigh, double delta_hull, force_set *theForceSet , double time_step, double force_factor, double *vmem, double *effm );
	void initializeDistances( double *r, surface *theSurface, double **M, int mlow, int mhigh );
	void tagParticlesForCollision( double *r, surface * theSurface, double delta_hull_collision, double **M, int mlow, int mhigh );
	void altTag( double *r, surface * theSurface, double delta_hull_collision, double **M, int mlow, int mhigh );
	void checkResolve( double *r, surface * theSurface, double delta_hull_collision, double **M, int mlow, int mhigh );
	int resolveCollision( double *r, double *g, double *vp, double *effm, double delta_hull_collision, surface *theSurface, double **M, int mlow, int mhigh, force_set *theForceSet, double dt, double force_factor );
	void subPos(double*dr);
} srd_integrator;

#endif
