#ifndef __interph__
#define __interph__
#include <stdlib.h>
#include <stdio.h>
#include "boundary_conditions.h"
#include "Bin.h"
#include "irr_kernel.h"
#include "lagrange.h"
#include "sparse.h"
#include "lipid_composition.h"
#ifdef FFTW
#include <fftw3.h>
#endif
#define SAFETY (1e-12)
#define MAX_VALENCE 15
#define CURVE_DATA 6
#define EDGE_DATA 11
#define WDATA_SIZE 4 
#define PATCH_DATA 21
#define MAX_DOMAIN 100

// forward declarations

struct pcomplex;
struct parameterBlock;
struct srd_integrator;
//const double KA = 0.215;
//const double kc = 7.2;

int bmap( int i, int j );
int loadLattice( const char *fileName, double noise );

// the derivatives.

#define N_BLOCK_QUANTITIES 1


typedef struct
{
	int np;
	int *cp;
	double *w;
	double *rpbc;
	double r0[3];
	double k;
} constraintPoint;


typedef struct
{	
	int id;
	int valence;
	int nfaces; // valence-1
	double r[3];
	double ro[3];
	double c0;


	/* For computing curvatures */
	double c1, c2;
	double cdir1[3],cdir2[3];
	double nav;

	int edgeIndices[MAX_VALENCE];
	int edges[MAX_VALENCE];
	int edge_rev[MAX_VALENCE];

	double edge_curvature[MAX_VALENCE];
	double patch_A0[MAX_VALENCE];
	double edge_PBC[MAX_VALENCE*3];
	char broken[MAX_VALENCE];
	char broken_chain[MAX_VALENCE];
	char broken_looped[MAX_VALENCE];
	int tag;
	int    fixed;
	int protein_pt;
	double tan[3];
	double der2[3];
	double midp[3];
	double midp_tan[3];
	double midp_der2[3];
	
	double opp_midp[3];
	double opp_midp_tan[3];
	double opp_midp_der2[3];
	double store_val;
	double store_val2;
	double store_vector[3];
	double store_c[2];
	char flag;

	double nrm[3];
	void nearPoint( double *r, int *i_out, int *e_out, double *u_out, double *v_out );
	int faces[MAX_VALENCE*2];
	int *irr_coord_set; // for irrational vertices, store the coordinate sets.
} vertex;

#define MASK_1  (1<<0)
#define MASK_2  ((1<<1)+(1<<2))
#define MASK_3  (1<<3)

typedef struct
{ // cannot have point members, I do a bitwise copy in removeEdge
	int ids[3];
	int sense;
	int permutation;
	int edges[3];
	int ft_vector;
	int f; // which face this corresponds to.

	double f_lipids_stashed;
	double f_lipids;
	double nlipids;

	int draw;

	int edge_type; // MASK_1 + MASK_2 + MASK_3
	int border_tri[3];

	int npSpace;
	int np;
	int *plist;
	double *pc0; // c0 of embedded particle.
	double *pa; // area of embedded particle.

	double fixed_tilt[3];
	double pbc1[3];
	double pbc2[3];	

	localLipidComposition composition;	
	double dp;
	double nrm[3];
	double A0;

	double area( double *r, double PBC_vec[3][3], double alpha );
	void d_area( double *g, double *r, double PBC_vec[3][3], double alpha, double *d_a_lambda ); 
	void cen( double *ce, double *r, double PBC_vec[3][3], double alpha, double fi = 1.0/3.0, double fj = 1.0/3.0 );
	void d_cen( double *g, double *r, double PBC_vec[3][3], double alpha, double *d_c_lambda, double fi = 1.0/3.0, double fj = 1.0/3.0 ); 
		
	void normal( double *nrm, double *r, double PBC_vec[3][3], double alpha );
	void d_normal( double *g, double *r, double PBC_vec[3][3], double alpha, double *d_n_lambda ); 
} triangle;

typedef struct
{	// no pointers can be put into the edge: currently I do a bitwise copy in removeEdge.
	int faces[3]; // extra space for problem edges
	int vertices[2];

	int sense;
	int fix_sense;
	int code[3];
} edge;

struct volume_formula
{
	int ni;
	int ncoor;
	int *cp;
	double *r_w; // generating formula for triangular vertexes
	double *r_pbc; // pbc info
};

struct formula
{
	int vertex;
	int edge;
	int ncoor;

	double undeformed_r[3];
	double c0;
	double g0;
	double g0_base;
	double RuRv0;
	double *r_pbc; // the base PBC vector for this point.
	double *r_w; // the weights that determine the position.

	double *r_u; // the weights that determine dr_du
	double *r_v; // the weights that determine dr_dv

	double *r_uu; // the weights that determine dr_du
	double *r_uv; // the weights that determine dr_du
	double *r_vv; // the weights that determine dr_dv
	double weight;
	double orig_u, orig_v;
	
	int tri;
	int *cp;
};

struct boxel
{
	int np;
	int npSpace;
	int *f; // face
	int *p; // subpoint
};
struct boxing
{

	int nx;
	int ny;
	int nz;
	
	int *pbox; // id of box for each particle.	
	boxel *boxes;
	
};


struct surface
{
	int surface_id; // global id.
	bilayerComposition bilayerComp;		
	boundary_condition *bcs;
	int nbc;
	int disable_PBC;
	double *cumulative_area;
	int *cumulative_face;
	int fix_sense;
	int nfaces;

	int nfaces6;
	formula *single_formulas;
	
	int **s_faces_for_vertex;
	int *s_nfaces_for_vertex;
	int *s_nfacesSpace_for_vertex;

	int nfacesI;
	formula *irr_formulas; // NYI.

	int nf_faces;
	int nf_g_q_p;

	int nf_irr_faces;
	int nf_irr_pts;

	int **faces_for_vertex;
	int *nfaces_for_vertex;
	int *nfacesSpace_for_vertex;

	volume_formula *theVolumeFormulas;
	formula *edgeFormulas;
	formula *theFormulas;
	formula *theIrregularFormulas;
	constraintPoint *cons;
	int ncons;

	double c0;
	double PBC_vec[3][3];
	double PBC[3];
	int nv;
	int nt;
	int nts; // ntri, space
	int nedges;
	int nedgesSpace;
	int opencl_init;
	int max_valence; // the maximum valence of any vertex.

	int npatches;
	double *coords;
	edge *theEdges;
	triangle *theTriangles;
	vertex *theVertices;
	boxing *ptBoxes;
	boxing *prBoxes;	
	irr_kernel *kernels[20];

	// BEGIN flags from input 
        int on_surface;
	// END flags from input

	Bin *** box;
	int xbox;
	int ybox;
	int zbox;	

	double max_width_x;
	double max_width_y;
	double max_width_z;

	int xlength;
	//int xBlength;
	int ylength;
	//int yBlength;
	int zlength;
	//int zBlength;

	double std_metric;

	int loadLattice( const char *fileName, double noise, surface *copyFrom=NULL );
	int loadAndCopyLattice( const char *fileName, surface *altSurface );
	void setg0(double *r, double reset_factor = 0);
	double energy( double *r, double *puv, int do_vertex=-1, int *plist = NULL, int *np_found = NULL, int doMonge=0 );
	void getModifiedFaces( int do_vertex, int *flist, int *nf, int *plist, int *np );
	double faceEnergy( int f, double *r, double *puv, int doMonge );
	double fenergy( int f, double *r, double *puv );
	double ifenergy( int f, double *r, double *puv );
        double fenergym( int f, double *r, double *puv );
	double penergy( int f, double *r, double *puv, int do_monge, double face_energy_density  );
	void area( double *r, int do_vertex, double *area, double *area0);
	double energyMonge( double *r, int do_vertex=-1 );
	double irregularEnergy( double *r );
	void grad(double *r, double *g, double *puv=NULL, double *pg=NULL );
	void pgrad(double *r, double *g, double *puv, double *pg=NULL);
	// single particle gradient function.
	void particle_H_grad( double *r, double *gr, int f, double u, double v, double p_area, double p_c0, double *p_uv_g );
	void igrad(double *r, double *g);
	void pointGradient( int face, double u, double v, double *r, double *gr, double *g_puv, double *de_dr, double *de_dnrm, double frac_mult = 1.0 );
	void saveSurface( const char *fileName );
	void saveLimitingSurface( const char *fileName );
	void readCoords( double *r );
	void writeCoords( double *r );
	void writeSurface( const char *name );
	void put( double *r );
	void get( double *r );
	
	double faceArea( int i, int e );
	void d_faceArea( double *g, int i, int e );
	void norm( double *nrm, int i, int e  );
	void d_norm( double *g, int i, int e );	
	void ijk( int i, int e, int *ijkv );
	
	void writeToSurface( surface *lowerSurface, double thickness );
	void writeToSurfaceGen( surface *lowerSurface, double thickness, const char *job_name );
	void writeAreaEnergyVMD( double *r );
	double area_energy( double *r );
	void area_energy_grad( double *g, double *r );
	void loadStartingConditions( double thick );
	void subdivideSurface( surface *copy_surface );
	int subdivideSurface3( surface *copy_surface, int refine_at );
	int  subdividePath3( surface *copy_surface, int *path, int npath );
	void subdivideTriangle3( surface *copy_surface, int center_i, int center_e, int *new_center, int *new_edge  );
	void writeXYZSurface( const char *fileNameXYZ, const char *fileNamePSF, surface *theSurface );
	void writeLimitingSurface( FILE *theFile, pcomplex **allComplexes = NULL, int ncomplexes = 0, double *alphas =NULL);
	void writeLimitingSurfacePSF(FILE*, pcomplex **allComplexes = NULL, int ncomplexes = 0);
	void generatePlan( void );
	void writeStructure( FILE *theFile );
	void writeLimitStructure( FILE *theFile, int *pf=NULL, double *puv=NULL, int np=0, double dist_nrm=0 );
	void writeBARLimitStructure( FILE *theFile, int *pf, double *puv, int np );
	void writeLimitPSF( FILE *theFile, int np=0, double dist_nrm = 0 );
	void writeBARPSF( FILE *theFile, int np );
	void writeLimitTriangles( FILE *theFile );
	double externalPotential( double *r );
	void externalPotentialGrad( double *r, double *g );
	void getHeight( double *h, double *r, int nx, int ny );
	void getHeight2( double *h, double *r, int nx, int ny );
	void getHeight3( double *h, double *r, int nx, int ny );
	void dft( double *fh, double *fq, int x, int ny, double Lx, double Ly );
	double volume( double *r);
	double dvolume( double *r, double *g, double scale = 1.0);
	void generateVolumePlan(void);
	void randomPointOnSurface( int *face, double *u, double *v );
	void updateFaceInfoForRandom( void );
	void evaluateRNRM( int f, double u, double v, double *rp, double *nrm, double *r );
	double totalArea( void );
	double get_target_z( double x, double y );
	double directFT( double *hq, double *hq2, double *r, int nx, int ny );

	void getAMAT( double *AMAT, double *ro, double *r0_pos );
	int getFormulaAMAT( double **AMAT, double *ro, double **r0_pos, double **weights );
	int getSphericalHarmonicModes( double *ro, int l_min, int l_max, double q_max, double **gen_transform, double **output_qvals, double **scaling_factors );
	int origSphericalHarmonicModes( double *ro, int l_min, int l_max, double **gen_transform, double **output_qvals, double **scaling_factors );
	int getPlanarHarmonicModes( double *ro, int mode_x, int mode_y, int mode_min, int mode_max, double q_max, double **gen_transform, double **output_qvals, double **scaling_factors );
	int origPlanarHarmonicModes( double *ro, int mode_x, int mode_y, int mode_min, int mode_max, double **gen_transform, double **output_qvals, double **scaling_factors );

	// THESE SHOULD EVENTUALLY BE DELETED	
	void addParticleToFace( int f, int pid, double c0, double p_area );
	int getNear( double *p_in, int id_in, double *all_p, double cut, double alpha_x, double alpha_y, double alpha_z, int *plist, double *rads );
	void addParticle( double *r, int id, double alpha_x, double alpha_y, double alpha_z );
	void updateParticle( double *r, int id, double alpha_x, double alpha_y, double alpha_z );
	void setupBoxing( double *r, double *lcoords, int npts=0, double rad=0, int do_pr=0  ); // setup surface or pt boxing.


	void mode_perturb( double *r, double *ro, int qi, int qj, int nx, int ny, double Lx, double Ly, double mag, int do_cosine);
	void setup_mode_perturb( double *ro, int qi, int qj, int nx, int ny, double Lx, double Ly );
	void setup_spherical_perturb( double *ro, int ql, int qm );
	void load_least_squares_fitting( force_set *theForceSet );
	int setup_spherical_set( double *ro, int mode_min, int mode_max, double **output_qvals );
	int  setup_planar_set( double *ro, int l_min, int l_max, double **output_qvals );
	void mode_perturb( double *r, double mag, int do_cosine);
	void project_grad_to_modes( double *r, double *g, double *mg );
	double openCLEnergy( double *r, double *oarr );
	double fixEdgeGrad( double *r, double *g);
	double fixEdgePotential( double *r );
	double rhoEnergy( double * r, double PBC_vec[3][3], double thick_i, double thick_o );
	double rhoGrad( double *r, double *g, double PBC_vec[3][3], double thick_i, double thick_o, double *thick_i_g, double *thick_o_g  );
	double rhoWorker( double * r, double *gr, double PBC_vec[3][3], int do_grad, double thick_i, double thick_o, double *thick_i_g, double *thick_o_g );
	void removeEdge( int e );
	int createEdge( int i, int j, int k ); // creates a single triangle.
	int sealEdge( int i, int j, int k ); // creates two triangles
	void assignEdgeIndices( void );
	void assignEdgePBC( void );
	void validateCodes( void );
	void constructTriangles( void );
	void getVPassData( double **verts, int *nvert, int **tris, int ** eft_tris, int *ntri, int **edges, int *nedges, int edge_dim );
	void clean( void ); // removes edges based on edges with three faces
	void scrub( void ); // removes edges based on cartesian coordinates
	void setEdgeRev(void);
	void flipOrientation( void );
	void orientSurface( void );
	void surfaceEliminateHighValence( void );
	void nearestPointOnSurface( double *r, int *vertex  );
	void spatialConstraintGrad( double *r, double *g );
	double spatialConstraintEnergy( double * r);	
	void writeXYZandPSFPeriodic( const char *baseName );
	void writeVertexXYZandPSFPeriodic( const char *baseName );
	void getReport( double *r, double cmin, double cmax, int nbins, double *J_hist, double *K_hist, double *jtot, double *ktot);
	void nearPt( double *p, double *lcoords, double cut, double alpha_x, double alpha_y, double alpha_z, int *bf, int *bp );
	void storeDipoleOP( double *r, int f, int p, double *dr );
	void removeParticleFromFace( int f, int pid, double c0, double p_area );
	void computeEdgeCurvature( double *r );
	void sortFaces( void );
	void constructIrregularKernels( void );

	// BOUNDARY CONDITIONS
	double BCEnergy( double *r );
	void BCGrad( double *r, double *g);
	void perturbPlanarMiddle( double dnrm, double k);
	void zeroPlanarEdge( double k );
	void clearBoundaryConditions( void );
	void computeModifiedVertices( void );
	void duplicate( surface *dupe, int nx, int ny, int nz );

	int nextFace( int f, double *u, double *v, double *du, double *dv, double *r, double *mom=NULL, double *coord_transform=NULL);
	int neighborList( int f, int *neighbors );
	int map( int face_from, int face_to, double u, double v );

	void generateBorderMappings(void);
	double g( int f, double u, double v, double *r );
	double metric( int f, double u, double v, double *r, double *gmat, double *gmat_u, double *gmat_v );
	void ru( int f, double u, double v, double *r, double *dr_u );
	void rv( int f, double u, double v, double *r, double *dr_v );
	void r2der( int f, double u, double v, double *r, double *dr_uu, double *dr_uv, double *dr_vv );
	double c( int f, double u, double v, double *r, double *k, double *c_vec_1=NULL, double *c_vec_2=NULL, double *c_val_1=NULL, double *c_val_2=NULL);
	double dG( int f, double u, double v, double grad[5], double *r );
	void get_pt_coeffs( int f, double u, double v, double *coeffs, int *coord_list, int *ncoords );
        void get_pt_dcoeffs( int f, double u, double v, double *coeffs, int *coord_list, int *ncoords );

	/* gradFetch gets the mesh derivatives of r, n, and c1+c2 */
	double gradFetch( int f, double u, double v, double *r, 
		double *rGrad, // dim 3 by nCoor
		double *nGrad, // dim 3 by nCoor
		double *hGrad, // dim nCoor	
		double *kGrad, // dim nCoor	
		int *nCoor,
		int *vecList, double *k );

	void fetchPuv( int f, double u, double v, double *P_uv, double *rsurf );
	int fetchdP_duv( int f, double u, double v, double **dP_duv, double *rsurf );
	int fetchCiu( int f, double u, double v, double **C_iu, double *rsurf );
	
	void localMove( int *f_in, double *u_in, double *v_in, double sigma, double *r, double *frc_duv=NULL, double dt=1, double *fstep =NULL, int max_corr_iter=0);
	void localMoveReflect( int *f_in, double *u_in, double *v_in, double sigma, double *r, double *frc_duv=NULL, double dt=1, double *fstep =NULL, int max_corr_iter=0);
	double trialMove( int *f_in, double *u_in, double *v_in, double duv[2], double lambda, double *r, double *base_position );

	// COMPUTING S(q)
	void sample_B_hist( double *rmesh, double *B_hist, double *A2dz2_sampled, int sample_type, int nsamples, double max_r, int nbins, int shape_correction );
	void processSANS( parameterBlock *block, double **qvals, int *nq );	

//	void ComputeSq( double *Sq, double *r, double * p_r_m, int np, double b_av, double b_p, double q_min, double q_max, int n_q, int Sq_res, int do_planar );
//	void binCurvature( double *cdist, double minc, double maxc, int nbins, double *r, int Sq_res );
//	double *Sq_cache; // this is the ``resting area'' of the pieces of the membrane.
//	void cacheSqG(int sq_res, double *r);
//	void EiQFace( int f, double *Sq, double *qvals, int nq, double *r, double b_val );
//	void generateSubdivisionMatrices( double **M, int mlow, int mhigh);

//	void nearPointOnBoxedSurface( double *pt, int *f, double *u, double *v, double **M, int mlow, int mhigh, double *distance, double L=-1 ); 
	void nearPointOnBoxedSurface( double *pt, int *f, double *u, double *v, double **M, int mlow, int mhigh, double *distance, double initial_distance=-1, int disable_PBC_z=0 ); 
	int linearCollisionPoint( double *pt1_in, double *pt2_in, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, int disable_PBC_z=0 );
	void assembleNearList( double *pt, int **f_list_io, double **puv_list_io, double **areas_io, int *npts, double **M, int mlow, int mhigh, double Rmax, int sub_limit, int disable_PBC_z=0 );
	void buildFaceData( double **, int *, int *);
	bool withinRadius( double *pt, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, double L, double radius, double *vertex_data, int *ptr_to_data, int *nump, int disable_PBC_z=0 );
	bool withinSurface( double *pt, int *f, double *u, double *v, double **M, int mlow, int mhigh, double *distance ); 
	bool withinBoxedSurface( double *pt, int *f, double *u, double *v, double **M, int mlow, int mhigh, double *distance, double L=-1, int inside_outside=-1, int disable_PBC_z=0, double reflecting_surface=0); 
	void rebox_system(void);
	void box_system( double edge_length=-1);
	void findBox(double* pt, int *f, double *u, double *v, double **M, int mlow, int mhigh, double *distance, double L);
	double returnRadius (double *pt, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, double distance, double currentR, double L, double defaultR, double *vertex_data, int *ptr_to_data, int *nump, int inside_outside=-1, int disable_PBC_z=0, double reflecting_surface=0);
	void returnSizeofBox(int* x, int* y, int* z);
	void moveParticleonSurface(int *f, double *u, double *v, double *p);	
	void checkCurvature(FILE* file);
	void test_force_set( force_set *theForceSet );
	void evaluate_fp( force_set *theForceSet, double *xp, double *r );
	double evaluate_T(  double *vq, double *p, double *vq_prev=NULL, double *p_prev=NULL, double *alphas=NULL);
	void dKE_dx_and2( force_set *theForceSet,SparseMatrix *effM, double *vp, double *unit_f_vec, int f, double u, double v, double *dKEdx, double *d2KEdx2 );
	void getEffectiveMass( force_set * theForceSet, double *effective_mass );
	void getSparseEffectiveMass( force_set * theForceSet, int *use_map, int *nuse, SparseMatrix **, double *gen_transform=NULL, int NQ=-1, double *mass_scaling=NULL);
	void getSparseRoot( force_set * theForceSet, SparseMatrix **theMatrix, double *gen_transform=NULL, int NQ=-1 );
	void approxSparseEffectiveMass( force_set * theForceSet, int *use_map, int *nuse, SparseMatrix **, double *gen_transform=NULL, int NQ=-1, double *mass_scaling=NULL);
	void applyForceAtPoint( int f, double u, double v, double *dp, double *force_vector, force_set *theForceSet );
	void velocityAtPoint( int f, double u, double v, double *vp, double *velocity_out );
	void fdiff_check_grad( double *r );
	void testFAP( int f, double u, double v, double *dp, double *effM );

	void wrapPBC( double *dr, double *alphas=NULL );
	void loadComplexes( pcomplex ***allComplexes, int *ncomplex, parameterBlock *block );
	void minimize( double *r, pcomplex **allComplexes, int ncomplex, int do_freeze_membrane = 0 );
	void debug_dynamics( double *r, force_set *theForceSet, double *Minv, pcomplex **allComplexes, int ncomplex );
	void timestep_analysis( double *r, force_set *theForceSet, double *Minv, pcomplex **allComplexes, int ncomplex, double approx_timestep );

	// restarts
	void saveRestart( FILE *theFile, double *rsurf, double *pp, pcomplex **allComplexes, int ncomplex, int NQ, int seed );
	void saveRestart( char **buf, double *rsurf, double *pp, pcomplex **allComplexes, int ncomplex, int NQ, int seed );

	double cellVolume( void );
	void evaluate_momentum( force_set *theForceSet, double *vq, double *pmesh, double *pout );

	/* lipid movement */
	void lipidSync( void );
	void lipidBroadcast( void );
	void set_g0_from_f(void);
	void set_g0_from_f(int f);
	void stashf(void);
	void unstashf(void);
	void local_lipidMCMove( double *r, pcomplex **allComplexes, int ncomplex, double dt, double beta, int swap_only=0);
	void measureLipidCurvature( double *r, int pre_equil /* don't do running average */ ); // part of lipid redistribution


	/* writing tachyon object file*/
	void printTachyonRadius( FILE *theFile, double scale, int tachyon_collision_level, double *pt_in, int tachyon_collision_draw_type, int *col_f, double *col_u, double *col_v, double **M, int mlow, int mhigh, double L, double radius,	
		double *vertex_data, int *ptr_to_data, int *nump, int disable_PBC_z=0  ); // utility for graphics
	int writeTachyon( const char *name,
				int grid_uv, // resolution, how many points we do on one side of a triangle 
				int nint, // number of interpolation frames
				double *r, // current coordinates				
				pcomplex **allComplexes,
				int ncomplex,
				parameterBlock *params,
				srd_integrator *srd,
				int face_center=-1
			);
	
	void getRegions( int *regions, int nregions );
	void readLipidComposition( FILE *inputFile );
	void debugDeformation( double *r);
	void createAllAtom( parameterBlock *block );
	int evaluate_at( double eval[3], double dr[3], int f, double *u, double *v, double *rsurf, int leaflet, double strain, 
	double *dx_duv=NULL, double *dy_duv=NULL, double w_use=1.0, double w_rim=0.0, double *rimp=NULL, double *rimn = NULL );
	int simple_evaluate_at( double eval[3], double dr[3], int f, double *u, double *v, double *rsurf, int leaflet, double strain, 
	double *dx_duv=NULL, double *dy_duv=NULL, double w_use=1.0, double w_rim=0.0, double *rimp=NULL, double *rimn = NULL );
	int getCoordinateSystem( int source_f,   double *source_u,  double *source_v, 
				double *dr, double strain, int leaflet,
				  double *dx_duv, double *dy_duv, double *rsurf, int *regional_face, int *ncrosses );

#ifdef FFTW
	fftw_complex *h_in;
	fftw_complex *h_out;
	fftw_plan fftw_p;
#endif
};

struct volel
{
	int rfetch[6];
	int tets[3][4];
	double sgn[3];
	double pbc1[3];
	double pbc2[3];
	double l_v0[3];
	double v0;
	double volume( double *r, double PBC_vec[3][3], double alpha, int t_select = -1 );
	void dvolume( double *g, double *r, double PBC_vec[3][3], double alpha, double *d_v_lambda, int t_select = -1 );
};

struct volumes
{
	surface *surface1;
	surface *surface2;
	
	int nvol;
	volel *vols;

	void setup( surface *surface1In, surface *surface2In, double thick );
	double volume_energy( double *r );
	void volume_energy_grad( double *g, double *r );
	void writeVolumes( const char *fileNameXYZ, const char *fileNamePSF, surface *upperSurface, surface *lowerSurface, double *coords, double alpha );
};


int collision( surface *surface1, surface *surface2, double **M, int mlow, int mhigh );
int checkCollision( double *r1, int nv1, double *r2, int nv2, double **M, int mlow, int mhigh, int level, int max_level );
bool gjk_algorithm( double *r1, int nv1, double *r2, int nv2 );
int checkCollision2( double *r1, int nv1, double *r2, double **M, int mlow, int mhigh, int level, int max_level, double radius, double scalef, double cur_u, double cur_v, double *col_u, double *col_v );
bool minimum_distance( double *r1, int nv1, double *r2, double radius );
//void box_system(double L);
int cmpfunc (const void* a, const void* b);


#endif
