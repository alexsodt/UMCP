#ifndef __inputh__
#define __inputh__

struct complex_record
{
	char *name;
	int inside_outside;
	int nbound;
	int nsolution;
	double coverage;
	double concentration;

	struct complex_record *next;
};

class parameterBlock
{
	public:

	int defaults_set;
	
	complex_record *complex_records;

	char *loadName;
	char *meshName;
	char *meshName2;
	char *jobName;
	char *qvals;
	char *fitRho;
	char *lipid_lib;
	char *rxnDiffusionInfoName;
	char *track_lipid_rho;
	double rho;

	// gather
	int do_gather;
	char *dcdName;
	char *structureName;

	int restrain_volume_inside;
	int restrain_volume_outside;
	int timestep_analysis;
	int mass_scaling;
	int disable_mesh;
	int mode_x;
	int mode_y;
	int mode_min;
	int mode_max;
	double mode_q_max;
	int monge;
	int nsteps;
	int nequil; // number of 
	int nmin; // number of steps of minimization

	int minimizeResetG; // area-based viscous regularization

	int z_only;
	int silent;
	int nruns;

	int outputMesh;

	int lipid_mc_period; // period which to do monte-carlo lipid movement
	int lipid_mc_swap_only; // only swap lipids between the mesh (don't allow total movement of particles)
	int npt_mc_period; // period which to do monte-carlo lipid movement
	int cyl_tension_mc_period;


	double alpha_restraint_x;
	double alpha_restraint_y;
	double alpha_restraint_z;
	double alpha_restraint_k;
	int write_alpha_period;

	// umbrella sampling particle g(r)
	int bin;
	double del;

	// kinetic MC	
	int debug_diffusion;
	int kinetics;
	int kinetic_corr_period;
	int kinetics_do_phase;
	double diffc; //particle diffusion constant 
	double aqueous_diffc; //particle diffusion constant 
	double time_step;
	double time_step_collision; // srd collision timestep
	int non_interacting;

	int record_curvature; // record average curvature experienced.

	int tachyon_pbc_x;
	int tachyon_pbc_y;
	int tachyon_pbc_z;
	
	double tachyon_view_x;
	double tachyon_view_y;
	double tachyon_view_z;
	
	int tachyon_flip_sense;
	int tachyon_overlay_mesh;
	int tachyon_interp;
	int tachyon_tri_center;
	int tachyon_face_box_spline;
	int tachyon_res;
	int tachyon;
	int tachyon_curvature;
	int tachyon_gauss;
	int tachyon_dull;
	int tachyon_clear;

	double tachyon_collision_point[3];
	double tachyon_collision_radius;
	int tachyon_collision_level;
	int tachyon_collision_draw_type;

	int movie;
	int random_seed;
	int debug;
	int sphere;
	int o_lim;
	double perturbCenter;
	int fixEdge;

	int track_rho;

	double T;
	double mode_KA;
	double KA;
	double kv;
	double kc;
	double kg;
	double radius1;
	double dist_nrm;
	double radius2;
	double dimer_radius;
	double dimer_eps;
	double c0;
	double dimer_c0;
	double footprint;
	double leaflet_fraction;
	double concentration;
	int    mean_field;
	int correlated;

	int do_vol_inside;
	int do_vol_outside;

	double k_off;
	double k_on;
	double sigma;

	char *betazFile;
	double b_particle; // scattering length of ``particle'' -- probably sum over all lipid atoms minus background.	
	double b_av; // scattering length per unit area.
	int s_q;     // should we compute s_q?
	int nse;     // should we compute the auto-correlation function to mimic spin echo?
	double q_min; // q_min to compute s_q
	double q_max; // q_max to compute s_q
	int s_q_res;
	int nq;
	double maxr;
	double binw;
	double max_time; // the max time before it resets computation of the correlation function.
	int s_q_period;
	int shape_correction;

	int ncorr;     // should we compute the auto-correlation function of the membrane normal?

	double hours;


	// Parameters for the BAR domain, one of our initial validation models 

	double mab_k_theta;
	double mab_d_theta;
	double mab_bond_length;

	int do_bar;
	double bar_bond_length;
	double bar_bond_k;
	double bar_theta_0;
	double bar_theta_k;
	double bar_phi_k;
	double bar_phi_0;

	// crowder params

	double crowder_d; // particle membrane distance
	double crowder_bond_k; // bond constant for distance
	double crowder_r; // particle radius
	double crowder_attraction;		
	double crowder_attraction_r;	
	double crowder_attraction_r_rep; // radius at which repulsion starts	
	double crowder_mass;

	double default_bond_k;
	double mab_bond_k;
	
	// thermostats, LD and SRD
	int sub_com_period;
	int do_ld;
	int do_bd_particles;
	int do_bd_membrane;
	int nve_switch; // outer loop iteration to switch to NVE dynamics
	int hard_z_boundary;
	double srd_M;
	double gamma_langevin;	
	double hull_fudge;
	// stochastic rotational dynamics
	double eta;
	int do_srd;
	int srd_collision_freq;
	int planar_topology;
	int collect_hk;
	//reaction diffusion
	int do_rd;
	double bound_sigma;

	int on_surface;

	// BEGIN section for creating all-atom structures

	int addSalt;
	double innerKCL;
	double outerKCL;

	char *addProteinPDB;
	char *addProteinPSF;
	char *solvatePDB;
	char *solvatePSF;

	double strainInner;
	double strainOuter;

	char *innerPatchPDB;
	char *innerPatchPSF;
	char *outerPatchPDB;
	char *outerPatchPSF;
	char *altPatchPDB;
	char *altPatchPSF;
	char *patchPSF;
	char *patchPDB;
	double neutral_surface;
	double scale_solvent_approach;
	int create_all_atom;
	int create_flip;
	double create_pore;
	int do_rim;
	int perfect_solvent_tiling;

	double shift[3];

	// END section for creating all-atom structures

	double fitCoupling;
	double fitThickness;

	parameterBlock(void);
};

int getInput( const char **argv, int argc, parameterBlock *block); 
void setDefaults( parameterBlock *block );

#endif
