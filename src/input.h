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
	char *jobName;
	double rho;


	int timestep_analysis;
	int mass_scaling;
	int mode_x;
	int mode_y;
	int mode_min;
	int mode_max;
	int monge;
	int nsteps;
	int nequil; // number of 
	int nmin; // number of steps of minimization
	int z_only;
	int silent;
	int nruns;

	int disable_mesh; // for debugging

	int lipid_mc_period; // period which to do monte-carlo lipid movement
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
	int kinetics;
	int kinetic_corr_period;
	int kinetics_do_phase;
	double diffc; //particle diffusion constant 
	double aqueous_diffc; //particle diffusion constant 
	double time_step;
	double time_step_collision; // srd collision timestep
	int non_interacting;

	int record_curvature; // record average curvature experienced.

	int tachyon_overlay_mesh;
	int tachyon_interp;
	int tachyon_tri_center;
	int tachyon_res;
	int tachyon;
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
	int do_bd;
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

	int on_surface;

	// BEGIN section for creating all-atom structures

	char *patchPDB;
	int create_all_atom;

	// END section for creating all-atom structures


	parameterBlock(void);
};

int getInput( const char **argv, int argc, parameterBlock *block); 
void setDefaults( parameterBlock *block );

#endif
