#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "input.h"
#include "util.h"
#include <math.h>
#include <sys/time.h>

parameterBlock::parameterBlock( void )
{
	defaults_set = 0;
}

void printParamBlock( parameterBlock *block );

void setDefaults( parameterBlock *block )
{
	block->defaults_set = 1;

	// defaults
	const char *defaultMesh = "planar.mesh";
	block->meshName = (char *)malloc( sizeof(char) * (strlen(defaultMesh)+1) );
	sprintf(block->meshName, "%s", defaultMesh );
	
	const char *defaultName = "default";
	block->jobName = (char *)malloc( sizeof(char) * (strlen(defaultName)+1) );
	sprintf(block->jobName, "%s", defaultName );

	block->loadName = NULL;
	block->concentration = 0;
	block->rho = 0;
	block->mode_x = -1;
	block->mode_y = -1;
	block->fix_alpha = 0;
	block->mode_KA = 0;
	block->nsteps = -1; // these are resolved later.
	block->nequil = 0;
	block->nmin = 0;
	block->KA = -1;
	block->kc = 14;
	block->kg = 0;
	block->mode_min = 2;
	block->mode_max = -1;
	block->T = 298;
	block->mab_k_theta = 1;
	block->mab_bond_length = 60;
	block->mab_d_theta = 10; // 10 degrees.

	block->sub_com_period = 10;
	block->complex_records = NULL;
	block->time_step_collision = -1;
	block->srd_M   = 2.0;
	block->hard_z_boundary = 0;
	block->do_srd = 0;
	block->do_ld = 0;
	block->nve_switch = -1; // switch to NVE dynamics at this timestep.
	block->gamma_langevin = 10; 
	block->planar_topology = 0;
	block->collect_hk = 0;
	block->nruns = 1;
	block->time_step = 1e-9;
	block->kinetics = 0;
	block->kinetics_do_phase = 0;
	block->diffc = 1e9; // angstroms^2/s
	block->aqueous_diffc = 1e10; // angstroms^2/s
	block->kinetic_corr_period = 10; 

	block->tachyon = 0;
	block->tachyon_res = 3;
	block->tachyon_interp = 1;

	block->lipid_mc_period = 0;
	block->npt_mc_period = 0;
	block->cyl_tension_mc_period = 0;
	block->alpha_restraint_x = -1;
	block->alpha_restraint_y = -1;
	block->alpha_restraint_z = -1;
	block->alpha_restraint_k = 0;

	block->leaflet_fraction = 0.5;

	block->track_rho = 0;

	block->monge = 0;
	block->z_only = 0;
	block->radius1 = sqrt(65/M_PI);
	block->dist_nrm = 0;
	block->radius2 = 0;
	block->dimer_radius = 0;
	block->dimer_eps = 0;
	block->c0 = 0;
	block->dimer_c0 = 0;
	block->footprint = 65;
	block->silent = 0;
	block->movie = 0;
	block->debug = 0;    
	block->sphere = 0;
	block->hours = -1; // use steps instead.
	block->o_lim = 10000;

	block->eta = 8.9e-4; // SI units
	block->hull_fudge = 1.0;
	block->fixEdge = 0;
	block->perturbCenter = 0;

	block->do_bar		= 0;
	block->bar_bond_length = 75;
	block->default_bond_k      = 1;
	
	block->bar_bond_k      = -1;
	block->mab_bond_k      = -1;
	block->bar_theta_0     = (360.0 - 137.0)/2;
	block->bar_theta_k     = 1000;
	block->bar_phi_k       = 1000;
	block->bar_phi_0       = 180;

	block->crowder_bond_k       = -1;
	block->crowder_r       = 25;
	block->crowder_d       = 25;
	block->crowder_attraction = 0;
	block->crowder_attraction_r = -1;

	block->mean_field  = 0;
	block->k_off = 1e4;
	block->k_on  = 1e11; 
	block->bin = 0;
	block->del = 3.0;
	block->sigma = 10.0;

	block->Sq_res = 5;	
	block->s_q = 0;
	block->b_particle = -2.325053524; // protiated POPC 
	block->b_av = 0.281376608; // perdeuterated POPC
	block->nse = 0;
	block->ncorr = 0;
	block->q_min = 0.001;
	block->q_max = 0.05;
	block->nq    = 100;
	block->non_interacting = 1;
	block->max_time = 1; // one second.
	
	// request a timestep analysis
	block->timestep_analysis = 0;

    	struct timeval tp;
        gettimeofday( &tp, NULL );
	block->random_seed = tp.tv_usec ;
}

int resolveParameters( parameterBlock *block )
{
	int warning = 0;

	if( block->mab_bond_k < 0 )
		block->mab_bond_k = block->default_bond_k;
	if( block->bar_bond_k < 0 )
		block->bar_bond_k = block->default_bond_k;
	if( block->crowder_bond_k < 0 )
		block->crowder_bond_k = block->default_bond_k;
	if( block->crowder_attraction != 0 && block->crowder_attraction_r < 0 )
		block->crowder_attraction_r = block->crowder_r;	
	if( block->npt_mc_period > 0 && block->cyl_tension_mc_period > 0 )
	{
		printf("As implemented, NPT and cylindrical tension are incompatible.\n");
		exit(1);
	}
	if( block->mode_x != -1 && block->mode_y == -1 )
		block->mode_y = 0;
	if( block->mode_y != -1 && block->mode_x == -1 )
		block->mode_x = 0;	

	if( block->mode_x >= 0 || block->mode_y >= 0 || block->mode_max >= 0 )
	{	// doing simple mode moves, should converge rapidly.
		if( block->nsteps == -1 )
			block->nsteps = 1000;
	}

	if( block->kinetics && (block->mode_x == -1 && block->mode_max == -1 ) )
	{
		printf("Kinetics with plain vertex moves not yet implemented.\n");
		exit(1);
	}

	if( block->nsteps == -1 )
		block->nsteps = 10000;

	if( block->nequil == -1 )
		block->nequil = block->nsteps * 0.1;

	if( block->KA < 0 && ((block->mode_x != -1 && block->mode_y != -1) || block->mode_max >0) )
		block->KA = 0;
	else if ( block->KA < 0 )
		block->KA = 0.215;
	
	if( block->dist_nrm > 0 && block->radius2 > 0 && (block->dist_nrm < block->radius1 + block->radius2) )
	{
		printf("ERROR: Particles 1 and 2 on the same inclusion overlap (r1+r2 < dist_nrm).\n");
		exit(1);
	}

	if( block->mode_x == -1 && block->mode_y == -1 && block->KA < 1e-8 && block->z_only == 0 && block->mode_max == -1)
	{
		printf("ERROR: KA is too small for lateral control point Monte Carlo moves.\n");
		exit(1);
		return 1;
	}

	if( block->sphere && !block->fix_alpha )
	{
		printf("ERROR: alpha is floating while simulating a sphere.\n");
		exit(1);
	}

	if( block->sphere && (block->mode_x >= 0 || block->mode_max >= 0) )
	{
		if( block->mode_x < 2 && block->mode_max < 2 )
		{
			printf("Mode spherical harmonics requires l > 1 (l=%d).\n", block->mode_x );
			exit(1);
		}

		if( block->mode_max >= 0 )
		{
		}
		else if( block->mode_y > block->mode_x || -block->mode_y > block->mode_x )
		{
			printf("Mode spherical harmonics requires m <= l (l=%d, m=%d).\n", block->mode_x, block->mode_y );
			exit(1);
		}
	}
	
	if(block->sphere && block->mode_max >= 0 && (block->mode_min < 2 || block->mode_min > block->mode_max) )
	{
		printf("Invalid mode_min/max: %d %d.\n", block->mode_min, block->mode_max );
		if( block->mode_min < 2 )
			printf("Min cannot be less than 2 for spherical harmonics.\n");
		exit(1);
	}

	if( block->sphere && block->monge )
	{
		printf("Error: The Monge gauge has been selected in sphere mode.\n");
		exit(1);
	}

	if( (block->mode_max >= 0 || block->mode_x >= 0) && block->KA != 0 )
	{
		printf("WARNING: Using modes but also applying a local KA.\n");
		warning += 1;
	}

	if( block->nse ) block->s_q = 1;

	if( block->time_step_collision < 0 )
		block->time_step_collision = block->time_step;

	if( block->do_ld && block->do_srd )
	{
		printf("ERROR: LD and SRD are both activated.\n");
		exit(1);
	}

	return warning;
}

#define FILE_PASS 		0
#define COMMAND_LINE_PASS 	1
#define END_PASS		2

int getInput( const char **argv, int argc, parameterBlock *block) 
{
	char *word1 = (char *)malloc( sizeof(char) * 4096 );
	char *word2 = (char *)malloc( sizeof(char) * 4096 );

	int is_file = 1;

	if( argc > 1 )
	{
		for( int t = 0; t < strlen(argv[1]); t++ )
		{
			if( argv[1][t] == '=' )
				is_file = 0;
		}
	}	

	FILE *theFile = NULL;
	int pass = FILE_PASS;

	int c = 2;
	if( argc <= 1 )
	{
		printf("Using defaults.\n");
		pass = END_PASS;
	}
	else if( !is_file )
	{
		printf("Interpreting the first argument as setting a parameter, not an input file.\n");
		pass = COMMAND_LINE_PASS;
		c = 1;
	}
	else
	{
		theFile = fopen(argv[1],"r");
		if( !theFile )
		{
			printf("Couldn't open input file '%s'\n", argv[1] );
			exit(1);
		}
	}
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	
	if( !block->defaults_set )
		setDefaults(block);

	if( !strcasecmp( argv[0] + strlen(argv[0])-3, "bar" ) )
		block->do_bar = 1;
	int ERROR = 0;

	while( pass != END_PASS )
	{
		const char *tbuf;

		if( pass == FILE_PASS )
		{
			getLine( theFile, buffer );
		
			if( feof(theFile) )
			{
				pass = COMMAND_LINE_PASS;
				continue;
			}
	
			if( strlen(buffer) > 4095 )
			{
				printf("Line %s too long.\n", buffer );
				ERROR = 1;
			}

			if( buffer[0] == '\0' )
				continue;

			int nr = sscanf( buffer, "%s %s", word1, word2 );
	
			if( !strcasecmp( word1, "lipid" ) ) continue; //special case handled in lipid_composition.C

			if( nr != 2 )
			{
				printf("Could not interpret input line '%s'.\n", buffer );
				ERROR = 1;
			}

			tbuf = buffer;
		
			const char *p = tbuf;
			while( *p == '\t' || *p == ' ' ) p += 1;

			if( *p == '#' ) 
			{
				continue;
			}
		}
		else if( pass == COMMAND_LINE_PASS )
		{
			if( c >= argc )
			{
				pass = END_PASS; 
				continue;
			}
			word1[0] = '\0';
			word2[0] = '\0';
			int x;
			for( x = 0; x < strlen(argv[c]); x++ )
			{
				if( argv[c][x] == '=' )
					break;
				word1[x] = argv[c][x];
				word1[x+1] = '\0';
			}
			x++;
			for( int y = 0; x < strlen(argv[c]); x++, y++ )
			{
				word2[y] = argv[c][x];
				word2[y+1] = '\0';
			}

			
			if( strlen(word1) < 1 || strlen(word2) < 1 )
			{
				printf("Couldn't parse command line argument '%s'.\n", argv[c] );
				ERROR = 1;
			}

			tbuf = argv[c];
			c++;
		}
		else break;


		if( !strcasecmp( word1, "lipid" ) )
		{
			 //special case handled in lipid_composition.C
		}
		else if( !strcasecmp( word1, "mesh" ) )
		{
			free(block->meshName);
			block->meshName = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->meshName, word2 );
		}
		else if( !strcasecmp( word1, "jobname" ) )
		{
			free(block->jobName);
			block->jobName = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->jobName, word2 );
		}
		else if( !strcasecmp( word1, "load" ) )
		{
			free(block->loadName);
			block->loadName = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->loadName, word2 );
		}
		else if( !strcasecmp( word1, "rho" ) )
			block->rho = atof( word2 );
		else if( !strcasecmp( word1, "concentration" ) )
		{
			block->concentration = atof( word2 );
			// input as molar.
			block->concentration *= (1e-30) * (6.022e23); // particles per cubic angstrom.
		}
		else if( !strcasecmp( word1, "lipid_mc_period" ) )
			block->lipid_mc_period = atoi( word2 );
		else if( !strcasecmp( word1, "npt_mc_period" ) )
			block->npt_mc_period = atoi( word2 );
		else if( !strcasecmp( word1, "cyl_tension_mc_period" ) )
			block->cyl_tension_mc_period = atoi( word2 );
		else if( !strcasecmp( word1, "alpha_restraint_x" ) )
			block->alpha_restraint_x = atof( word2 );
		else if( !strcasecmp( word1, "alpha_restraint_y" ) )
			block->alpha_restraint_y = atof( word2 );
		else if( !strcasecmp( word1, "alpha_restraint_z" ) )
			block->alpha_restraint_z = atof( word2 );
		else if( !strcasecmp( word1, "alpha_restraint_k" ) )
			block->alpha_restraint_k = atof( word2 );
		else if( !strcasecmp( word1, "nruns" ) )
			block->nruns = atoi( word2 );
		else if( !strcasecmp( word1, "mode_x" ) || !strcasecmp( word1, "mode_l" ))
			block->mode_x = atoi( word2 );
		else if( !strcasecmp( word1, "mode_min" ) )
			block->mode_min = atoi( word2 );
		else if( !strcasecmp( word1, "mode_max" ) )
			block->mode_max = atoi( word2 );
		else if( !strcasecmp( word1, "mode_y" ) || !strcasecmp( word1, "mode_m" ) )
			block->mode_y = atoi( word2 );
		else if( !strcasecmp( word1, "mode_KA" ) )
			block->mode_KA = atof( word2 );
		else if( !strcasecmp( word1, "nsteps" ) )
			block->nsteps = atoi( word2 );
		else if( !strcasecmp( word1, "nequil" ) )
			block->nequil = atoi( word2 );
		else if( !strcasecmp( word1, "nmin" ) )
			block->nmin = atoi( word2 );
		else if( !strcasecmp( word1, "o_lim" ) )
			block->o_lim = atoi( word2 );
		else if( !strcasecmp( word1, "seed" ) )
			block->random_seed = atoi( word2 );
		else if( !strcasecmp( word1, "T" ) )
			block->T = atof( word2 );
		else if( !strcasecmp( word1, "KA" ) )
			block->KA = atof( word2 );
		else if( !strcasecmp( word1, "kc" ) )
			block->kc = atof( word2 );
		else if( !strcasecmp( word1, "kg" ) )
			block->kg = atof( word2 );
		else if( !strcasecmp( word1, "perturb_center" ) )
			block->perturbCenter = atof( word2 );
		else if( !strcasecmp( word1, "radius1" ) )
			block->radius1 = atof( word2 );
		else if( !strcasecmp( word1, "dist_nrm" ) )
			block->dist_nrm = atof( word2 );
		else if( !strcasecmp( word1, "radius2" ) )
			block->radius2 = atof( word2 );
		else if( !strcasecmp( word1, "dimer_radius" ) )
			block->dimer_radius = atof( word2 );
		else if( !strcasecmp( word1, "dimer_eps" ) )
			block->dimer_eps = atof( word2 );
		else if( !strcasecmp( word1, "c0" ) )
			block->c0 = atof( word2 );
		else if( !strcasecmp( word1, "footprint" ) )
			block->footprint = atof( word2 );
		else if( !strcasecmp( word1, "dimer_c0" ) )
			block->dimer_c0 = atof( word2 );
		else if( !strcasecmp( word1, "hours" ) )
			block->hours = atof( word2 );
		else if( !strcasecmp( word1, "eta" ) )
			block->eta = atof( word2 );
		else if( !strcasecmp( word1, "hull_fudge" ) )
			block->hull_fudge = atof( word2 );
		else if( !strcasecmp( word1, "sub_com_period" ) )
			block->sub_com_period = atoi( word2 );
		// Crowder parameters
		else if( !strcasecmp( word1, "crowder_r" ) )
			block->crowder_r = atof( word2 );
		else if( !strcasecmp( word1, "crowder_attraction_r" ) )
			block->crowder_attraction_r = atof( word2 );
		else if( !strcasecmp( word1, "crowder_attraction" ) )
			block->crowder_attraction = atof( word2 );
		else if( !strcasecmp( word1, "crowder_d" ) )
			block->crowder_d = atof( word2 );
		else if( !strcasecmp( word1, "crowder_bond_k" ) )
			block->crowder_bond_k = atof( word2 );
		// BAR domain parameters
		else if( !strcasecmp( word1, "bar_bond_length" ) )
			block->bar_bond_length = atof( word2 );
		else if( !strcasecmp( word1, "mab_bond_length" ) )
			block->mab_bond_length = atof( word2 );
		else if( !strcasecmp( word1, "bar_bond_k" ) )
			block->bar_bond_k = atof( word2 );
		else if( !strcasecmp( word1, "default_bond_k" ) )
			block->default_bond_k = atof( word2 );
		else if( !strcasecmp( word1, "mab_bond_k" ) )
			block->mab_bond_k = atof( word2 );
		else if( !strcasecmp( word1, "mab_d_theta" ) )
			block->mab_d_theta = atof( word2 );
		else if( !strcasecmp( word1, "mab_k_theta" ) )
			block->mab_k_theta = atof( word2 );
		else if( !strcasecmp( word1, "bar_theta_0" ) )
			block->bar_theta_0 = atof( word2 );
		else if( !strcasecmp( word1, "bar_theta_k" ) )
			block->bar_theta_k = atof( word2 );
		else if( !strcasecmp( word1, "bar_phi_k" ) )
			block->bar_phi_k = atof( word2 );
		else if( !strcasecmp( word1, "bar_phi_0" ) )
			block->bar_phi_0 = atof( word2 );
		// end of BAR domain parameters
		else if( !strcasecmp( word1, "bin" ) )
			block->bin = atoi( word2 );
		else if( !strcasecmp( word1, "del" ) )
			block->del = atof( word2 );
		// end of umbrella sampling.
		else if( !strcasecmp( word1, "leaflet_fraction" ) )
			block->leaflet_fraction = atof( word2 );
		else if( !strcasecmp( word1, "gamma_langevin" ) )
			block->gamma_langevin = atof( word2 );
		else if( !strcasecmp( word1, "srd_M" ) )
			block->srd_M = atof( word2 );
		else if( !strcasecmp( word1, "nve_switch" ) )
			block->nve_switch = atoi( word2 );
		// SRD input parameters
		else if( !strcasecmp( word1, "do_srd" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->do_srd = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->do_srd = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "hard_z_boundary" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->hard_z_boundary = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->hard_z_boundary = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "do_ld" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->do_ld = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->do_ld = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "timestep_analysis" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->timestep_analysis = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->timestep_analysis = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "planar_topology" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->planar_topology = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->planar_topology = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "collect_hk" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->collect_hk = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->collect_hk = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		// end of SRD parameters
		else if( !strcasecmp( word1, "track_rho" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->track_rho = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->track_rho = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "non_interacting" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->non_interacting = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->non_interacting = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "s_q" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->s_q = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->s_q = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "max_time" ) )
			block->max_time = atof( word2 );
		else if( !strcasecmp( word1, "ncorr" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->ncorr = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->ncorr = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "nse" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->nse = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->nse = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "q_min" ) )
			block->q_min = atof( word2 );
		else if( !strcasecmp( word1, "q_max" ) )
			block->q_max = atof( word2 );
		else if( !strcasecmp( word1, "nq" ) )
			block->nq = atoi( word2 );
		else if( !strcasecmp( word1, "kinetic_corr_period" ) )
			block->kinetic_corr_period = atoi(word2);
		else if( !strcasecmp( word1, "diffc" ) )
			block->diffc = atof(word2);
		else if( !strcasecmp( word1, "Sq_res" ) )
			block->Sq_res = atoi( word2 );
		else if( !strcasecmp( word1, "b_av" ) )
			block->b_av = atof( word2 );
		else if( !strcasecmp( word1, "b_particle" ) )
			block->b_particle = atof( word2 );
		else if( !strcasecmp( word1, "time_step" ) )
			block->time_step = atof( word2 );
		else if( !strcasecmp( word1, "time_step_collision" ) )
			block->time_step_collision = atof( word2 );
		else if( !strcasecmp( word1, "k_off" ) )
			block->k_off = atof( word2 );
		else if( !strcasecmp( word1, "k_on" ) )
			block->k_on = atof( word2 );
		else if( !strcasecmp( word1, "sigma" ) )
			block->sigma = atof( word2 );
		else if( !strcasecmp( word1, "kinetics" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->kinetics = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->kinetics = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "do_phase" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->kinetics_do_phase = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->kinetics_do_phase = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "silent" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->silent = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->silent = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "debug" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->debug = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->debug = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "movie" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->movie = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->movie = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_res" ) )
			block->tachyon_res = atoi( word2 );
		else if( !strcasecmp( word1, "tachyon_interp" ) )
			block->tachyon_interp = atoi( word2 );
		else if( !strcasecmp( word1, "sphere" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->sphere = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->sphere = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "fix_edge" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->fixEdge = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->fixEdge = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "fix_alpha" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->fix_alpha = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->fix_alpha = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "monge" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->monge = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->monge = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "mean_field" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->mean_field = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->mean_field = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "add") )
		{
			// format: add COMPLEX_NAME nbound        %d  [inside|outside]
			// format: add COMPLEX_NAME nsolution     %d  [inside|outside]
			// format: add COMPLEX_NAME coverage      %lf [inside|outside]
			// format: add COMPLEX_NAME concentration %lf [inside|outside]

			complex_record *rec = (complex_record *)malloc( sizeof(complex_record)  );

			rec->next = block->complex_records;
			block->complex_records = rec;

			rec->inside_outside = 0;
			rec->nbound = 0;
			rec->nsolution = 0;
			rec->coverage = 0;
			rec->concentration = 0;
			rec->name = (char *)malloc( sizeof(char) * (1+ strlen(word2) ) );
			strcpy( rec->name, word2 );

			char *word3 = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );	
			char *word4 = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );	
			char *word5 = (char *)malloc( sizeof(char) * (strlen(buffer)+1) );	

			int nr = sscanf( buffer, "%s %s %s %s %s", word1, word2, word3, word4, word5 );

			if( nr < 4 )
			{
				printf("Syntax error in ``add'' command.\n");
				printf("add COMPLEX nbound|nsolution|coverage|concentration value [inside|outside]\n");
					exit(1);
			}

			double value = atof(word4);
			double uvalue = value;

			if( !strcasecmp(word3, "nbound") ) 
			{
				uvalue = lround(value);
				rec->nbound = value;
			}
			else if( !strcasecmp( word3, "nsolution") )
			{ 
				uvalue = lround(value);
				rec->nsolution = uvalue;
			}
			else if( !strcasecmp( word3, "coverage") ) 
				rec->coverage = value;
			else if( !strcasecmp( word3, "concentration") ) 
				rec->concentration = value;
			
			if( nr == 5 )
			{
				if( !strcasecmp( word5, "inside") ) 
					rec->inside_outside = -1;
				else if( !strcasecmp( word5, "outside" ) )
					rec->inside_outside = 1;
				else
				{
					printf("Error interpreting optional inside/outside fifth argument in command '%s'.\n", buffer );
					exit(1);
				}
			}

			if( fabs( value -uvalue) > 1e-10 )
			{
				printf("ERROR. likely floating point/integer error for command '%s'.\n", buffer );
				exit(1);
			}
		}
		else
		{
			printf("Could not interpret input line '%s'.\n", tbuf );
			ERROR = 1;
		}
	}
	if ( ERROR )
	{
		exit(1);
	}
	int were_warnings = resolveParameters(block);

	printParamBlock(block);

	free(word1);
	free(word2);
	free(buffer);
	
	return were_warnings;
}

void printParamBlock( parameterBlock *block )
{
	printf("Parameters:\n");
	printf("/****************************/\n");
	printf("General:\n");
	printf("\tnequil:    %d\n", block->nequil );
	if( block->hours > 0 )
		printf("\thours:    %lf\n", block->hours );
	else
		printf("\tnsteps:    %d\n", block->nsteps );
	printf("Surface:\n");
	printf("\tMesh:      %s\n", block->meshName );	
	printf("\tmode_x:    %d\n", block->mode_x );	
	printf("\tmode_y:    %d\n", block->mode_y );	
	printf("\tfix_alpha  %s\n", (block->fix_alpha ? "yes" : "no" ) );
	printf("\tMonge      %s\n", (block->monge ? "yes" : "no" ) );
	printf("\tmode_KA:   %lf kcal/mol/A^2\n", block->mode_KA );
	printf("\tKA:        %lf kcal/mol/A^2\n", block->KA );
	printf("\tkc:        %lf kcal/mol\n", block->kc );
	printf("Particle:\n");
	printf("\trho:       %lf\n", block->rho );	
	printf("\tradius1:   %lf A\n", block->radius1 );
	printf("\tc0:        %lf A\n", block->c0 );
	printf("\tfootprint: %lf A^2\n", block->footprint );

	if( block->fixEdge )
		printf("Zeroing the surrounding edge of the bilayer.\n");
	if( fabs(block->perturbCenter) > 1e-10 )
		printf("Perturbing the middle of the bilayer by %lf Angstroms.\n", block->perturbCenter );


	if( block->dist_nrm <= 0 || block->radius2 < 1e-7)
	{
		printf("Not using a second particle off the membrane surface.\n");
	}
	else
	{
		printf("dist_nrm:  %lf A\n", block->dist_nrm );
		printf("radius2 :  %lf A\n", block->dist_nrm );
	}

	if( block->dimer_radius <= block->radius1 )
	{
		printf("Not modeling particle dimerization [dimer radius %lf less than particle radius %lf].\n", block->dimer_radius, block->radius1);
	}
	else
	{
		printf("Modeling particle dimerization:\n");
		printf("dimer radius:  %lf A\n", block->dimer_radius );
		printf("dimer eps:     %lf kcal/mol\n", block->dimer_eps );
	}
	if(  block->do_bar )
	{
		printf("/****************************/\n");
		printf("BAR domain parameters\n");
		printf("	bar_bond_length %lf Angstroms\n", block->bar_bond_length ); 
		printf("	bar_bond_k      %lf kcal/mol/Angstroms^2\n", block->bar_bond_k ); 
		printf("	bar_theta_0     %lf degrees\n", block->bar_theta_0 ); 
		printf("    bar_phi_0       %lf degrees\n", block->bar_phi_0 );
		printf("	bar_theta_k     %lf kcal/mol/degrees^2\n", block->bar_theta_k ); 
		printf("	bar_phi_k       %lf kcal/mol/degrees^2\n", block->bar_phi_k ); 

		printf("/****************************/\n");
	}
	if( block->kinetics )
	{
		printf("/****************************/\n");
		printf("Doing kinetics analysis.\n");
		if( block->kinetics_do_phase )
			printf("Doing phase analysis (sin and cos together).\n");
		else
			printf("Only doing the cos() mode.\n");
		printf("Kinetics diffusion coefficient %le\n", block->diffc );
		printf("Kinetics period %d.\n", block->kinetic_corr_period ); 
		printf("/****************************/\n");
		
	}
	if( block->hours > 0 )
	{
		printf("Running for %lf hours (instead of %d steps).\n", block->hours, block->nsteps );
	}
	printf("Random seed: %d\n", block->random_seed );

	if( block->lipid_mc_period > 0 )
		printf("Lipid Monte Carlo move period: %d\n", block->lipid_mc_period );
	if( block->npt_mc_period > 0 )
		printf("NPT Monte Carlo move period: %d\n", block->npt_mc_period );	
	
}




