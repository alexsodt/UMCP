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

	const char *defaultBZ = "beta.z";
	block->betazFile = (char *)malloc( sizeof(char) * ( strlen(defaultBZ)+1) );
	sprintf(block->betazFile, "%s", defaultBZ );

	block->fitRho = NULL;
	block->fitCoupling = 1.0;
	block->fitThickness = 15.0;

	block->rxnDiffusionInfoName = NULL;

	block->qvals = NULL;
	block->loadName = NULL;
	block->concentration = 0;
	block->rho = 0;
	block->lipid_lib = NULL;
	block->meshName2 = NULL;

	block->record_curvature = 0;
	block->create_all_atom = 0;	
	block->create_flip = 0;
	block->do_rim = 0;
	block->perfect_solvent_tiling = 0;

	block->mode_x = -1;
	block->mode_y = -1;
	block->mode_KA = 0;
	block->nsteps = -1; // these are resolved later.
	block->nequil = 0;
	block->nmin = 0;
	block->minimizeResetG = 0;

	block->KA = -1;
	block->kc = 14;
	block->kg = 0;
	block->mode_min = -1;
	block->mode_max = -1;
	block->mode_q_max = -1;
	block->T = 298;
	block->mab_k_theta = 1;
	block->mab_bond_length = 60;
	block->mab_d_theta = 10; // 10 degrees.

	block->shift[0] = 0;
	block->shift[1] = 0;
	block->shift[2] = 0;

	block->outputMesh = 0;

	block->lipid_lib = NULL;

	block->sub_com_period = 10;
	block->complex_records = NULL;
	block->time_step_collision = -1;
	block->srd_M   = 2.0;
	block->hard_z_boundary = 0;
	block->do_srd = 0;
	block->do_ld = 0;
	block->do_bd_membrane = 0;
	block->do_bd_particles = 0;


	block->do_rd = 0;
	block->bound_sigma = 10.0;

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

	block->disable_mesh = 0;

	block->tachyon_collision_draw_type = 0;
	block->tachyon_collision_point[0] = 0;
	block->tachyon_collision_point[1] = 0;
	block->tachyon_collision_point[2] = 0;
	block->tachyon_collision_radius = -1;
	block->tachyon_collision_level = 3;

	block->tachyon_face_box_spline = -1;
	block->tachyon_dull = 0;
	block->tachyon_clear = 0;
	block->tachyon_flip_sense = 0;
	block->tachyon_clear = 0;
	block->tachyon_curvature = 0;
	block->tachyon_gauss = 0;
	block->tachyon_view_x = 0;
	block->tachyon_view_y = 0;
	block->tachyon_view_z = 0;
	
	block->tachyon_pbc_x = 0;
	block->tachyon_pbc_y = 0;
	block->tachyon_pbc_z = 0;

	block->tachyon = 0;
	block->tachyon_overlay_mesh = 0;
	block->tachyon_tri_center = -1;
	block->tachyon_res = 3;
	block->tachyon_interp = 1;

	block->lipid_mc_period = 0;
	block->lipid_mc_swap_only = 0;
	block->npt_mc_period = 0;
	block->cyl_tension_mc_period = 0;
	block->alpha_restraint_x = -1;
	block->alpha_restraint_y = -1;
	block->alpha_restraint_z = -1;
	block->alpha_restraint_k = 0;
	block->write_alpha_period = -1;

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
	block->crowder_attraction_r_rep = -1;
	block->crowder_mass	= 2e5;

	block->mean_field  = 0;
	block->correlated  = 0;
	block->k_off = 1e4;
	block->k_on  = 1e11; 
	block->bin = 0;
	block->del = 3.0;
	block->sigma = 10.0;

	block->s_q_res = 5;	
	block->s_q = 0;
	block->b_particle = -2.325053524; // protiated POPC 
	block->b_av = 0.281376608; // perdeuterated POPC
	block->nse = 0;
	block->ncorr = 0;
	block->q_min = 0.001;
	block->q_max = 0.5;
	block->nq    = 1000;
	block->non_interacting = 1;
	block->max_time = 1; // one second.
	block->s_q_period = 1000; // 1000 steps per S_q	
	block->shape_correction = 0;

	block->on_surface = 0;
	
	block->addSalt = 0;
	block->innerKCL = 0;
	block->outerKCL = 0;

	block->addProteinPDB   = NULL;
	block->addProteinPSF   = NULL;

	block->strainInner = 0;
	block->strainOuter = 0;

	block->innerPatchPDB = NULL;
	block->innerPatchPSF = NULL;
	
	block->outerPatchPDB = NULL;
	block->outerPatchPSF = NULL;
	
	block->altPatchPDB = NULL;
	block->altPatchPSF = NULL;

	block->patchPDB = NULL;
	block->patchPSF = NULL;

	block->solvatePDB = NULL;
	block->solvatePSF = NULL;

	block->do_gather = 0;
	block->structureName = NULL;
	block->dcdName = NULL;


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
	if( block->crowder_attraction_r > 0 && block->crowder_attraction_r_rep < 0 )
		block->crowder_attraction_r_rep = 0.75 * block->crowder_attraction_r;	
	if( block->npt_mc_period > 0 && block->cyl_tension_mc_period > 0 )
	{
		printf("As implemented, NPT and cylindrical tension are incompatible.\n");
		exit(1);
	}
	if( block->mode_x != -1 && block->mode_y == -1 )
		block->mode_y = 0;
	if( block->mode_y != -1 && block->mode_x == -1 )
		block->mode_x = 0;	

	if( block->mode_x >= 0 || block->mode_y >= 0 || block->mode_max >= 0 || block->mode_q_max >= 0 )
	{	// doing simple mode moves, should converge rapidly.
		if( block->nsteps == -1 )
			block->nsteps = 1000;
	}

	if( block->kinetics && (block->mode_x == -1 && block->mode_max == -1 && block->mode_q_max < 0 ) )
	{
		printf("Kinetics with plain vertex moves not yet implemented.\n");
		exit(1);
	}

	if( block->nsteps == -1 )
		block->nsteps = 10000;

	if( block->nequil == -1 )
		block->nequil = block->nsteps * 0.1;

	if( block->KA < 0 )
		block->KA = 0.215;

//	if( block->KA < 0 && ((block->mode_x != -1 && block->mode_y != -1) || block->mode_max >0) )
//		block->KA = 0;
//	else if ( block->KA < 0 )
//		block->KA = 0.215;
	
	if( block->dist_nrm > 0 && block->radius2 > 0 && (block->dist_nrm < block->radius1 + block->radius2) )
	{
		printf("ERROR: Particles 1 and 2 on the same inclusion overlap (r1+r2 < dist_nrm).\n");
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
	
	if( block->sphere && block->mode_max >= 0 && block->mode_min < 0 )
		block->mode_min = 2;
	if( !block->sphere && block->mode_max >= 0 && block->mode_min < 0 )
		block->mode_min = 0;

	if( block->mode_min > block->mode_max )
	{
		printf("Invalid mode_min/max: %d %d.\n", block->mode_min, block->mode_max );
		if( block->mode_min < 2 )
			printf("Min cannot be less than 2 for spherical harmonics.\n");
		exit(1);
	}


	if(block->sphere && block->mode_max >= 0 && (block->mode_min < 2 || block->mode_min > block->mode_max) )
	{
		printf("Invalid mode_min/max: %d %d.\n", block->mode_min, block->mode_max );
		exit(1);
	}

	if( block->sphere && block->monge )
	{
		printf("Error: The Monge gauge has been selected in sphere mode.\n");
		exit(1);
	}

	if( (block->mode_max >= 0 || block->mode_x >= 0 || block->mode_q_max >= 0 ) && block->KA != 0 )
	{
		printf("WARNING: Using modes but also applying a local KA.\n");
		warning += 1;
	}

	if( block->nse ) block->s_q = 1;

	if( block->time_step_collision < 0 )
		block->time_step_collision = block->time_step;
	if( block->do_ld && (block->do_bd_membrane && block->do_bd_particles) )
	{
		printf("ERROR: Langevin and Brownian dynamics are both activated.\n");
		exit(1);
	}

	if( block->do_ld && block->do_srd )
	{
		printf("ERROR: LD and SRD are both activated.\n");
		exit(1);
	}


	if( !block->outerPatchPDB )
	{
		if( block->patchPDB )
		{
			block->outerPatchPDB = (char *)malloc( sizeof(char) * (1+strlen(block->patchPDB) ) );
			strcpy( block->outerPatchPDB, block->patchPDB );

			if( block->patchPSF )
			{
				block->outerPatchPSF = (char *)malloc( sizeof(char) * (1+strlen(block->patchPSF) ) );
				strcpy( block->outerPatchPSF, block->patchPSF );
			}
		}
	}
	
	if( !block->innerPatchPDB )
	{
		if( block->patchPDB )
		{
			block->innerPatchPDB = (char *)malloc( sizeof(char) * (1+strlen(block->patchPDB) ) );
			strcpy( block->innerPatchPDB, block->patchPDB );

			if( block->patchPSF )
			{
				block->innerPatchPSF = (char *)malloc( sizeof(char) * (1+strlen(block->patchPSF) ) );
				strcpy( block->innerPatchPSF, block->patchPSF );
			}
		}
	}
	
	if( block->do_rim && !block->altPatchPDB )
	{
		if( block->patchPDB )
		{
			block->altPatchPDB = (char *)malloc( sizeof(char) * (1+strlen(block->patchPDB) ) );
			strcpy( block->altPatchPDB, block->patchPDB );

			if( block->patchPSF )
			{
				block->altPatchPSF = (char *)malloc( sizeof(char) * (1+strlen(block->patchPSF) ) );
				strcpy( block->altPatchPSF, block->patchPSF );
			}
		}
		else
		{
			printf("Adding a rim patch to a membrane structure requires either altPatchPDB or patchPDB to be set.\n");
			exit(1);
		}
	}
	
	if( block->do_gather && (!block->structureName) )
	{
		printf("Gathering requires a structure (e.g., PSF) file.\n");
		exit(1);
	}
	
	if( block->do_gather && (!block->dcdName) )
	{
		printf("Gathering requires a coordinate (e.g., dcd) file.\n");
		exit(1);
	}

	return warning;
}

#define FILE_PASS 		0
#define COMMAND_LINE_PASS 	1
#define END_PASS		2

int getInput( const char **argv, int argc, parameterBlock *block) 
{
	printf("Input line:\n");
	for( int x = 0; x < argc; x++ )
	{
		printf("%s ", argv[x] );
	}
	printf("\n");

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
		printf("Input file:\n");
		printf("#####################################\n");
	}
	char *buffer = (char *)malloc( sizeof(char) * 100000 );
	
	if( !block->defaults_set )
		setDefaults(block);

	int ERROR = 0;

	while( pass != END_PASS )
	{
		const char *tbuf;

		if( pass == FILE_PASS )
		{
			getLine( theFile, buffer );
		
			if( feof(theFile) )
			{
				printf("#####################################\n");
				pass = COMMAND_LINE_PASS;
				continue;
			}
			printf("%s\n", buffer );	
			if( strlen(buffer) > 4095 )
			{
				printf("Line %s too long.\n", buffer );
				ERROR = 1;
			}

			if( buffer[0] == '\0' )
				continue;
			
			tbuf = buffer;
			const char *p = tbuf;
			while( *p == '\t' || *p == ' ' ) p += 1;

			if( *p == '#' ) 
			{
				continue;
			}
				
			int nr = sscanf( buffer, "%s %s", word1, word2 );
	
			if( !strcasecmp( word1, "lipid" ) ) continue; //special case handled in lipid_composition.C
			
			if( nr != 2 )
			{
				printf("Could not interpret input line '%s'.\n", buffer );
				ERROR = 1;
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
		else if( !strcasecmp( word1, "structure" ) ) // for gathering.
		{
			if( block->structureName ) free(block->structureName);
			block->structureName = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->structureName, word2 );
		}
		else if( !strcasecmp( word1, "dcd" ) ) // for gathering.
		{
			if( block->dcdName ) free(block->dcdName);
			block->dcdName = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->dcdName, word2 );
		}
		else if( !strcasecmp( word1, "fitRho" ) )
		{
			if( block->fitRho )
				free(block->fitRho);
			block->fitRho = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->fitRho, word2 );
		}
		else if( !strcasecmp( word1, "mesh2" ) )
		{
			if( block->meshName2 )
				free(block->meshName2);
			block->meshName2 = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->meshName2, word2 );
		}
		else if( !strcasecmp( word1, "rxn_diffusion" ) )
		{
			if( block->rxnDiffusionInfoName )
				free(block->rxnDiffusionInfoName);
			block->rxnDiffusionInfoName = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->rxnDiffusionInfoName, word2 );
		}
		else if( !strcasecmp( word1, "lipid_lib" ) )
		{
			if( block->lipid_lib )
				free(block->lipid_lib);
			block->lipid_lib = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->lipid_lib, word2 );
		}
		else if( !strcasecmp( word1, "betaz" ) )
		{
			free(block->betazFile);
			block->betazFile = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->betazFile, word2 );
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
		else if( !strcasecmp( word1, "qvals" ) )
		{
			free(block->qvals);
			block->qvals = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->qvals, word2 );
		}
		else if( !strcasecmp( word1, "rho" ) )
			block->rho = atof( word2 );
		else if( !strcasecmp( word1, "concentration" ) )
		{
			block->concentration = atof( word2 );
			// input as molar.
//			block->concentration *= (1e-30) * (6.022e23); // particles per cubic angstrom.
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
		else if( !strcasecmp( word1, "write_alpha_period" ) )
			block->write_alpha_period = atoi( word2 );
		else if( !strcasecmp( word1, "nruns" ) )
			block->nruns = atoi( word2 );
		else if( !strcasecmp( word1, "lipid_mc_swap_only" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->lipid_mc_swap_only = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->lipid_mc_swap_only = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "minimizeResetG" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->minimizeResetG = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->minimizeResetG = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "gather" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->do_gather = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->do_gather = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "outputMesh" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->outputMesh = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->outputMesh = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "mass_scaling" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->mass_scaling = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->mass_scaling = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "record_curvature" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->record_curvature = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->record_curvature = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "disable_mesh" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->disable_mesh = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->disable_mesh = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "mode_x" ) || !strcasecmp( word1, "mode_l" ))
			block->mode_x = atoi( word2 );
		else if( !strcasecmp( word1, "mode_min" ) )
			block->mode_min = atoi( word2 );
		else if( !strcasecmp( word1, "mode_max" ) )
			block->mode_max = atoi( word2 );
		else if( !strcasecmp( word1, "mode_q_max" ) )
			block->mode_q_max = atof( word2 );
		else if( !strcasecmp( word1, "mode_y" ) || !strcasecmp( word1, "mode_m" ) )
			block->mode_y = atoi( word2 );
		else if( !strcasecmp( word1, "mode_KA" ) )
			block->mode_KA = atof( word2 );
		else if( !strcasecmp( word1, "nsteps" ) || !strcasecmp( word1, "nouter")  )
			block->nsteps = atoi( word2 );
		else if( !strcasecmp( word1, "nequil" ) )
			block->nequil = atoi( word2 );
		else if( !strcasecmp( word1, "nmin" ) )
			block->nmin = atoi( word2 );
		else if( !strcasecmp( word1, "o_lim" ) ||  !strcasecmp( word1, "ninner") )
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
		else if( !strcasecmp( word1, "crowder_attraction_r_rep" ) )
			block->crowder_attraction_r_rep = atof( word2 );
		else if( !strcasecmp( word1, "crowder_attraction" ) )
			block->crowder_attraction = atof( word2 );
		else if( !strcasecmp( word1, "crowder_d" ) )
			block->crowder_d = atof( word2 );
		else if( !strcasecmp( word1, "crowder_bond_k" ) )
			block->crowder_bond_k = atof( word2 );
		else if( !strcasecmp( word1, "crowder_mass" ) )
			block->crowder_mass = atof( word2 );
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
		else if( !strcasecmp( word1, "do_bd_particles" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
			{
				block->do_bd_particles = 1;
			}
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
			{
				block->do_bd_particles = 0;
			}
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}
		}
		else if( !strcasecmp( word1, "do_bd_membrane" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
			{
				block->do_bd_membrane = 1;
			}
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
			{
				block->do_bd_membrane = 0;
			}
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
		else if( !strcasecmp( word1, "do_rd" ) )
                {
                        if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
                                block->do_rd = 1;
                        else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
                                block->do_rd = 0;
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
		else if( !strcasecmp( word1, "shape_correction" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->shape_correction = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->shape_correction = 0;
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
		else if( !strcasecmp( word1, "aqueous_diffc" ) )
			block->aqueous_diffc = atof(word2);
		else if( !strcasecmp( word1, "s_q_res" ) )
			block->s_q_res = atoi( word2 );
		else if( !strcasecmp( word1, "s_q_period" ) )
			block->s_q_period = atoi( word2 );
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
		else if( !strcasecmp( word1, "tachyon_clear" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_clear = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_clear = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_dull" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_dull = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_dull = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_flip_sense" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_flip_sense = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_flip_sense = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_clear" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_clear = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_clear = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_pbc_x" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_pbc_x = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_pbc_x = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_pbc_y" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_pbc_y = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_pbc_y = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_pbc_z" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_pbc_z = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_pbc_z = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_curvature" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_curvature = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_curvature = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_gauss" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_gauss = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_gauss = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_overlay_mesh" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->tachyon_overlay_mesh = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->tachyon_overlay_mesh = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "tachyon_res" ) )
			block->tachyon_res = atoi( word2 );
		else if( !strcasecmp( word1, "tachyon_tri_center" ) )
			block->tachyon_tri_center = atoi( word2 );
		else if( !strcasecmp( word1, "tachyon_face_box_spline" ) )
			block->tachyon_face_box_spline = atof( word2 );
		else if( !strcasecmp( word1, "tachyon_view_x" ) )
			block->tachyon_view_x = atof( word2 );
		else if( !strcasecmp( word1, "tachyon_view_y" ) )
			block->tachyon_view_y = atof( word2 );
		else if( !strcasecmp( word1, "tachyon_view_z" ) )
			block->tachyon_view_z = atof( word2 );
		else if( !strcasecmp( word1, "tachyon_interp" ) )
			block->tachyon_interp = atoi( word2 );
		else if( !strcasecmp( word1, "tachyon_collision_level" ) )
			block->tachyon_collision_level = atoi( word2 );
		else if( !strcasecmp( word1, "tachyon_collision_radius" ) )
			block->tachyon_collision_radius = atof( word2 );
		else if( !strcasecmp( word1, "tachyon_collision_draw_type" ) )
			block->tachyon_collision_draw_type = atoi( word2 );
		else if( !strcasecmp( word1, "tachyon_collision_x" ) )
			block->tachyon_collision_point[0] = atof( word2 );
		else if( !strcasecmp( word1, "tachyon_collision_y" ) )
			block->tachyon_collision_point[1] = atof( word2 );
		else if( !strcasecmp( word1, "tachyon_collision_z" ) )
			block->tachyon_collision_point[2] = atof( word2 );
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
		else if( !strcasecmp( word1, "correlated" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->correlated = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->correlated = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "on_surface" ) )
		  {
		    if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
		      block->on_surface = 1;
		    else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
		      block->on_surface = 0;
		    else
		      {
			printf("Could not interpret input line '%s'.\n", tbuf );
			ERROR = 1;
		      }
		  }
		else if( !strcasecmp( word1, "shiftx" ) )
			block->shift[0] = atof(word2);
		else if( !strcasecmp( word1, "shifty" ) )
			block->shift[1] = atof(word2);
		else if( !strcasecmp( word1, "shiftz" ) )
			block->shift[2] = atof(word2);
		else if( !strcasecmp( word1, "do_rim" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->do_rim = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->do_rim = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "create_all_atom" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->create_all_atom = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->create_all_atom = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "perfect_solvent_tiling" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->perfect_solvent_tiling = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->perfect_solvent_tiling = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "create_flip" ) )
		{
			if( !strcasecmp( word2, "TRUE" ) || !strcasecmp( word2, "yes") || !strcasecmp( word2, "on" ) )
				block->create_flip = 1;
			else if( !strcasecmp( word2, "FALSE" ) || !strcasecmp( word2, "no") || !strcasecmp( word2, "off" ) )
				block->create_flip = 0;
			else
			{
				printf("Could not interpret input line '%s'.\n", tbuf );
				ERROR = 1;
			}	
		}
		else if( !strcasecmp( word1, "strainOuter") )
			block->strainOuter = atof(word2);
		else if( !strcasecmp( word1, "strainInner") )
			block->strainInner = atof(word2);
		else if( !strcasecmp( word1, "innerKCL" ) )
		{
			block->addSalt = 1;
			block->innerKCL = atof( word2 );
		}
		else if( !strcasecmp( word1, "outerKCL" ) )
		{
			block->addSalt = 1;
			block->outerKCL = atof( word2 );
		}
		else if( !strcasecmp( word1, "altPatchPDB" ) )
		{
			if( block->altPatchPDB )
				free(block->altPatchPDB);
			block->altPatchPDB = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->altPatchPDB, word2 );
		}
		else if( !strcasecmp( word1, "altPatchPSF" ) )
		{
			if( block->altPatchPSF )
				free(block->altPatchPSF);
			block->altPatchPSF = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->altPatchPSF, word2 );
		}
		else if( !strcasecmp( word1, "outerPatchPDB" ) )
		{
			if( block->outerPatchPDB )
				free(block->outerPatchPDB);
			block->outerPatchPDB = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->outerPatchPDB, word2 );
		}
		else if( !strcasecmp( word1, "outerPatchPSF" ) )
		{
			if( block->outerPatchPSF )
				free(block->outerPatchPSF);
			block->outerPatchPSF = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->outerPatchPSF, word2 );
		}
		else if( !strcasecmp( word1, "addProteinPDB" ) )
		{
			if( block->addProteinPDB )
				free(block->addProteinPDB);
			block->addProteinPDB = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->addProteinPDB, word2 );
		}
		else if( !strcasecmp( word1, "addProteinPSF" ) )
		{
			if( block->addProteinPSF )
				free(block->addProteinPSF);
			block->addProteinPSF = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->addProteinPSF, word2 );
		}
		else if( !strcasecmp( word1, "solvatePDB" ) )
		{
			if( block->solvatePDB )
				free(block->solvatePDB);
			block->solvatePDB = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->solvatePDB, word2 );
		}
		else if( !strcasecmp( word1, "solvatePSF" ) )
		{
			if( block->solvatePSF )
				free(block->solvatePSF);
			block->solvatePSF = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->solvatePSF, word2 );
		}
		else if( !strcasecmp( word1, "innerPatchPDB" ) )
		{
			if( block->innerPatchPDB )
				free(block->innerPatchPDB);
			block->innerPatchPDB = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->innerPatchPDB, word2 );
		}
		else if( !strcasecmp( word1, "innerPatchPSF" ) )
		{
			if( block->innerPatchPSF )
				free(block->innerPatchPSF);
			block->innerPatchPSF = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->innerPatchPSF, word2 );
		}
		else if( !strcasecmp( word1, "patchPDB" ) )
		{
			if( block->patchPDB )
				free(block->patchPDB);
			block->patchPDB = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->patchPDB, word2 );
		}
		else if( !strcasecmp( word1, "patchPSF" ) )
		{
			if( block->patchPSF )
				free(block->patchPSF);
			block->patchPSF = (char *)malloc( sizeof(char) * (1 + strlen(word2) ) );
			strcpy( block->patchPSF, word2 );
		}
		else if( !strcasecmp( word1, "fitCoupling") )
			block->fitCoupling = atof(word2);
		else if( !strcasecmp( word1, "fitThickness") )
			block->fitThickness = atof(word2);
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
	printf("\tKA:        %lf kcal/mol/A^2\n", block->KA );
	printf("\tkc:        %lf kcal/mol\n", block->kc );
	printf("Random seed: %d\n", block->random_seed );
	printf("Timestep:    %le seconds.\n", block->time_step );

	if( block->do_ld )
	{
		printf("Doing Langevin dynamics with gamma=%lf.\n", block->gamma_langevin );
	}
	
	if( block->do_ld == 0 && block->nequil > 0 )
	{
		printf("Doing Langevin dynamics only during the %d equilibration steps with gamma=%lf\n", block->nequil,	block->gamma_langevin );
	}

	if( block->lipid_mc_period > 0 )
		printf("Lipid Monte Carlo move period: %d\n", block->lipid_mc_period );
	if( block->npt_mc_period > 0 )
		printf("NPT Monte Carlo move period: %d\n", block->npt_mc_period );	
	
}




