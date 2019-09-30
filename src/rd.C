#include "simulation.h"
#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include "parallel.h"
#include <string.h>
#include "global_boxing.h"
#include "units.h"
#include "mutil.h"
#include "rd.h"
#include "fpr_subroutines/2D.h"
#include "random_global.h"
#include "util.h"
#include "uv_map.h"

static global_boxing *boxing = NULL;

void RD::init( Simulation * theSimulation, double time_step_in, parameterBlock *block )
{
	params = block;

	nreactants = 0;
	nreactantsSpace = 10;
	reactants = (ReactantType *)malloc( sizeof(ReactantType) * nreactantsSpace );
	rxn_lookup_table = NULL;

	nreactionsSpace = 10;
	nreactions = 0;

	allReactions = (ReactionInformation *)malloc( sizeof(ReactionInformation) * nreactionsSpace );

	if( block->rxnDiffusionInfoName )
		parseRDFile( block->rxnDiffusionInfoName );	

        // binding radius and the rates will eventually be read from a database or user input.
        // for now K debugs on these great values.
	dt	       	       = time_step_in;

	nsites_tracked = 0;
	nsites_tracked_space = theSimulation->ncomplex;

	tracked = (RD_tracked **)malloc( sizeof(RD_tracked *) * nsites_tracked_space );

	for( int p = 0; p < theSimulation->ncomplex; p++ )
	{
		// if there is a reaction associated with this site then we need to track it. This will be called repeatedly as sites are created by reactions.
		for( int s = 0; s < theSimulation->allComplexes[p]->nsites; s++ )
			registerSite( theSimulation->allComplexes[p], p, s );
	}

	double max_binding_radius = 0;

	double max_diff_c = 0;

	for( int t = 0; t < nsites_tracked; t++ )
	{
		tracked[t]->reacted = 0;
		int pid = tracked[t]->pid;
		int sid = tracked[t]->sid;

		if( theSimulation->allComplexes[pid]->DC[sid] > max_diff_c )
			max_diff_c = theSimulation->allComplexes[pid]->DC[sid];
	}

	for( int r = 0; r < nreactions; r++ )
	{
		if( allReactions[r].binding_radius > max_binding_radius )
			max_binding_radius = allReactions[r].binding_radius;
	}

	max_Rmax = max_binding_radius + 3 * sqrt( 2 * max_diff_c * dt ); 
}

void RD::get_tracked( Simulation *theSimulation  )
{
	pcomplex **allComplexes = theSimulation->allComplexes;
	int ncomplex = theSimulation->ncomplex;


	// set up the system boxing
	double *alphas = theSimulation->alpha;
        int id = 0;

	// reset info flag on previously tracked data.
	for( int t = 0; t < nsites_tracked; t++ )
	{
		for( int n = 0; n < tracked[t]->ntracked; n++ ) 
			tracked[t]->tracked_info[n].info = 0;
	}

	double n_to_track = 0;

	for( int t = 0; t < nsites_tracked; t++ )
	{
		int p = tracked[t]->pid;
		int s = tracked[t]->sid;

		if( allComplexes[p]->disabled ) continue;

		int npairs = 0;
		int nnew = 0;
		/////////// TO-DO: use max Rmax from this site's reactions.

		int near = boxing->getNearPts(allComplexes[p]->rall+3*s, nearlist, max_Rmax );
	
		for(int np = 0; np < near; np++)
		{
			int id = nearlist[np];
			
			int s2 = subp_for_id[id];
			int p2 = complex_for_id[id];
			if( allComplexes[p2]->disabled ) continue;

			if( p == p2 )
				continue;
			if( p2 < p )
				continue;

			int rxn = rxnLookup( allComplexes[p]->stype[s], allComplexes[p2]->stype[s2] );
		
			if( rxn == -1 ) continue;

			double binding_radius = allReactions[rxn].binding_radius;

			double D1 = allComplexes[p]->DC[s];
			double D2 = allComplexes[p2]->DC[s2];
			if( s <  allComplexes[p]->nattach ) D1 *= (2.0/3.0);
			if( s2 <  allComplexes[p2]->nattach ) D2 *= (2.0/3.0);
			double Rmax = binding_radius + 3 * sqrt( (D1+D2) * dt );

			double dr[3] = {
				allComplexes[p]->rall[3*s+0] - allComplexes[p2]->rall[3*s2+0],
				allComplexes[p]->rall[3*s+1] - allComplexes[p2]->rall[3*s2+1],
				allComplexes[p]->rall[3*s+2] - allComplexes[p2]->rall[3*s2+2] };

			theSimulation->wrapPBC( dr, alphas );
			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

			if( r < Rmax)
			{
				n_to_track++;
				npairs++;
				int old = 0;
				int old_id;
				// Need to determine if newly tracked or already tracked
				for(int n = 0; n < tracked[t]->ntracked; n++)
				{
					// Does this pair already exist?
					if(tracked[t]->tracked_info[n].id == p2 && 
					   tracked[t]->tracked_info[n].sid2 == s2 &&
					   tracked[t]->tracked_info[n].sid1 == s )
					{
						old = 1;
						old_id = n;
					}

				}
				//If already tracked need to update prev_sep, prevnorm, curr_sep, etc.
				if(old)
				{
					tracked[t]->tracked_info[old_id].info = 1;
					tracked[t]->tracked_info[old_id].prevsep = tracked[t]->tracked_info[old_id].curr_sep;
					tracked[t]->tracked_info[old_id].prevnorm = tracked[t]->tracked_info[old_id].currnorm;
					tracked[t]->tracked_info[old_id].curr_sep = r;
					tracked[t]->tracked_info[old_id].currnorm = 1.0;
				}
				// If newly tracked need to add to list
				else
				{
					if(nnew == tracked[t]->ntracked_space )
					{
						tracked[t]->ntracked_space *= 2;
						tracked[t]->tracked_new = (RD_tracked_info *)realloc(tracked[t]->tracked_new,  sizeof(RD_tracked_info) * tracked[t]->ntracked_space);
						tracked[t]->tracked_info = (RD_tracked_info *)realloc(tracked[t]->tracked_info,  sizeof(RD_tracked_info) * tracked[t]->ntracked_space);
					}

					//////// TO-DO: Store reaction information here to be retrieved later (in do_rd)
					tracked[t]->tracked_new[nnew].id = p2;
					tracked[t]->tracked_new[nnew].sid1 = s;
					tracked[t]->tracked_new[nnew].sid2 = s2;
					tracked[t]->tracked_new[nnew].curr_sep = r;
					tracked[t]->tracked_new[nnew].currnorm = 1.0;
					tracked[t]->tracked_new[nnew].prevnorm = 1.0;
					tracked[t]->tracked_new[nnew].prevsep= -1;
					tracked[t]->tracked_new[nnew].ps_prev = 1.0;
			
					
					nnew++;
				}
			}
		}
			
		while(npairs >= tracked[t]->ntracked_space )
		{
			tracked[t]->ntracked_space *= 2;
			tracked[t]->tracked_new = (RD_tracked_info *)realloc(tracked[t]->tracked_new,  sizeof(RD_tracked_info) * tracked[t]->ntracked_space);
			tracked[t]->tracked_info = (RD_tracked_info *)realloc(tracked[t]->tracked_info,  sizeof(RD_tracked_info) * tracked[t]->ntracked_space);
		}

		int count = nnew;

		for(int n = 0; n < tracked[t]->ntracked; n++)
		{
			if(tracked[t]->tracked_info[n].info)
			{
				tracked[t]->tracked_new[count].id = tracked[t]->tracked_info[n].id;
				tracked[t]->tracked_new[count].prevsep = tracked[t]->tracked_info[n].prevsep;
				tracked[t]->tracked_new[count].ps_prev = tracked[t]->tracked_info[n].ps_prev;
				tracked[t]->tracked_new[count].prevnorm = tracked[t]->tracked_info[n].prevnorm;
				tracked[t]->tracked_new[count].curr_sep = tracked[t]->tracked_info[n].curr_sep;
				tracked[t]->tracked_new[count].prevnorm = tracked[t]->tracked_info[n].currnorm;
				tracked[t]->tracked_new[count].sid1 = tracked[t]->tracked_info[n].sid1;
				tracked[t]->tracked_new[count].sid2 = tracked[t]->tracked_info[n].sid2;
				tracked[t]->tracked_new[count].info = 0;	
				count++;
			}
		}

		tracked[t]->ntracked = count;

		// swap new <-> old
		RD_tracked_info *temp = tracked[t]->tracked_info;
		tracked[t]->tracked_info = tracked[t]->tracked_new;
		tracked[t]->tracked_new = temp;
	}

}

void RD::do_rd( Simulation *theSimulation )
{
	/*
		theSimulation->allComplexes can get re-alloc'd here. don't save it in a local variable.

	*/
	double prob;
	int ncomplex = theSimulation->ncomplex;
	
	static int *bin_counts = NULL;
	static double *p_sampled = NULL;
	static int nbins = 100;
	static int ncount=0;
	static int delay = 0;
	static int delay_start = 0;
	static double use_binding_radius = 5.0;
	static double use_Rmax = 18;
	static int nreact = 0, ndissoc = 0;

	static double p_no_reaction = 1.0;
	static double sum_prob = 0;

	if( !bin_counts )
	{
		bin_counts = (int *)malloc( sizeof(int) * nbins );
		memset( bin_counts, 0, sizeof(int) * nbins );
		p_sampled = (double *)malloc( sizeof(double) * nbins );
		memset( p_sampled, 0, sizeof(double) * nbins );
	}

	// box_reactants MUST be called before get_tracked, sets up boxing across function calls.
	get_tracked( theSimulation  );

	double nfound = 0;
	
	for( int p = 0; p < theSimulation->ncomplex; p++ ) 	
	for( int s = 0; s < theSimulation->allComplexes[p]->nsites; s++ )
		theSimulation->allComplexes[p]->rd_timestep_disabled[s] = 0;

	for( int t = 0; t < nsites_tracked; t++ )
	{
		tracked[t]->reacted = 0;
		int p = tracked[t]->pid;
		int s = tracked[t]->sid;
		
			
		if( theSimulation->allComplexes[p]->disabled ) continue;
	
		for( int n = 0; n < tracked[t]->ntracked; n++ )
		{
			int p2 = tracked[t]->tracked_info[n].id;	
			int s2 = tracked[t]->tracked_info[n].sid2;

			if( theSimulation->allComplexes[p2]->disabled ) continue;
			int rxn = rxnLookup( theSimulation->allComplexes[p]->stype[s], theSimulation->allComplexes[p2]->stype[s2] );
			if( rxn == -1 ) continue; // should never be -1..
	
			double prevnorm = tracked[t]->tracked_info[n].prevnorm;
			double prevsep  = tracked[t]->tracked_info[n].prevsep;
			double ps_prev  = tracked[t]->tracked_info[n].ps_prev;

			double k_on = allReactions[rxn].k_on;
			double k_off = allReactions[rxn].k_off;
			double binding_radius = allReactions[rxn].binding_radius;
			double D1 = theSimulation->allComplexes[p]->DC[s];
			double D2 = theSimulation->allComplexes[p]->DC[s2];
			if( s <  theSimulation->allComplexes[p]->nattach ) D1 *= (2.0/3.0);
			if( s2 <  theSimulation->allComplexes[p2]->nattach ) D2 *= (2.0/3.0);
			double Rmax = binding_radius + 3 * sqrt( (D1+D2) * dt );
//			double rn = gsl_rng_uniform(rng_x);
			double rn = rand() / (double)RAND_MAX;
			// TO DO: retrieve stored reaction information from tracked_info.			

			double p0_ratio = 1.0;

#ifndef DISABLE_RD
			double pre_prob = get_2D_2D_rxn_prob(tracked[t]->tracked_info[n].curr_sep, k_on*2*binding_radius, binding_radius, D1+D2, dt, Rmax,
				prevsep, ps_prev, &p0_ratio);
#else
			prevsep=0;
			double pre_prob = get_2D_2D_rxn_prob(tracked[t]->tracked_info[n].curr_sep, k_on*2*binding_radius, binding_radius, D1+D2, dt, Rmax,
				prevsep, ps_prev, &p0_ratio);
//			double pre_prob = 0;
#endif	


			double currnorm = tracked[t]->tracked_info[n].prevnorm * p0_ratio;
			double prob = pre_prob * currnorm;


			tracked[t]->tracked_info[n].ps_prev = 1.0 - pre_prob * currnorm;
			tracked[t]->tracked_info[n].prevnorm = currnorm;

			use_Rmax = Rmax;
			use_binding_radius = binding_radius;


			nfound++;
			if( delay > delay_start )
			{
				int rb = nbins*(tracked[t]->tracked_info[n].curr_sep-binding_radius)/(Rmax-binding_radius);	
				if( rb >= nbins ) rb = nbins-1;
				if( rb < 0 )
				{	
					 rb = 0;
				}
				bin_counts[rb] += 1;
				p_sampled[rb] += prob;
			}
#ifndef DISABLE_RD
			if(prob > rn )
			{
				tracked[t]->reacted = 1;
				sum_prob = 0;
				p_no_reaction = 1.0;
//				printf("Reacted! p: %le\n", prob );
				// binding reaction between these two complexes... 
				// possible outcomes are to add to a previous complex (likely case with Actin polymerization) or create a new one.	
				// for now: create new complex.
	
				pcomplex *product = loadComplex( allReactions[rxn].productName );
				product->loadParams(params);			
				struct surface_record *sRec = theSimulation->fetch( theSimulation->allComplexes[p]->sid[s] );
				product->init( sRec->theSurface, sRec->r, theSimulation->allComplexes[p]->fs[s], theSimulation->allComplexes[p]->puv[2*s+0], theSimulation->allComplexes[p]->puv[2*s+1] );
				product->copyParentParameters( theSimulation->allComplexes[p] );
				for( int s = 0; s < product->nattach; s++ )
					product->sid[s] = sRec->id;
				int padd = theSimulation->AddComplex( product );
				for( int s = 0; s < theSimulation->allComplexes[padd]->nsites; s++ )
					registerSite( theSimulation->allComplexes[padd], padd, s );
	 
				// this stops the complex from being propagated, simulated, etc. leaves it for garbage collection later.
				theSimulation->RemoveComplexDelayed(p);
				theSimulation->RemoveComplexDelayed(p2);			

				product->refresh(theSimulation);
				
				for( int s = 0; s < product->nsites; s++ )
					product->rd_timestep_disabled[s] = 1;

				nreact+=1;
			
				break;
			}
			else
			{
				p_no_reaction *= (1.0-prob); 
				sum_prob += prob;
			}
#endif
		}
	}

//#define DEBUG_RD

#ifdef DEBUG_RD
	printf("p_no_reaction: %le sum_prob: %le\n", p_no_reaction, sum_prob );

	if( delay > delay_start )
		ncount++;
	delay++;
	for( int b = 0; b < nbins; b++ )
		printf("%le %le %le\n",
			use_binding_radius+(b+0.5)*(use_Rmax-use_binding_radius)/nbins, bin_counts[b]/(double)ncount, p_sampled[b]/ncount);
#endif
//	printf("After rd, ncomplex: %d\n", theSimulation->ncomplex );

	// dissociation.

	for( int p = 0; p < theSimulation->ncomplex; p++ )
	for( int r = 0; r < nreactions; r++ )
	{   
		int reacted = 0;
		for( int s = 0; s < theSimulation->allComplexes[p]->nsites; s++ )
			if( theSimulation->allComplexes[p]->rd_timestep_disabled[s] ) reacted = 1;
		if( reacted ) continue;
		if( !strcasecmp( theSimulation->allComplexes[p]->complex_name, allReactions[r].productName ) ) 
		{   
			double pr = allReactions[r].k_off * dt; 
			double rn = gsl_rng_uniform(rng_x);

#ifndef DISABLE_RD
			if( rn < pr ) 
			{
				// creates the reactants.

				double binding_radius = allReactions[r].binding_radius;

				// HACK: now only working for single-site attachment.
				int s = 0;
	
//				printf("Dissociated! p: %le\n", pr );
				// binding reaction between these two complexes... 
				// possible outcomes are to add to a previous complex (likely case with Actin polymerization) or create a new one.	
				// for now: create new complex.
	
				pcomplex *reactant1 = loadComplex( allReactions[r].pcomplex_name1 );
				reactant1->loadParams(params);			
				pcomplex *reactant2 = loadComplex( allReactions[r].pcomplex_name2 );
				reactant2->loadParams(params);			
	
				struct surface_record *sRec = theSimulation->fetch( theSimulation->allComplexes[p]->sid[s] );

				int trial_f = theSimulation->allComplexes[p]->fs[s];
				double trial_u = theSimulation->allComplexes[p]->puv[2*s+0];
				double trial_v = theSimulation->allComplexes[p]->puv[2*s+1];

				int rd_blocked = 0;

				reactant1->init( sRec->theSurface, sRec->r, theSimulation->allComplexes[p]->fs[s], theSimulation->allComplexes[p]->puv[2*s+0], theSimulation->allComplexes[p]->puv[2*s+1] );
				reactant2->init( sRec->theSurface, sRec->r, theSimulation->allComplexes[p]->fs[s], theSimulation->allComplexes[p]->puv[2*s+0], theSimulation->allComplexes[p]->puv[2*s+1] );

				reactant1->copyParentParameters( theSimulation->allComplexes[p] );
				reactant2->copyParentParameters( theSimulation->allComplexes[p] );

				for( int s = 0; s < reactant1->nattach; s++ )
					reactant1->sid[s] = sRec->id;
				int padd1 = theSimulation->AddComplex( reactant1 );
				for( int s = 0; s < theSimulation->allComplexes[padd1]->nsites; s++ )
					registerSite( theSimulation->allComplexes[padd1], padd1, s );
				
				// this stops the complex from being propagated, simulated, etc. leaves it for garbage collection later.
				theSimulation->RemoveComplexDelayed(p);
					
				// Eventually we will place the reactants based on Boltzmann probabilities.
				// for now, just don't block RD.

				theSimulation->allComplexes[padd1]->refresh(theSimulation);
				do {
					rd_blocked = check_RD_blocked( theSimulation, padd1, 0 ); 	

					if( rd_blocked )
						theSimulation->allComplexes[padd1]->propagate_surface_q( theSimulation, dt );

				} while( rd_blocked );
				
				for( int s = 0; s < reactant2->nattach; s++ )
					reactant2->sid[s] = sRec->id;
				int padd2 = theSimulation->AddComplex( reactant2 );
				for( int s = 0; s < theSimulation->allComplexes[padd2]->nsites; s++ )
					registerSite( theSimulation->allComplexes[padd2], padd2, s );

				double duv[2];
				randomDirection( sRec->theSurface, sRec->r, reactant1->fs[0], reactant1->puv[0], reactant1->puv[1], binding_radius*1.03, duv );

				int f_1 = reactant1->fs[0];
				int nf = f_1;
				double uv1[2] = { reactant1->puv[0], reactant1->puv[1] };
				do {
					f_1 = nf;
					nf = sRec->theSurface->nextFace( f_1, uv1+0, uv1+1, duv+0, duv+1, sRec->r ); 
				} while( nf != f_1 );


				theSimulation->allComplexes[padd2]->fs[0] = f_1;
				theSimulation->allComplexes[padd2]->puv[0] = uv1[0];
				theSimulation->allComplexes[padd2]->puv[1] = uv1[1];
				theSimulation->allComplexes[padd2]->refresh(theSimulation);
				do {
					rd_blocked = check_RD_blocked( theSimulation, padd2, 0 ); 	

					double dr[3] = { 
							theSimulation->allComplexes[padd2]->rall[0] - theSimulation->allComplexes[padd1]->rall[0],
							theSimulation->allComplexes[padd2]->rall[1] - theSimulation->allComplexes[padd1]->rall[1],
							theSimulation->allComplexes[padd2]->rall[2] - theSimulation->allComplexes[padd1]->rall[2] };

					double lr =normalize(dr);

					if( lr < allReactions[r].binding_radius )
						rd_blocked = 1;
					
					if( rd_blocked )
						theSimulation->allComplexes[padd2]->propagate_surface_q( theSimulation, dt );

				} while( rd_blocked );
				
				// TO-DO: needs to disable only the reaction site.
				for( int s = 0; s < reactant1->nsites; s++ )
					reactant1->rd_timestep_disabled[s] = 1;
				for( int s = 0; s < reactant2->nsites; s++ )
					reactant2->rd_timestep_disabled[s] = 1;
				ndissoc+=1;

				double dr[3] = { 
					reactant2->rall[0] - reactant1->rall[0],
					reactant2->rall[1] - reactant1->rall[1],
					reactant2->rall[2] - reactant1->rall[2] };
				double lr = normalize(dr);
//				printf("PLACED AT %le\n", lr );
			}
#endif
		}   
	}   

#define RD_GM1_HACK // REMOVE THIS CODE FOR DISTRIBUTION

#ifdef RD_GM1_HACK
	int nMonomer = 0;
	int nDimer = 0;
	static int counter = 0;
	for( int c = 0; c < theSimulation->ncomplex; c++ )
	{
		if( theSimulation->allComplexes[c]->disabled ) continue;

		if( !strcasecmp( theSimulation->allComplexes[c]->complex_name, "simpleLipid" ) )
			nMonomer++;	
		if( !strcasecmp( theSimulation->allComplexes[c]->complex_name, "simpleDimer" ) )
			nDimer++;	
	}	
		
	if( counter % 100 == 0 )
		printf("RDOUT %le %d %d nreact: %d ndissoc: %d\n", theSimulation->current_time, nMonomer, nDimer, nreact, ndissoc );
	counter++;
#endif
	
}

void RD::parseRDFile( const char *fileName )
{
	FILE *theFile = fopen(fileName,"r");

	if( !theFile )
	{
		printf("Couldn't open reaction-diffusion input file '%s'.\n", fileName );
		exit(1);
	}

	char *buffer = (char *)malloc( sizeof(char) * 100000 );

	while( !feof(theFile) )
	{
		getLine( theFile, buffer );

		if( feof(theFile) ) break;

		// Syntax: pcomplex_name1 site_type1 pcomplex_name2 site_type2 k_on k_off radius
		const char *syntax = "pcomplex_name1 site_type1 pcomplex_name2 site_type2 k_on k_off radius productName";
		int len = strlen(buffer);

		char *parse = buffer;

		while( *parse == ' ' || *parse == '\t' ) parse += 1;

		if( *parse == '#' ) continue;


		char *pn1 = (char *)malloc( sizeof(char) * ( 1 + len ) );
		char *pn2 = (char *)malloc( sizeof(char) * ( 1 + len ) );
		char *st1 = (char *)malloc( sizeof(char) * ( 1 + len ) );
		char *st2 = (char *)malloc( sizeof(char) * ( 1 + len ) );
		char *productName = (char *)malloc( sizeof(char) * ( 1 + len ) );		

		double k_on, k_off, binding_radius;
		int nr = sscanf( parse, "%s %s %s %s %lf %lf %lf %s",
			pn1, st1, pn2, st2, &k_on, &k_off, &binding_radius, productName );

		if( nr != 8 )
		{
			printf("Failed to read eight fields from '%s'.\n", buffer );
			printf("Syntax: %s\n", syntax );
			exit(1); 
		}

		int site_code_1 = SITE_TYPE_ANY;
		if( !strcasecmp( st1, "any" ) )
			site_code_1 = SITE_TYPE_ANY;
		else  // put other cases you can dream up and use... SITE_TYPE_SH3, etc.
			site_code_1 = atoi(st1);	
		
		int site_code_2 = SITE_TYPE_ANY;
		if( !strcasecmp( st2, "any" ) )
			site_code_2 = SITE_TYPE_ANY;
		else  // put other cases you can dream up and use... SITE_TYPE_SH3, etc.
			site_code_2 = atoi(st2);	
	
		registerReaction( site_code_1, pn1, site_code_2, pn2, k_on, k_off, binding_radius, productName );

		free(productName);	
		free(pn1);
		free(pn2);
		free(st1);
		free(st2);
	}
}

void RD::registerSite( pcomplex *theComplex, int pid, int sid )
{
	// does this site have any reactions? has it been entered yet?

	// TO-DO: read this from theComplex.

	int site_type1 = SITE_TYPE_ANY;
	int entry = -1;

	for( int ts = 0; ts < nreactants; ts++ )
	{
		if( reactants[ts].stype == site_type1 &&
		    !strcasecmp( reactants[ts].complex_name, theComplex->complex_name) ) 
		{
			entry = ts;	
		}
	}

	double max_Rmax = 0;

	if( entry == -1 )
	{
		// not entered, does it have any reactions?

		int got_reaction = 0;

		for( int r = 0; r < nreactions; r++ )
		{
			int l_got_reaction = 0;
			if( allReactions[r].site1 == site_type1 && !strcasecmp( allReactions[r].pcomplex_name1, theComplex->complex_name) )
				l_got_reaction = 1;			
			if( allReactions[r].site2 == site_type1 && !strcasecmp( allReactions[r].pcomplex_name2, theComplex->complex_name) )
				l_got_reaction = 1;			
			if( l_got_reaction ) 
				got_reaction = 1;
		}

		if( got_reaction )
		{
			// create an entry in the reactant list.

			if( nreactants == nreactantsSpace )
			{
				nreactantsSpace *= 2;
				reactants = (ReactantType *)realloc( reactants, sizeof(ReactantType) * nreactantsSpace );
			}

			reactants[nreactants].stype = site_type1;
			reactants[nreactants].complex_name = (char *)malloc( sizeof(char) * (1+strlen(theComplex->complex_name) ) );
			strcpy( reactants[nreactants].complex_name, theComplex->complex_name);

			entry = nreactants;

			nreactants++;

			// new entry, rebuild the rxn_lookup table.

			int *new_rxn_table = (int *)malloc( sizeof(int) * nreactants * nreactants );
			for( int r1 = 0; r1 < nreactants; r1++ )
			for( int r2 = 0; r2 < nreactants; r2++ )
			{		
				new_rxn_table[r1*nreactants+r2] = -1;

				if( r1 < nreactants-1 && r2 < nreactants-1 )
					new_rxn_table[r1*nreactants+r2] = rxn_lookup_table[r1*(nreactants-1)+r2]; 
			} 
			if( rxn_lookup_table ) free(rxn_lookup_table);
			rxn_lookup_table = new_rxn_table;

			for( int r = 0; r < nreactions; r++ )
			{
				if( allReactions[r].site1 == site_type1 &&
				  !strcasecmp(  allReactions[r].pcomplex_name1, theComplex->complex_name ) )
				{
					for( int r2 = 0; r2 < nreactants; r2++ )
					{
						if( allReactions[r].site2 == reactants[r2].stype &&
							  !strcasecmp(  allReactions[r].pcomplex_name2, reactants[r2].complex_name ) )
						{
							rxn_lookup_table[entry*nreactants+r2] = r;
							rxn_lookup_table[r2*nreactants+entry] = r;
						}
					}
				} 
				
				if( allReactions[r].site2 == site_type1 &&
				  !strcasecmp(  allReactions[r].pcomplex_name2, theComplex->complex_name ) )
				{
					for( int r2 = 0; r2 < nreactants; r2++ )
					{
						if( allReactions[r].site1 == reactants[r2].stype &&
							  !strcasecmp(  allReactions[r].pcomplex_name1, reactants[r2].complex_name ) )
						{
							rxn_lookup_table[entry*nreactants+r2] = r;
							rxn_lookup_table[r2*nreactants+entry] = r;
						}
					}
				} 
			}
		}
	}

	theComplex->stype[sid] = -1;

	if( entry >= 0 )
	{
		theComplex->stype[sid] = entry;

		if( nsites_tracked == nsites_tracked_space )
		{
			nsites_tracked_space *= 2;
			tracked = (RD_tracked **)realloc( tracked, sizeof(RD_tracked *) * nsites_tracked_space );
		}	

		tracked[nsites_tracked] = (RD_tracked *)malloc( sizeof(RD_tracked) );

		tracked[nsites_tracked]->ntracked = 0;
		tracked[nsites_tracked]->ntracked_space = 10;
		tracked[nsites_tracked]->pid = pid; 
		tracked[nsites_tracked]->sid = sid; 
                tracked[nsites_tracked]->tracked_info = (RD_tracked_info *)malloc(sizeof(RD_tracked_info) * tracked[nsites_tracked]->ntracked_space);
                tracked[nsites_tracked]->tracked_new = (RD_tracked_info *)malloc(sizeof(RD_tracked_info) * tracked[nsites_tracked]->ntracked_space);

		nsites_tracked++;

		theComplex->watch();
	}
}
 

int RD::rxnLookup( int code1, int code2 )
{
	if( code1 >=0 && code2 >= 0 ) 
	{
		int entry = rxn_lookup_table[code1*nreactants+code2];
 	
		return entry;
	}
	else
		return -1;
}

void RD::registerReaction( int site_type1, const char *pcomplex_name1,
			   int site_type2, const char *pcomplex_name2,
				double k_on, double k_off, double binding_radius, const char *productInstructions )
{
	if( nreactions == nreactionsSpace )
	{
		nreactionsSpace *= 2;
		allReactions = (ReactionInformation *)realloc( allReactions, sizeof(ReactionInformation) * nreactionsSpace );
	}		

	allReactions[nreactions].k_off = k_off;
	allReactions[nreactions].k_on = k_on;
	allReactions[nreactions].binding_radius = binding_radius;
	allReactions[nreactions].Dtot = 0; // these depend on simulation setup, specific particle identity
	allReactions[nreactions].Rmax = 0;
	allReactions[nreactions].site1 = site_type1;
	allReactions[nreactions].site2 = site_type2;

	allReactions[nreactions].pcomplex_name1 = (char *)malloc( sizeof(char) * (1+strlen(pcomplex_name1) ) );
	strcpy( allReactions[nreactions].pcomplex_name1, pcomplex_name1);
	allReactions[nreactions].pcomplex_name2 = (char *)malloc( sizeof(char) * (1+strlen(pcomplex_name2) ) );
	strcpy( allReactions[nreactions].pcomplex_name2, pcomplex_name2);
	allReactions[nreactions].productName = (char *)malloc( sizeof(char) * (1+strlen(productInstructions) ) );
	strcpy( allReactions[nreactions].productName, productInstructions);

	nreactions++;
}


void RD::box_reactants( Simulation *theSimulation )
{
	pcomplex **allComplexes = theSimulation->allComplexes;
	int ncomplex = theSimulation->ncomplex;
	// set up the system boxing
	double *alphas = theSimulation->alpha;
	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->setrall(theSimulation);
	if( !boxing )
	{
		boxing = (global_boxing *)malloc( sizeof(global_boxing) );
		boxing->setup_boxing(max_Rmax, theSimulation->PBC_vec); //box the system
	}

	ntotp = nsites_tracked;

	boxing->setPBC(theSimulation->PBC_vec, alphas); // add periodic boundary conditions
	to_clear = (int *)malloc( sizeof(int) * ntotp );
	complex_for_id = (int *)malloc( sizeof(int) * ntotp );
        subp_for_id = (int *)malloc( sizeof(int) * ntotp );
        int id = 0;
	for( int t = 0; t < nsites_tracked; t++ )
	{
		int p = tracked[t]->pid;
		int s = tracked[t]->sid;
			
		complex_for_id[t] = p;
		subp_for_id[t] = s;
		to_clear[t] = boxing->addp(allComplexes[p]->rall+3*s, t );
	} 
	
	nearlist = (int *)malloc( sizeof(int) * ntotp );
}

void RD::unbox_reactants( void )
{
	boxing->clearBoxing(to_clear, ntotp);
	free(to_clear);
	free(complex_for_id);
	free(subp_for_id);
	free(nearlist);

	to_clear = NULL;
	complex_for_id = NULL;
	subp_for_id = NULL;
	nearlist = NULL;
	ntotp = 0;
}

int RD::check_RD_blocked( Simulation *theSimulation, int p, int s)
{
	double *alphas = theSimulation->alpha;
	int near = boxing->getNearPts(theSimulation->allComplexes[p]->rall+3*s, nearlist, max_binding_radius );


	for(int np = 0; np < near; np++)
	{
		int id = nearlist[np];
		
		int s2 = subp_for_id[id];
		int p2 = complex_for_id[id];

		if( theSimulation->allComplexes[p2]->disabled ) continue;

		if( p == p2 )
			continue;
//		if( p2 < p )
//			continue;

		int rxn = rxnLookup( theSimulation->allComplexes[p]->stype[s], theSimulation->allComplexes[p2]->stype[s2] );
	
		if( rxn == -1 ) continue;

		double binding_radius = allReactions[rxn].binding_radius;

		double dr[3] = {
			theSimulation->allComplexes[p]->rall[3*s+0] - theSimulation->allComplexes[p2]->rall[3*s2+0],
			theSimulation->allComplexes[p]->rall[3*s+1] - theSimulation->allComplexes[p2]->rall[3*s2+1],
			theSimulation->allComplexes[p]->rall[3*s+2] - theSimulation->allComplexes[p2]->rall[3*s2+2] };

		theSimulation->wrapPBC( dr, alphas );

		double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

		if( r < binding_radius )
		{
//			printf("%d and %d blocked, r: %.12le\n", p, p2, r ); 
			return 1;	
		}
	}	

	return 0;
}


