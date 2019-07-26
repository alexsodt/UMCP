#include "simulation.h"
#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include "parallel.h"
#include <string.h>
#include "global_boxing.h"
#include "units.h"
#include "mutil.h"
#include "rd_kayla.h"
#include "fpr_subroutines/2D.h"
#include "random_global.h"

static global_boxing *boxing = NULL;

void RD::init( Simulation * theSimulation, double time_step_in)
{
	int ncomplex =theSimulation->ncomplex;

        // binding radius and the rates will eventually be read from a database or user input.
        // for now K debugs on these great values.
        binding_radius = 10.0;
	Dtot	       = 1e10; // for initial development.
        k_off          = 1e4; // per second.
        k_on           = 1e10; // per second times A^3 I guess?
	dt	       = time_step_in;
	Rmax 	       = binding_radius + 3 * sqrt( 6 * Dtot * dt );
        // setup tracked: allocate tracked
        tracked = (RD_tracked **)malloc(sizeof(RD_tracked)*ncomplex);
        // for each particle, set ntracked_space = 10, set ntracked=0
        for(int p = 0; p < ncomplex; p++)
        {
                tracked[p] = (RD_tracked *)malloc(sizeof(RD_tracked));
                tracked[p]->ntracked = 0;
                tracked[p]->ntracked_space = 10;
                tracked[p]->tracked_info = (RD_tracked_info *)malloc(sizeof(RD_tracked_info) * tracked[p]->ntracked_space);
                tracked[p]->tracked_new = (RD_tracked_info *)malloc(sizeof(RD_tracked_info) * tracked[p]->ntracked_space);
        }
}

void RD::get_tracked( Simulation *theSimulation  )
{
	pcomplex **allComplexes = theSimulation->allComplexes;
	int ncomplex = theSimulation->ncomplex;

	// set up the system boxing
	double *alphas = theSimulation->alpha;
	double Rmax = 3 * binding_radius; //some relationship to binding_radius
	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->setrall(theSimulation);
	if( !boxing )
	{
		boxing = (global_boxing *)malloc( sizeof(global_boxing) );
		boxing->setup_boxing(Rmax, theSimulation->PBC_vec); //box the system
	}
	boxing->setPBC(theSimulation->PBC_vec, alphas); // add periodic boundary conditions

	int ntotp = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		ntotp += allComplexes[c]->nsites;
	}
	int *to_clear = (int *)malloc( sizeof(int) * ntotp );
	int *complex_for_id = (int *)malloc( sizeof(int) * ntotp );
        int *subp_for_id = (int *)malloc( sizeof(int) * ntotp );

        int id = 0;

	//////////////// TO DO:
	// from list of particles that can do a binding reaction, only add those to this list:
	// currently we are adding all particles (suitable for initial test of lipid-lipid dimerization)

        for( int c1 = 0; c1 < ncomplex; c1++ )
	{
		for( int p = 0; p < allComplexes[c1]->nsites; p++ )
		{
			complex_for_id[id] = c1;
			subp_for_id[id] = p;
			to_clear[id] = boxing->addp(allComplexes[c1]->rall+3*p, id );
			id++;
		}
	}
	
	// find particles that are near you, keep track of the "new" tracked particles and update the "previously" tracked particles in the tracked structure
	int *nearlist = (int *)malloc( sizeof(int) * ntotp );

	// reset info flag on previously tracked data.
	for( int p = 0; p < ncomplex; p++ )
	{
		for( int n = 0; n < tracked[p]->ntracked; n++ ) 
			tracked[p]->tracked_info[n].info = 0;
	}

	for(int p = 0; p < ncomplex; p++)
	{
		/////////// TO DO: only do sites that can do binding reaction.
		for( int s = 0; s < allComplexes[p]->nsites; s++ )
		{
			int npairs = 0;
			int nnew = 0;
			/////////// TO-DO: use max Rmax from this site's reactions.
			int near = boxing->getNearPts(allComplexes[p]->rall+3*s, nearlist, Rmax);
		
			for(int np = 0; np < near; np++)
			{
				int id = nearlist[np];
				
				int s2 = subp_for_id[id];
				int p2 = complex_for_id[id];

				if( p == p2 )
					continue;
				if( p2 < p )
					continue;
	
				///////// TO-DO verify that the pair found can undergo a reaction.

				double dr[3] = {
					allComplexes[p]->rall[3*s+0] - allComplexes[p2]->rall[3*s2+0],
					allComplexes[p]->rall[3*s+1] - allComplexes[p2]->rall[3*s2+1],
					allComplexes[p]->rall[3*s+2] - allComplexes[p2]->rall[3*s2+2] };
	
				theSimulation->wrapPBC( dr, alphas );
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				if( r < binding_radius )
				{
					npairs++;
					int old = 0;
					int old_id;
					// Need to determine if newly tracked or already tracked
					for(int n = 0; n < tracked[p]->ntracked; n++)
					{
						// Does this pair already exist?
						if(tracked[p]->tracked_info[n].id == p2 && 
						   tracked[p]->tracked_info[n].sid2 == s2 &&
						   tracked[p]->tracked_info[n].sid1 == s )
						{
							old = 1;
							old_id = n;
						}
	
					}
					//If already tracked need to update prev_sep, prevnorm, curr_sep, etc.
					if(old)
					{
						tracked[p]->tracked_info[old_id].info = 1;
						tracked[p]->tracked_info[old_id].prev_sep = tracked[p]->tracked_info[old_id].curr_sep;
						tracked[p]->tracked_info[old_id].prev_norm = tracked[p]->tracked_info[old_id].curr_norm;
						tracked[p]->tracked_info[old_id].curr_sep = r;
						tracked[p]->tracked_info[old_id].curr_norm = 1.0;
					}
					// If newly tracked need to add to list
					else
					{
						if(nnew == tracked[p]->ntracked_space )
						{
							tracked[p]->ntracked_space *= 2;
							tracked[p]->tracked_new = (RD_tracked_info *)realloc(tracked[p]->tracked_new,  sizeof(RD_tracked_info) * tracked[p]->ntracked_space);
							tracked[p]->tracked_info = (RD_tracked_info *)realloc(tracked[p]->tracked_info,  sizeof(RD_tracked_info) * tracked[p]->ntracked_space);
						}

						//////// TO-DO: Store reaction information here to be retrieved later (in do_rd)
						tracked[p]->tracked_new[nnew].id = p2;
						tracked[p]->tracked_new[nnew].sid1 = s;
						tracked[p]->tracked_new[nnew].sid2 = s2;
						tracked[p]->tracked_new[nnew].curr_sep = r;
						tracked[p]->tracked_new[nnew].curr_norm = 1.0;
						
						nnew++;
					}
				}
			}
				
			while(npairs >= tracked[p]->ntracked_space )
			{
				tracked[p]->ntracked_space *= 2;
				tracked[p]->tracked_new = (RD_tracked_info *)realloc(tracked[p]->tracked_new,  sizeof(RD_tracked_info) * tracked[p]->ntracked_space);
				tracked[p]->tracked_info = (RD_tracked_info *)realloc(tracked[p]->tracked_info,  sizeof(RD_tracked_info) * tracked[p]->ntracked_space);
			}

			int count = nnew;

			for(int n = 0; n < tracked[p]->ntracked; n++)
			{
				if(tracked[p]->tracked_info[n].info)
				{
					tracked[p]->tracked_new[count].id = tracked[p]->tracked_info[n].id;
					tracked[p]->tracked_new[count].prev_sep = tracked[p]->tracked_info[n].prev_sep;
					tracked[p]->tracked_new[count].prev_norm = tracked[p]->tracked_info[n].prev_norm;
					tracked[p]->tracked_new[count].curr_sep = tracked[p]->tracked_info[n].curr_sep;
					tracked[p]->tracked_new[count].prev_norm = tracked[p]->tracked_info[n].curr_norm;
					tracked[p]->tracked_new[count].info = 0;	
					count++;
				}
			}

			tracked[p]->ntracked = count;

			// swap new <-> old
			RD_tracked_info *temp = tracked[p]->tracked_new;
			tracked[p]->tracked_info = tracked[p]->tracked_new;
			tracked[p]->tracked_new = temp;
		}
	}
	
	boxing->clearBoxing(to_clear, ntotp);
	free(to_clear);
	free(complex_for_id);
	free(subp_for_id);
	free(nearlist);
}

void RD::do_rd( Simulation *theSimulation )
{
	double prob;
	int ncomplex = theSimulation->ncomplex;

	get_tracked( theSimulation  );

	for(int p = 0; p < ncomplex; p++)
        {
		for(int n = 0; n < tracked[p]->ntracked; n++)
		{
			double rn = gsl_rng_uniform(rng_x);
		
			// TO DO: retrieve stored reaction information from tracked_info.			
			prob = get_2D_2D_rxn_prob(tracked[p]->tracked_info[n].curr_sep, k_on, binding_radius, Dtot, dt, Rmax);

			if(prob > rn)
			{
				// binding reaction between these two complexes.				
			}
		}
	}
}
