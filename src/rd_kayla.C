#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include "parallel.h"
#include <string.h>
#include "global_boxing.h"
#include "units.h"
#include "mutil.h"
#include "rd_kayla.h"

static global_boxing *boxing = NULL;

void RD::get_tracked(surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, RD **tracked)
{
	// set up the system boxing
	double *alphas = rsurf + theSurface->nv*3;
	double Rmax = 3 * binding_radius; //some relationship to binding_radius
	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->setrall(theSurface, rsurf);
	if( !boxing )
	{
		boxing = (global_boxing *)malloc( sizeof(global_boxing) );
		boxing->setup_boxing(Rmax, theSurface->PBC_vec); //box the system
	}
	boxing->setPBC(theSurface->PBC_vec, alphas); // add periodic boundary conditions

	int ntotp = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		ntotp += allComplexes[c]->nsites;
	}
	int *to_clear = (int *)malloc( sizeof(int) * ntotp );
	int *complex_for_id = (int *)malloc( sizeof(int) * ntotp );
        int *subp_for_id = (int *)malloc( sizeof(int) * ntotp );

        int id = 0;
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
	int *nearlist = (int *)malloc( sizeof(int) * ncomplex );

	for(int p = 0; p < ncomplex; p++)
	{
		int npairs;
		int near = boxing->getNearPts(allComplexes[p]->rall+3*p, nearlist, binding_radius);
		for(int n = 0; n < near; n++)
		{
			int p2 = nearlist[n];

			if( p == p2 )
				continue;
			if( p2 < p )
				continue;
			// we have the list, now we need the particles within the binding radius that are being tracked -- kept in the tracked structure
			double dr[3] = {
				allComplexes[p]->rall[3*p+0] - allComplexes[p2]->rall[3*p2+0],
				allComplexes[p]->rall[3*p+1] - allComplexes[p2]->rall[3*p2+1],
				allComplexes[p]->rall[3*p+2] - allComplexes[p2]->rall[3*p2+2] };

			theSurface->wrapPBC( dr, alphas );
			double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
			if( r < binding_radius )
			{
				npairs++;
				if(npairs == tracked[p]->ntracked)
				{
					tracked[p]->tracked_id = (int *)realloc(tracked[p]->tracked_id,  sizeof(int) * npairs);
				}
				tracked[p]->tracked_id[npairs] = p2;
			}
		}
		int nnew = npairs - tracked[p]->ntracked;
		tracked[p]->ntracked += nnew;
	}
}

//void RD() //use the structure from get_tracked and perfrom RD, update pcomplex structure
//{

//}

//void garbage_clean() //cleanup the garbage in pcomplex after performing x # of iterations of RD
//{

//}
