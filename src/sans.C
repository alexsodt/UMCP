
#include "interp.h"
#include "sans.h"
#include "spline.h"
#include "util.h"
#include "input.h"
#include <math.h>
#include "parallel.h"
#include <string.h>
void loadBetaZ( const char *fileName, int load_to_spline );

static double MINZ = -75;
static double MAXZ =  75;

void surface::processSANS( parameterBlock *block )
{
	if( block->s_q )
	{
		loadBetaZ( block->betazFile, SANS_SPLINE );
	}
}

// a measurement of Sq is added into Sq_inst (it is not zero'd).

void surface::sample_B_hist( double *rmesh, double *B_hist, 
			int sample_type, 
			int nsamples, double maxr, int nbins)
{

	int samples_per_proc = ceil( nsamples/ (double)par_info.nprocs );

#ifdef PARALLEL
	double *local_sampler = (double *)malloc( sizeof(double) * nbins );
	memset( local_sampler, 0, sizeof(double) * nbins );
#endif

	if( sample_type == SANS_SAMPLE_NRM )
	{	
		// this is sampling of the full scattering signal using a laterally averaged spline of beta(z)
		// sampled stochastically. "nsamples" is relevant here.

		// we want to draw points on the surface randomly so we need to know the area of different faces right now:
//		updateFaceInfoForRandom();

		for( int smpl = 0; smpl < samples_per_proc; smpl++ )
		{
			// compute two points on the membrane at random.
			
			int f1 = rand() % nt;
			double u1 = rand() / (double)RAND_MAX, v1 = rand() / (double)RAND_MAX;
			
			int f2 = rand() % nt;
			double u2 = rand() / (double)RAND_MAX, v2 = rand() / (double)RAND_MAX;

			while( u1 + v1 >= 1.0 )
			{
				u1 = rand() / (double)RAND_MAX;
				v1 = rand() / (double)RAND_MAX;
			}
			
			while( u2 + v2 >= 1.0 )
			{
				u2 = rand() / (double)RAND_MAX;
				v2 = rand() / (double)RAND_MAX;
			}

//			randomPointOnSurface(&f1,&u1,&v1);
//			randomPointOnSurface(&f2,&u2,&v2);

			// drawing points is somewhat expensive so let's do a bunch of z at this value.
				
			double r1[3], nrm1[3];
			double r2[3], nrm2[3];

			evaluateRNRM( f1, u1, v1, r1, nrm1, rmesh); 
			evaluateRNRM( f2, u2, v2, r2, nrm2, rmesh); 

			double w1 = g(f1,u1,v1,rmesh);
			double w2 = g(f2,u2,v2,rmesh);
			
			for( int nz = 0; nz < 20; nz++ )
			{
				double z1 = MINZ + (MAXZ-MINZ) * rand() / (double)RAND_MAX;
				double z2 = MINZ + (MAXZ-MINZ) * rand() / (double)RAND_MAX;
	
				double p1[3] = { r1[0] + z1 * nrm1[0],
						 r1[1] + z1 * nrm1[1],
						 r1[2] + z1 * nrm1[2] };
				
				double p2[3] = { r2[0] + z2 * nrm2[0],
						 r2[1] + z2 * nrm2[1],
						 r2[2] + z2 * nrm2[2] };

				double dr[3] = { p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2]};

				double rv = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

				int rbin = nbins * rv / maxr;

				double b1 = evaluateSpline( z1, SANS_SPLINE );
				double b2 = evaluateSpline( z2, SANS_SPLINE );
				
				double val = b1 * b2 * w1 * w2;

				if( rbin < nbins )
#ifdef PARALLEL					
					local_sampler[rbin] += val;
#else
					B_hist[rbin] += val;
#endif
			}
		}
	}	
	else if( sample_type == SANS_SAMPLE_LAT )
	{
		// this is sampling of the signal using the local scattering per unit area [the projection of beta(z) onto the surface]
		
	}

#ifdef PARALLEL
	ParallelSum( local_sampler, nbins );
	for( int b = 0; b < nbins; b++ )
		B_hist[b] += local_sampler[b];
	free(local_sampler);
#endif

}



void loadBetaZ( const char *fileName, int load_to_spline )
{
	double solvent_rho = 0;
	double nsolvent_rho = 0;
	double average_over_dist_from_box = 5.0;
	FILE *bzFile = fopen(fileName, "r");
   	char *buffer = (char *)malloc( sizeof(char) * 100000 );
 
	if( !bzFile )
	{   
	        printf("Requested beta.z file \"%s\" not found.\n", fileName );
	        exit(1);
	}   
	
	getLine( bzFile, buffer );

	int N_Z_BINS;
	double Lx, Ly, Lz;

	
	int nr = sscanf(buffer, "%d %lf %lf %lf\n", &N_Z_BINS, &Lx, &Ly, &Lz );

	MINZ = -Lz/2;
	MAXZ = Lz/2;
	
	getLine( bzFile, buffer );
	getLine( bzFile, buffer );

	int fp = ftell(bzFile);
	
	for( int z = 0; z < N_Z_BINS; z ++ )
	{
		double val = 0;
		double zbin;
        	getLine( bzFile, buffer );
        	double zv, rv;
        	sscanf( buffer, "%lf %lf", &zv, &rv );

		if( zv < -Lz/2 + average_over_dist_from_box || 
		    zv > Lz/2 - average_over_dist_from_box )
		{
			solvent_rho += rv;
			nsolvent_rho += 1;
		}
	}

	fseek( bzFile, fp, SEEK_SET );	
        
	solvent_rho /= nsolvent_rho;

	setupSpline( -Lz/2, Lz/2, N_Z_BINS, load_to_spline, 1 );

	for( int z = 0; z < N_Z_BINS; z ++ )
	{
		double val = 0;
		double zbin;
	
		double t_vals[2];
		
		getLine( bzFile, buffer );
		int nr = readNDoubles( buffer, t_vals, 2 );
			
		if( nr == 2 )
			AddPointToSpline( t_vals[0], t_vals[1] - solvent_rho, 0 );
		else
		{
			printf("Failed to read two numbers from line '%s' of file '%s'.\n", buffer, fileName );
			exit(1);
		}
	}

	SolveSpline(load_to_spline);

	printf("beta(z) loaded into spline %d.\n", load_to_spline );

	fflush(stdout);
	
}

void writeSq( char *fileName, double *B_hist, double sans_max_r, int nsans_bins, double qmin, double qmax, int nq )
{
#ifdef PARALLEL
	if( par_info.my_id != BASE_TASK ) return;
#endif
	FILE *theFile = fopen(fileName, "w");
	
	if( !theFile )
	{
		printf("Failed to open file '%s' for writing.\n", fileName );
		exit(1);
	}

	double dr = sans_max_r/nsans_bins;
	for( int ir = 0; ir < nsans_bins; ir++ )
		fprintf(theFile, "%le %le\n", (ir+0.5)*dr, B_hist[ir] );

	double eps = 1e-3;
	for( int iq = -1; iq < nq; iq++ )
	{
		double q = (qmin+qmax)/2; 

		if( nq > 1 )
			q = qmin + iq * (qmax-qmin)/(nq-1);

		double val = 0;

		if( iq == -1 )
		{
			q = 0;
			for( int ir = 0; ir < nsans_bins; ir++ )
			{
				double r = (ir+0.5)*dr;
	
				val += dr * B_hist[ir];
			}		
		}
		else
		{
			for( int ir = 0; ir < nsans_bins; ir++ )
			{
				double r = (ir+0.5)*dr;
	
				val += dr * B_hist[ir] * sin(q*r)/(q*r);
			}		
		}
	
		fprintf(theFile, "%le %le\n", q, val );
	}

	fflush(theFile);
	fclose(theFile);
}

