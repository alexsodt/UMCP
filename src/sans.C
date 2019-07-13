
//#define DEBUG_PTS

#include "interp.h"
#include "sans.h"
#include "spline.h"
#include "util.h"
#include "input.h"
#include <math.h>
#include "parallel.h"
#include <string.h>
#include "simulation.h"
void loadBetaZ( const char *fileName, int load_to_spline );

static int N_Z_BINS = 100;
static double MINZ = -75;
static double MAXZ =  75;
double Ap = 1.0; // square angstroms

void initSANS( parameterBlock *block, double **qvals, int *nq )
{
	if( block->s_q )
	{
		loadBetaZ( block->betazFile, SANS_SPLINE );
		*qvals = NULL;

		if( block->qvals )
		{
			FILE *qvalsFile = fopen(block->qvals, "r" );

			if( !qvalsFile )
			{
				printf("Failed to open q value file \"%s\".\n", block->qvals );
				exit(1);
			}

			char *buffer = (char *)malloc( sizeof(char)*100000 );
			int qvalSP = 10;
			int nq_local = 0;
			*qvals = (double *)malloc( sizeof(double) * qvalSP );
			while( !feof(qvalsFile) )
			{
				getLine(qvalsFile, buffer );
				if( feof(qvalsFile) ) break;
				double qv;
				int nr = sscanf(buffer, "%lf", & qv );
				if( nr == 1 )
				{
					if( qvalSP == nq_local )
					{
						qvalSP *= 2;
						*qvals = (double *)realloc( *qvals, sizeof(double) * qvalSP );
					}
					(*qvals)[nq_local] = qv;
					nq_local++;
				} 
			}
			*nq = nq_local;
			free(buffer);
		}
		else
			*nq = block->nq;
	}
}

// a measurement of Sq is added into Sq_inst (it is not zero'd).

void Simulation::sample_B_hist( double *B_hist, double *A2dz2_sampled,
			int sample_type, 
			int nsamples, double maxr, int nbins, int shape_correction )
{
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		surface *theSurface = sRec->theSurface;
		int nt = theSurface->nt;
		double *rmesh = sRec->r;

		static double av_b_sampled = 0;
		static double av_b2_sampled = 0;
		static double n_samples = 0;

		double special_q = 6.494494e-03;
		static double av_r[3] = {0,0,0};
		static double av_c[3] = {0,0,0};
		
		double local_av_r[3] = {0,0,0};
		double local_av_c[3] = {0,0,0};

#ifdef DEBUG_PTS
		FILE *dbgFile = fopen("debug.xyz","w");
		double max_area = 100000;
#endif
		int samples_per_proc = ceil( nsamples/ (double)par_info.nprocs );

#ifdef PARALLEL
		double *local_sampler = (double *)malloc( sizeof(double) * nbins );
		memset( local_sampler, 0, sizeof(double) * nbins );
#endif
		

		double PP = 12.0; // approximately.

		// for clarity in accounting for actual volume contribution we use a real value of dz/Ap:
		double dz = (MAXZ-MINZ)/N_Z_BINS;

		double local_A2dz2_sampled = 0;

		int draw_z = 20;

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

				theSurface->evaluateRNRM( f1, u1, v1, r1, nrm1, rmesh); 
				theSurface->evaluateRNRM( f2, u2, v2, r2, nrm2, rmesh); 

				double w1 = theSurface->g(f1,u1,v1,rmesh);
				double w2 = theSurface->g(f2,u2,v2,rmesh);
				double c1  = 0;
				double c2  = 0;

				if( shape_correction )
				{
					c1 = theSurface->c(f1,u1,v1,rmesh);
                        	c2 = theSurface->c(f2,u2,v2,rmesh);			
				}
				for( int nz = 0; nz < draw_z; nz++ )
				{
					double z1 = MINZ + (MAXZ-MINZ) * rand() / (double)RAND_MAX;
					double z2 = MINZ + (MAXZ-MINZ) * rand() / (double)RAND_MAX;

					double alpha1 = 1.0;
					double alpha2 = 1.0;

					if( shape_correction )
					{
						if( z1 > 0 )
							alpha1 = exp( -(z1-PP)*c1);
						else
							alpha1 = exp( -(z1+PP)*c1);

						if( z2 > 0 )
							alpha2 = exp( -(z2-PP)*c2);
						else
							alpha2 = exp( -(z2+PP)*c2);
					}

					double p1[3] = { r1[0] + z1 * nrm1[0],
							 r1[1] + z1 * nrm1[1],
							 r1[2] + z1 * nrm1[2] };
					
					double p2[3] = { r2[0] + z2 * nrm2[0],
							 r2[1] + z2 * nrm2[1],
							 r2[2] + z2 * nrm2[2] };

					double dr[3] = { p1[0]-p2[0],p1[1]-p2[1],p1[2]-p2[2]};

					double rv = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

					int rbin = nbins * rv / maxr;

					double b1 = alpha1*evaluateSpline( z1*alpha1, SANS_SPLINE );
					double b2 = alpha2*evaluateSpline( z2*alpha2, SANS_SPLINE );

					av_b_sampled += (b1+b2)/2;
					av_b2_sampled += b1*b2;
					n_samples += 1;

					// number of volume ``units'' we put into histogram;
					double vol_squared_weight = (w1) * (w2) * dz * dz / (Ap*dz) / (Ap*dz);				
					double val = b1 * b2 * vol_squared_weight;
					(local_A2dz2_sampled) += w1 * w2;
		
					av_r[0] += b1 * cos( p1[0] * special_q) * w1;
					av_r[1] += b1 * cos( p1[1] * special_q) * w1;
					av_r[2] += b1 * cos( p1[2] * special_q) * w1;
					
					av_c[0] += b1 * sin( p1[0] * special_q) * w1;
					av_c[1] += b1 * sin( p1[1] * special_q) * w1;
					av_c[2] += b1 * sin( p1[2] * special_q) * w1;
					
					av_r[0] += b2 * cos( p2[0] * special_q) * w2;
					av_r[1] += b2 * cos( p2[1] * special_q) * w2;
					av_r[2] += b2 * cos( p2[2] * special_q) * w2;
					
					av_c[0] += b2 * sin( p2[0] * special_q) * w2;
					av_c[1] += b2 * sin( p2[1] * special_q) * w2;
					av_c[2] += b2 * sin( p2[2] * special_q) * w2;


#ifdef DEBUG_PTS	
					if( nz == 0 )
					{
						double rn = rand()/(double)RAND_MAX;

						if( rn < w1 / max_area )
							fprintf(dbgFile, "C %lf %lf %lf %d %lf %lf %lf\n", p1[0], p1[1], p1[2], f1, u1, v1, w1 );
						if( rn < w2 / max_area )
							fprintf(dbgFile, "C %lf %lf %lf %d %lf %lf %lf\n", p2[0], p2[1], p2[2], f2, u2, v2, w2 );
					}
#endif
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

#ifdef DEBUG_PTS
		fclose(dbgFile);
		exit(1);
#endif

#ifdef PARALLEL
		ParallelSum( local_av_r, 3 );
		ParallelSum( local_av_c, 3 );

		av_r[0] += local_av_r[0];
		av_r[1] += local_av_r[1];
		av_r[2] += local_av_r[2];
		
		av_c[0] += local_av_c[0];
		av_c[1] += local_av_c[1];
		av_c[2] += local_av_c[2];

		printf("Specialq: %le %le %le\n", 
				av_r[0]*av_r[0]+av_c[0]*av_c[0],
				av_r[1]*av_r[1]+av_c[1]*av_c[1],
				av_r[2]*av_r[2]+av_c[2]*av_c[2] );

		ParallelSum( &local_A2dz2_sampled, 1 );
		ParallelSum( local_sampler, nbins );
		for( int b = 0; b < nbins; b++ )
			B_hist[b] += local_sampler[b];
		free(local_sampler);
#endif

		*A2dz2_sampled += local_A2dz2_sampled;
 	}
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

void writeSq( char *fileName, double *B_hist, double A2dz2_sampled, double sans_max_r, int nsans_bins, double qmin, double qmax, int nq, double *qvals )
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
//	for( int ir = 0; ir < nsans_bins; ir++ )
//		fprintf(theFile, "%le %le\n", (ir+0.5)*dr, B_hist[ir] );


	// update histogram B0.

	double B0 = 0;

	

	double eps = 1e-3;
	for( int iq = -1; iq < nq; iq++ )
	{
		double q;

		if( qvals )
		{
			if( iq == -1 ) 
				q = 0;
			else
				q = qvals[iq];
		}
		else
		{
			q = (qmin+qmax)/2; 
	
			if( nq > 1 )
				q = qmin + iq * (qmax-qmin)/(nq-1);
		}

		double val = B0;

		if( iq == -1 )
		{
			q = 0;
			for( int ir = 0; ir < nsans_bins; ir++ )
			{

				double r = (ir+0.5)*dr;
	
				val += B_hist[ir];
			}		
		}
		else
		{
			for( int ir = 0; ir < nsans_bins; ir++ )
			{
				double r = (ir+0.5)*dr;
	
				val += B_hist[ir] * sin(q*r)/(q*r);
			}		
		}
	
		fprintf(theFile, "%le %le\n", q, val );
	}

	fflush(theFile);
	fclose(theFile);
}

