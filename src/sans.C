
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

	int nSurfaces = 0;
			
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
		nSurfaces++;

	// we will draw points from the surfaces such that we target their total areas.
	double tarea[nSurfaces*2];
	memset( tarea, 0, sizeof(double) * nSurfaces * 2);
	// compute each total area.

	int n_draw = 1000;

	double **draw_pts = (double **)malloc( sizeof(double) * nSurfaces );

	for( int x = 0; x < nSurfaces; x++ )
		draw_pts[x] = (double *)malloc( sizeof(double) * 4 * n_draw * draw_z );
 
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		double area0;
		sRec->theSurface->area( sRec->r, -1, tarea+sRec->id*2, &area0 ); 
				
		surface *theSurface = sRec->theSurface;
		int nt = theSurface->nt;
		double *rmesh = sRec->r;
		double total_area = 0;
		int cur_pt = 0;
		for( int d = 0; d < n_draw; d++ )
		{
			int f1 = rand() % nt;
			double u1 = rand() / (double)RAND_MAX, v1 = rand() / (double)RAND_MAX;

			while( u1 + v1 >= 1.0 )
			{
				u1 = rand() / (double)RAND_MAX;
				v1 = rand() / (double)RAND_MAX;
			}
				
			double r1[3], nrm1[3];

			theSurface->evaluateRNRM( f1, u1, v1, r1, nrm1, rmesh); 

			double w1 = theSurface->g(f1,u1,v1,rmesh);
			double c1  = 0;

			double k;
			if( shape_correction )
				c1 = theSurface->c(f1,u1,v1,rmesh,&k);

			total_area += w1;
			for( int nz = 0; nz < draw_z; nz++)
			{
				double z1 = MINZ + (MAXZ-MINZ) * rand() / (double)RAND_MAX;

				double alpha1 = 1.0;

				if( shape_correction )
				{
					if( z1 > 0 )
						alpha1 = exp( -(z1-PP)*c1);
					else
						alpha1 = exp( -(z1+PP)*c1);
				}

				double p1[3] = { r1[0] + z1 * nrm1[0],
						 r1[1] + z1 * nrm1[1],
						 r1[2] + z1 * nrm1[2] };
				
				double b1 = alpha1*evaluateSpline( z1*alpha1, SANS_SPLINE );

				draw_pts[sRec->id][cur_pt*4+0] = p1[0];	
				draw_pts[sRec->id][cur_pt*4+1] = p1[1];	
				draw_pts[sRec->id][cur_pt*4+2] = p1[2];	
				draw_pts[sRec->id][cur_pt*4+3] = b1 * w1 * dz;

				cur_pt++;
			}
		}
		tarea[2*sRec->id+1] = total_area;
	}	

#ifdef TOTAL_DRAW
	// do every single pair we've got. 
	// not necessarily the "more correct" thing to do, but guaranteed to be positive.
	for( int s1 = 0; s1 < nSurfaces; s1++ )
	for( int s2 = s1; s2 < nSurfaces; s2++ )
	{
		double *draw1 = draw_pts[s1];
		double *draw2 = draw_pts[s2];
		double fac = 1.0;

		if( s1 != s2 )
			fac = 2.0;

		for( int p1 = 0; p1 < n_draw * draw_z; p1++ )
		{
			for( int p2 = 0; p2 < n_draw * draw_z; p2++ )
			{
				double *rp1 = draw1 + p1 * 4;	
				double *rp2 = draw2 + p2 * 4;	
		
				double dr[3] = { rp1[0]-rp2[0],rp1[1]-rp2[1],rp1[2]-rp2[2]};
	
				double rv = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

				int rbin = nbins * rv / maxr;
					
				if( rbin < nbins )
				{
					double val = fac * rp1[3] * rp2[3] * (tarea[s1*2+0]/tarea[s1*2+1])*(tarea[s2*2+0]/tarea[s2*2+1]);
#ifdef PARALLEL					
					local_sampler[rbin] += val;
#else
					B_hist[rbin] += val;
#endif
				}
			}
		}
	}
#else
	int n_pairs = 10000000;

	for( int p = 0; p < n_pairs; p++ )
	{
		int s1 = rand() % nSurfaces;
		int s2 = rand() % nSurfaces;

		int p1 = rand() % n_draw;
		int z1 = rand() % draw_z;
		int p2 = rand() % n_draw;
		int z2 = rand() % draw_z;

		if( p1 == p2 ) continue;

		double *rp1 = draw_pts[s1] + (p1*draw_z+z1) * 4;	
		double *rp2 = draw_pts[s2] + (p2*draw_z+z2) * 4;	

		double dr[3] = { rp1[0]-rp2[0],rp1[1]-rp2[1],rp1[2]-rp2[2]};

		double rv = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

		int rbin = nbins * rv / maxr;
			
		if( rbin < nbins )
		{
			double val = rp1[3] * rp2[3] * (tarea[s1*2+0]/tarea[s1*2+1])*(tarea[s2*2+0]/tarea[s2*2+1]);
#ifdef PARALLEL					
			local_sampler[rbin] += val;
#else
			B_hist[rbin] += val;
#endif
		}
	
	}
#endif

	for( int x = 0; x < nSurfaces; x++ )
		free(draw_pts[x]);

	free(draw_pts);

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

