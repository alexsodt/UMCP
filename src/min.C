//#include <OpenCL/opencl.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "l-bfgs.h"
#include "m_triangles.h"
#include "parallel.h"
//#define SKIP_OPT

#define FIX_VOLUME

surface *minSurface = NULL;
double thick = 13.8;
extern int do_print_area_ps;
extern double av_vol_err;
extern double KV_scale;
extern double KA ;
void fdiff_check( double *coords );
int print_energy = 0;
extern double tilt_scale;
surface *upperSurface;
surface *lowerSurface;
volumes *theVolumes; 
 void amoeba(double **p, double y[], int ndim, double ftol,
               double (*funk)(double []), int *nfunk);
void amoebaOpt( double *r, int nv, double *rc, int nrc );
double surface_grad( double *r, double *g );
double surface_e( double *r );
	void curvatureGrad( double *g, surface *surface1, surface *surface2, double *r, double thick );
	double curvatureEnergy( surface *surface1, surface *surface2, double *r, double thick);
	void tiltGrad( double *g, surface *surface1, surface *surface2, double *r );
	void tiltGradExt( double *g, surface *surface1, surface *surface2, double *r );
	double tiltEnergy( surface *surface1, surface *surface2, double *r);
	double tiltEnergyExt( surface *surface1, surface *surface2, double *r, double *cons_energy);
double f(double *r);
double fdf( double *r, double *g);
void writeDualSurface( const char *fileNameXYZ, const char *fileNamePSF, surface *upperSurface, surface *lowerSurface, double *coords, double alpha );

int fix_alpha = 0;

int main( int argc, char **argv )
{
#ifdef PARALLEL
	int ierr;
	ierr = MPI_Init(&argc,&argv);
	if( ierr ) { exit(1); }
	quietParallel();
#endif
#ifdef FIXED_SEED
	srand(1);
#else
          struct timeval tp;
 
          gettimeofday( &tp, NULL );
 
          srand((long)tp.tv_usec);
#endif

	char buffer[4096];

	if( argc < 2 )
	{
		printf("Syntax: min lattice\n");
		return -1;
	}
	surface *theSurface =(surface *)malloc( sizeof(surface) );
	theSurface->loadLattice( argv[1], 0. );

	theSurface->generatePlan();
	
	double *r = (double *)malloc( sizeof(double) * (3 * theSurface->nv+3) );
	double *g = (double *)malloc( sizeof(double) * (3 * theSurface->nv+3) );
	theSurface->get(r);
	r[theSurface->nv*3] = 1.0;
	r[theSurface->nv*3+1] = 1.0;
	r[theSurface->nv*3+2] = 1.0;
	theSurface->setg0(r);

	KA = 1.215;

//	theSurface->fdiff_check_grad(r);
	double vol = theSurface->volume(r);

	double area0;
	double cur_area;
	theSurface->area(r, -1, &cur_area, &area0 );
	printf("area: %le area0: %le\n", cur_area, area0 );

	printf("vol: %lf\n", vol );
	
	if( argc > 2 )
	{
		FILE *loadFile = fopen(argv[2], "r");
	
		for( int x = 0; x < theSurface->nv; x++ )
		{
			getLine(loadFile, buffer );

			sscanf( buffer, "%lf %lf %lf\n",
				r+3*x+0, r+3*x+1, r+3*x+2 );
		}
		fclose(loadFile);
	}
	double rand_mag = 0;
	for( int v = 0; v < theSurface->nv; v++ )
	{
		r[3*v+0] += rand_mag * 2 * ( rand() / (double)RAND_MAX - 0.5);
		r[3*v+1] += rand_mag * 2 * ( rand() / (double)RAND_MAX - 0.5);
		r[3*v+2] += rand_mag * 2 * ( rand() / (double)RAND_MAX - 0.5);
	}
	
	minSurface = theSurface;
	FILE *tFile = fopen("min.xyz","w");
//	fdiff_check(r);

	 FILE *tpsf = fopen("min.psf","w");
	minSurface->writeLimitPSF(tpsf);
	fclose(tpsf);

	minSurface = theSurface;
	for( int o = 0; o < 5; o++ )
	{
		l_bfgs_setup( 20, 3*minSurface->nv+3, r, 0, f, fdf); 


		int n_iter = 0;

		int x = 0;
		do
		{
//			print_energy = 1;
//			f(coords);
//			print_energy = 0;

			n_iter++;
			if( n_iter > 50 ) 
				break;

			if( x % 10 == 0 )
			{	
				theSurface->put(r);
				theSurface->writeLimitStructure(tFile);
				fflush(tFile);
			}
			x++;

		} while ( l_bfgs_iteration( r ) );
		
//		fdiff_check(r);
		
	
//		fdiff_check(coords);
		
		l_bfgs_clear();

	}
//	fdiff_check(r);
	fclose(tFile);

	FILE *saveFile = fopen("file.save", "w");

	for( int x = 0; x < theSurface->nv; x++ )
	{
		fprintf(saveFile, "%.14le %.14le %.14le\n",
			r[3*x+0], r[3*x+1], r[3*x+2] );
	}
	fclose(saveFile);
	
	theSurface->saveSurface("min.mesh");

}

double vext = 0;
extern double VA, VC;
double f( double *r )
{
	vext = 0;
#ifdef EXTERNAL
	vext = minSurface->externalPotential(r);
#endif
	return minSurface->energy(r,NULL) + vext;
}

double fdf( double *r, double *g )
{
	double en = f(r);
	static int xx = 0;
//	if( xx % 100 == 0 )
		printf("fdf e: %.14le VE: %.14le AE: %.14le CE: %.14le\n", en, vext, VA, VC );
	xx++;
	memset( g, 0, sizeof(double) * (3 * minSurface->nv+3) );
	minSurface->grad(r,g);
#ifdef FIX_VOLUME
	g[3*minSurface->nv+0] = 0;
	g[3*minSurface->nv+1] = 0;
	g[3*minSurface->nv+2] = 0;
#endif
#ifdef EXTERNAL
	minSurface->externalPotentialGrad(r,g);
#endif
	return en;
}
	
void fdiff_check( double *r )
{
	double *g = (double *)malloc( sizeof(double) * (3 * minSurface->nv+3) );
	memset( g, 0, sizeof(double) * (3 * minSurface->nv+3) );
//	minSurface->grad(r,g);

	fdf( r, g );

	printf("fdiff check\n");
	double eps = 1e-5;
	for( int i = 0; i < minSurface->nv; i++ )
	{
		for( int c = 0; c < 3; c++ )
		{
			double ens[2];
			for( int di = 0; di < 2; di++ )
			{
				r[i*3+c] += -eps/2 + di * eps;
				ens[di] = f(r);
				r[i*3+c] -= -eps/2 + di * eps;
			}
			
			printf("i: %d c: %d fd: %.14le g: %.14le del: %.14le ",
				i, c, (ens[1]-ens[0])/eps, g[i*3+c], fabs((ens[1]-ens[0])/eps-g[i*3+c]) );
			if( fabs((ens[1]-ens[0])/eps-g[i*3+c]) < 1e-3 )
				printf(" OK \n");
			else
				printf(" BAD \n");
		}
	}
	free(g);
}
