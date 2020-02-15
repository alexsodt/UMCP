#define __gatherc__
/*
 * Routines for gathering nearest-point curvature from a surface.
 * Initial purpose is to gather lipid's sampled curvature given a fit to the surface.
 * */
#include "simulation.h"
#include "interp.h"
#include "gather.h"
#include "M_matrix.h"
#include "meshCollisionLib.h"
#include "io_mol_read.h"
#include "input.h"
#include "mutil.h"

double Simulation::nearCurvature( double *rpt, double *c_out, double *k_out, double *dp_out, double *dz_out, int *leaflet_out )
{
	double **M;
	int mlow,mhigh;

	getM( &M, &mlow, &mhigh );

	double min_dist = 1e10;

	double min_c = 0;
	double min_k = 0;
	double min_dp = 0;
	double min_dz = 0;
	int min_leaflet = 0; // 0 is lower.
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		int f;
		double u, v;
		double distance;

		sRec->theSurface->nearPointOnBoxedSurface( rpt, &f, &u, &v, M, mlow, mhigh, &distance );					

		if( distance < min_dist )
		{
			min_c = sRec->theSurface->c( f, u, v, sRec->r, &min_k ); 
			min_dist = distance;

			// evaluate upper/lower

			double fpt[3],nrm[3];
			sRec->theSurface->evaluateRNRM( f, u, v, fpt, nrm, sRec->r );		

			if( sRec->gather_flip )
			{
				min_c *= -1;
				nrm[0] *= -1;
				nrm[1] *= -1;
				nrm[2] *= -1;
			}

			double dr[3] = { rpt[0] - fpt[0], rpt[1] - fpt[1], rpt[2] - fpt[2] };
			wrapPBC( dr, alpha );	

			double l = normalize(dr);
			min_dp = dr[0]*nrm[0]+dr[1]*nrm[1]+dr[2]*nrm[2];
			min_dz = l;
	
			if( min_dp < 0 )
			{
				min_leaflet = 0;
				min_c *= -1;
			}
			else
				min_leaflet = 1;	 
		}
	}	

	if( c_out ) *c_out = min_c;
	if( k_out ) *k_out = min_k;
	if( dp_out ) *dp_out = min_dp;
	if( dz_out ) *dz_out = min_dz;
	if( leaflet_out ) *leaflet_out = min_leaflet;

	return min_c;
}

void Simulation::gather( parameterBlock *block )
{
	for( surface_record *sRec = allSurfaces; sRec; sRec = sRec->next )
	{
		double av_c = 0;
		double av_k = 0;
		double sum_g = 0;	

		for( int f = 0; f < sRec->theSurface->nt; f++ )
		{
			double k;
			double g = sRec->theSurface->g(f, 1.0/3.0, 1.0/3.0, sRec->r ); 
			av_c += g * sRec->theSurface->c(f, 1.0/3.0, 1.0/3.0, sRec->r, &k ); 
			av_k += k * g;
			sum_g += g;
		}

		av_c /= sum_g;
		av_k /= sum_g;
		printf("Surface %d has <c> %lf <k> %lf ", sRec->id, av_c, av_k );
		if( av_c < 0 )
			printf("(flipping).\n");
		else
			printf("(not flipping).\n");

		if( av_c < 0 )
			sRec->gather_flip = 1;
		else
			sRec->gather_flip = 0;
	}	


        char *gathered = (char *)malloc( sizeof(char) * (1+strlen(block->jobName) + strlen("_gathered.xyz")) );
	sprintf(gathered, "%s_gathered.xyz", block->jobName );
	FILE *gatheredXYZ = fopen(gathered,"w");
	if( gatheredXYZ )		
	{
		writeLimitingSurface(gatheredXYZ);
		fclose(gatheredXYZ);
	}
	if( !block->structureName )
	{
		printf("Gathering requires a structure (e.g., PSF) file.\n");
		exit(1);
	}

	io_readStructureFile( block->structureName );	
	io_initialize_read( block->dcdName );

	char *gatherName = (char *)malloc( sizeof(char) * ( 1 + strlen(block->jobName) + strlen("_gather.txt") ) );

	sprintf(gatherName, "%s_gather.txt", block->jobName );

	FILE *gatherFile = fopen(gatherName, "w");

	if( !gatherFile )
	{
		printf("Couldn't open gather file '%s' for writing.\n", gatherName );
		exit(1);
	}

	printf("Starting gather.\n");

	int nframes = io_nframes();
	for( int f = 0; f < nframes; f++ )
	{
		io_readFrame( block->dcdName );
		printf("Gathering from frame %d.\n", f );
		fflush(stdout);
		int ns = io_nNSAtoms();
		io_align();
 
		struct atom_rec *at = io_getCurrentFrame();

		for( int ax = 0; ax < ns; ax ++ )
		{
			int a = io_getNSAtom(ax);
			double r[3] = { at[a].x, at[a].y, at[a].z };
			double cout,kout,dp_out,dz_out;
			int leaflet_out;
			double c = nearCurvature( r,&cout,&kout,&dp_out,&dz_out,&leaflet_out );
			
			fprintf(gatherFile, "%s %s %d %le %le %lf %lf %s\n", at[a].segid, at[a].resname, at[a].res, cout, kout, dp_out, dz_out, (leaflet_out == 0 ? "lower" : "upper" ) );	
			fflush(gatherFile);
		}	
	}	
		
	printf("Finished with gather.\n" );
	fflush(stdout);
} 

