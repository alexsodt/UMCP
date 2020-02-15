#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "clusterv.h"
#include "mesh_utility.h"

void getAltCylinderMesh( double **rvals, int **edges, int **nedges, int *nv, int grain, double R, double LC)
{ 
	double APL = 65.0;
	//double R = (nlipids*APL)/(2*M_PI*LC);

	double Circumf = 2 * M_PI * R;
	double Height  = LC;

	char buffer[4096];
	double *points;

	double hex_spacing = Circumf / grain;
	


//	LC / ny_int;
//	double hex_spacing = sqrt(Circumf * LC / (nlipids/2) / sqrt(3));
	double pt_spacing = 0.5 * hex_spacing / (sqrt(3.0)/2);
	double special_L = cos(30.0*M_PI/180.0) * 2 * hex_spacing;

	int nx_int = grain;
//	int ny_int = lround(LC/special_L);
	int ny_int = lround(LC/pt_spacing);


	double LA = nx_int * hex_spacing;
	double LB = ny_int * pt_spacing;

	double Ly_special = hex_spacing;

	int n_special_points = nx_int * ny_int;
	points = (double *)malloc( sizeof(double) * n_special_points * 3 );
	
	double PeriodicLengthA = LA;
	double PeriodicLengthB = LB;		

	n_special_points = 0;

	R = LA / (2*M_PI);

	for( int iy = 0; iy < ny_int; iy++ )
	{ // y is the axial direction.
		for( int ix = 0; ix < nx_int; ix++ )
		{ // x is the cirumferential direction.
			double x_val0 = ix * hex_spacing                 ;//                   - LA/2;
			double x_val1 = (ix+0.5) * hex_spacing		;//		    - LA/2;
//			double x_val1 = ix * special_L +  hex_spacing * cos(30.0*M_PI/180.0)- LA/2; 

			double theta0 = x_val0 / (R);
			double theta1 = x_val1 / (R);

			double r[3];
			if( iy % 2 == 1 )
			{
				r[0] = R * sin(theta0);
				r[1] = R * cos(theta0);
				r[2] = iy * pt_spacing;
			}
			else
			{
				r[0] = R * sin(theta1);
				r[1] = R * cos(theta1);
				r[2] = iy * pt_spacing;
			}
			
			//	printf("ok %d %d %lf %lf %lf\n", ix, iy, r[0], r[1], r[2] );
			for( int xp = 0; xp < 1; xp++ )
			{	
				int skip = 0;
				points[n_special_points*3+0] = r[xp*3+0];
				points[n_special_points*3+1] = r[xp*3+1];
				points[n_special_points*3+2] = r[xp*3+2];
				n_special_points++;
			}
		}
	}

	int npoints = n_special_points;
	double *r = points;
		
	int max_bonds = 20;

	int *bonds = (int *)malloc( sizeof(int) * npoints * max_bonds );
	int *nbonds = (int*)malloc( sizeof(int) * npoints );

	memset( nbonds, 0, sizeof(int) * npoints );
	double BoxL[2] = { PeriodicLengthA, PeriodicLengthB };
	double *rn_vec = (double *)malloc( sizeof(double) * npoints * 3 );
	for( int x = 0; x < npoints; x++ )
	{
		rn_vec[3*x+0] = (1e-2) * rand() / (double)RAND_MAX;
		rn_vec[3*x+1] = (1e-2) * rand() / (double)RAND_MAX;
		rn_vec[3*x+2] = (1e-2) * rand() / (double)RAND_MAX;

		r[3*x+0] += rn_vec[3*x+0];
		r[3*x+1] += rn_vec[3*x+1];
		r[3*x+2] += rn_vec[3*x+2];
	}
	
	{
		double BoxL[2] = { PeriodicLengthB, 2 * M_PI * R };
		double *r_mod = (double *)malloc( sizeof(double) * 3 * npoints );
		for( int x = 0; x < npoints; x++ )
		{
			r_mod[3*x+0] = r[3*x+2];
			double th = atan2( r[3*x+1], r[3*x+0] );
			r_mod[3*x+1] = R * th;
		}
		get2DVoronoiConnectivity(r_mod, npoints, "unique", bonds, nbonds, max_bonds, BoxL );
		free(r_mod);
	}

	for( int x = 0; x < npoints; x++ )
	{
		r[3*x+0] -= rn_vec[3*x+0];
		r[3*x+1] -= rn_vec[3*x+1];
		r[3*x+2] -= rn_vec[3*x+2];
	}

	freopen("/dev/tty", "w", stdout );
	
	int *tri = (int *)malloc( sizeof(int) * npoints * 3 * 6 );
	int ntri = 0; 

	int *new_bonds = (int *)malloc( sizeof(int) * npoints * max_bonds );
	int *n_new_bonds = (int*)malloc(sizeof(int) * npoints );
	memset( n_new_bonds, 0, sizeof(int) * npoints );


	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < nbonds[i]; jv++ )
		{		
			int consist = 0;
			int j = bonds[i*max_bonds+jv];
		
			for( int kv = 0; kv < nbonds[j]; kv++ )
			{
				if( bonds[j*max_bonds+kv] == i )
					consist = 1;
			}

			if(! consist )
			{
				printf("Inconsistent: %d bonded to %d but not %d to %d.\n", i, j, j, i );
			}
		}
	}
	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < nbonds[i]; jv++ )
		{
			int j = bonds[i*max_bonds+jv];
			if( j <= i ) continue;

			for( int kv = 0; kv < nbonds[j]; kv++ )
			{
				int k = bonds[j*max_bonds+kv];

				if( k <= j ) continue;
				
				int gotit = 0;
			
				for( int iv = 0; iv < nbonds[k]; iv++ )
				{
					int ic = bonds[k*max_bonds+iv];

					if( ic == i ) gotit = 1;
				}		
			
				if( gotit )
				{
					int gb = 0;
		
					int ipairs[6][2] = 
					{
						{i,j}, {i,k},
						{j,i}, {j,k},
						{k,i}, {k,j},
					};

					for( int pairs = 0; pairs< 6; pairs++ )
					{		
						int check_i = ipairs[pairs][0];
						int check_j = ipairs[pairs][1];	
						gb = 0;
						for( int p = 0; p < n_new_bonds[check_i]; p++ )
						{
							if( new_bonds[check_i*max_bonds+p] == check_j )
								gb = 1;
						}
	
						if( !gb )
						{
							new_bonds[check_i*max_bonds+n_new_bonds[check_i]] = check_j;
							n_new_bonds[check_i]++;
						}
					}
					tri[3*ntri+0] = i;
					tri[3*ntri+1] = j;
					tri[3*ntri+2] = k;
					ntri++;
				}
			} 
		}
	}

	// remove points that have no neighbors.

	int npts_use = 0;
	int pt_map[npoints];

	for( int i = 0; i < npoints; i++ )
	{
		pt_map[i] = -1;
		if( n_new_bonds[i] > 0 )
		{
			pt_map[i] = npts_use;
			npts_use++;
		}
/*		else
		{
			printf("Shouldn't happen now that we are reading points.\n");
			exit(1);
		}*/

	}
	

	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < n_new_bonds[i]; jv++ )
		{		
			int consist = 0;
			int j = new_bonds[i*max_bonds+jv];
		
			for( int kv = 0; kv < n_new_bonds[j]; kv++ )
			{
				if( new_bonds[j*max_bonds+kv] == i )
					consist = 1;
			}

			if(! consist )
			{
				printf("Inconsistent: %d bonded to %d but not %d to %d.\n", i, j, j, i );
			}
		}
	}

	int nit = 0;
	int op;

	int p = 0;

	int npt = 0;
	int net = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		npt++;

		net += n_new_bonds[p];
	}

	(*rvals)  = (double *)malloc( sizeof(double) * 3 * npt );
	(*edges)  = (int *)malloc( sizeof(int) * MAX_EDGES * npt );
	(*nedges) = (int *)malloc( sizeof(int) * npt );
	
	npt = 0;

	p = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		(*rvals)[3*npt+0] = r[3*p+0];
		(*rvals)[3*npt+1] = r[3*p+1];
		(*rvals)[3*npt+2] = r[3*p+2];

		for( int px = 0; px < n_new_bonds[p]; px++ )
			(*edges)[npt*MAX_EDGES+px] = pt_map[new_bonds[max_bonds*p+px]];

		(*nedges)[npt] = n_new_bonds[p];

		npt++;
	}
	
	*nv = npt;	

	free(points);
	free(bonds);
	free(new_bonds);
	free(n_new_bonds);
	free(tri);
/*
	char fileName[256];
	
	sprintf(fileName,"cylindrical.mesh");

	FILE *theLattice = fopen(fileName,"w");

	fprintf(theLattice, "3d R = %lf\n", R );
	fprintf(theLattice, "%lf 0.0 0.0\n", 4 * R );
	fprintf(theLattice, "0.0 %lf 0.0\n", 4 * R );
	fprintf(theLattice, "0.0 0.0 %lf\n", PeriodicLengthB );

	int p = 0;
	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		fprintf(theLattice, "%d %lf %lf %lf %d", pt_map[p], r[3*p+0], r[3*p+1], r[3*p+2] , n_new_bonds[p] );

		for( int b = 0; b < n_new_bonds[p]; b++ )
			fprintf(theLattice, " %d", pt_map[new_bonds[max_bonds*p+b]] );

		fprintf(theLattice, "\n"); 
	}

	fprintf(theLattice, "ntri %d\n", ntri );
*/
}

void getCylinderMesh( double **rvals, int **edges, int **nedges, int *nv, int grain, double R, double LC)
{ 
	int ny_int = grain;

	double APL = 65.0;
	//double R = (nlipids*APL)/(2*M_PI*LC);

	double Circumf = 2 * M_PI * R;
	double Height  = LC;

	char buffer[4096];
	double *points;

	double hex_spacing = LC / ny_int;
//	double hex_spacing = sqrt(Circumf * LC / (nlipids/2) / sqrt(3));
	double pt_spacing = 0.5 * hex_spacing / (sqrt(3.0)/2);
	double special_L = cos(30.0*M_PI/180.0) * 2 * hex_spacing;

	int nx_int = lround(Circumf / special_L/2)*2;
//	int ny_int = lround(LC/hex_spacing);



	double LA = nx_int * special_L;
	double LB = ny_int * hex_spacing;

	double Ly_special = hex_spacing;

	int n_special_points = nx_int * ny_int *2;
	points = (double *)malloc( sizeof(double) * n_special_points * 3 );
	
	double PeriodicLengthA = LA;
	double PeriodicLengthB = LB;		

	n_special_points = 0;

	R = LA / (2*M_PI);

	for( int ix = 0; ix < nx_int; ix++ )
	{ // x is the cirumferential direction.
		for( int iy = 0; iy < ny_int; iy++ )
		{ // y is the axial direction.
			double x_val0 = ix * special_L                                     - LA/2;
			double x_val1 = ix * special_L +  hex_spacing * cos(30.0*M_PI/180.0)- LA/2; 

			double theta0 = x_val0 / (R);
			double theta1 = x_val1 / (R);
			
			double r[6] = { R * sin(theta0), R*cos(theta0), iy * Ly_special              ,
		                        R * sin(theta1), R*cos(theta1), iy * Ly_special - hex_spacing * sin(30.0*M_PI/180.0)  };	

			for( int xp = 0; xp < 2; xp++ )
			{	
				int skip = 0;

				points[n_special_points*3+0] = r[xp*3+0];
				points[n_special_points*3+1] = r[xp*3+1];
				points[n_special_points*3+2] = r[xp*3+2];
				n_special_points++;
			}
		}
	}

	int npoints = n_special_points;
	double *r = points;
		
	int max_bonds = 20;

	int *bonds = (int *)malloc( sizeof(int) * npoints * max_bonds );
	int *nbonds = (int*)malloc( sizeof(int) * npoints );

	memset( nbonds, 0, sizeof(int) * npoints );
	double BoxL[2] = { PeriodicLengthA, PeriodicLengthB };
	double *rn_vec = (double *)malloc( sizeof(double) * npoints * 3 );
	for( int x = 0; x < npoints; x++ )
	{
		rn_vec[3*x+0] = (1e-2) * rand() / (double)RAND_MAX;
		rn_vec[3*x+1] = (1e-2) * rand() / (double)RAND_MAX;
		rn_vec[3*x+2] = (1e-2) * rand() / (double)RAND_MAX;

		r[3*x+0] += rn_vec[3*x+0];
		r[3*x+1] += rn_vec[3*x+1];
		r[3*x+2] += rn_vec[3*x+2];
	}
	
	{
		double BoxL[2] = { PeriodicLengthB, 2 * M_PI * R };
		double *r_mod = (double *)malloc( sizeof(double) * 3 * npoints );
		for( int x = 0; x < npoints; x++ )
		{
			r_mod[3*x+0] = r[3*x+2];
			double th = atan2( r[3*x+1], r[3*x+0] );
			r_mod[3*x+1] = R * th;
		}
		get2DVoronoiConnectivity(r_mod, npoints, "unique", bonds, nbonds, max_bonds, BoxL );
		free(r_mod);
	}

	for( int x = 0; x < npoints; x++ )
	{
		r[3*x+0] -= rn_vec[3*x+0];
		r[3*x+1] -= rn_vec[3*x+1];
		r[3*x+2] -= rn_vec[3*x+2];
	}

	freopen("/dev/tty", "w", stdout );
	
	int *tri = (int *)malloc( sizeof(int) * npoints * 3 * 6 );
	int ntri = 0; 

	int *new_bonds = (int *)malloc( sizeof(int) * npoints * max_bonds );
	int *n_new_bonds = (int*)malloc(sizeof(int) * npoints );
	memset( n_new_bonds, 0, sizeof(int) * npoints );


	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < nbonds[i]; jv++ )
		{		
			int consist = 0;
			int j = bonds[i*max_bonds+jv];
		
			for( int kv = 0; kv < nbonds[j]; kv++ )
			{
				if( bonds[j*max_bonds+kv] == i )
					consist = 1;
			}

			if(! consist )
			{
				printf("Inconsistent: %d bonded to %d but not %d to %d.\n", i, j, j, i );
			}
		}
	}
	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < nbonds[i]; jv++ )
		{
			int j = bonds[i*max_bonds+jv];
			if( j <= i ) continue;

			for( int kv = 0; kv < nbonds[j]; kv++ )
			{
				int k = bonds[j*max_bonds+kv];

				if( k <= j ) continue;
				
				int gotit = 0;
			
				for( int iv = 0; iv < nbonds[k]; iv++ )
				{
					int ic = bonds[k*max_bonds+iv];

					if( ic == i ) gotit = 1;
				}		
			
				if( gotit )
				{
					int gb = 0;
		
					int ipairs[6][2] = 
					{
						{i,j}, {i,k},
						{j,i}, {j,k},
						{k,i}, {k,j},
					};

					for( int pairs = 0; pairs< 6; pairs++ )
					{		
						int check_i = ipairs[pairs][0];
						int check_j = ipairs[pairs][1];	
						gb = 0;
						for( int p = 0; p < n_new_bonds[check_i]; p++ )
						{
							if( new_bonds[check_i*max_bonds+p] == check_j )
								gb = 1;
						}
	
						if( !gb )
						{
							new_bonds[check_i*max_bonds+n_new_bonds[check_i]] = check_j;
							n_new_bonds[check_i]++;
						}
					}
					tri[3*ntri+0] = i;
					tri[3*ntri+1] = j;
					tri[3*ntri+2] = k;
					ntri++;
				}
			} 
		}
	}

	// remove points that have no neighbors.

	int npts_use = 0;
	int pt_map[npoints];

	for( int i = 0; i < npoints; i++ )
	{
		pt_map[i] = -1;
		if( n_new_bonds[i] > 0 )
		{
			pt_map[i] = npts_use;
			npts_use++;
		}
/*		else
		{
			printf("Shouldn't happen now that we are reading points.\n");
			exit(1);
		}*/

	}
	

	for( int i = 0; i < npoints; i++ )
	{
		for( int jv = 0; jv < n_new_bonds[i]; jv++ )
		{		
			int consist = 0;
			int j = new_bonds[i*max_bonds+jv];
		
			for( int kv = 0; kv < n_new_bonds[j]; kv++ )
			{
				if( new_bonds[j*max_bonds+kv] == i )
					consist = 1;
			}

			if(! consist )
			{
				printf("Inconsistent: %d bonded to %d but not %d to %d.\n", i, j, j, i );
			}
		}
	}

	int nit = 0;
	int op;

	int p = 0;

	int npt = 0;
	int net = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		npt++;

		net += n_new_bonds[p];
	}

	(*rvals)  = (double *)malloc( sizeof(double) * 3 * npt );
	(*edges)  = (int *)malloc( sizeof(int) * MAX_EDGES * npt );
	(*nedges) = (int *)malloc( sizeof(int) * npt );
	
	npt = 0;

	p = 0;

	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		(*rvals)[3*npt+0] = r[3*p+0];
		(*rvals)[3*npt+1] = r[3*p+1];
		(*rvals)[3*npt+2] = r[3*p+2];

		for( int px = 0; px < n_new_bonds[p]; px++ )
			(*edges)[npt*MAX_EDGES+px] = pt_map[new_bonds[max_bonds*p+px]];

		(*nedges)[npt] = n_new_bonds[p];

		npt++;
	}
	
	*nv = npt;	
	printf("npt: %d\n", npt );

	free(points);
	free(bonds);
	free(new_bonds);
	free(n_new_bonds);
	free(tri);
/*
	char fileName[256];
	
	sprintf(fileName,"cylindrical.mesh");

	FILE *theLattice = fopen(fileName,"w");

	fprintf(theLattice, "3d R = %lf\n", R );
	fprintf(theLattice, "%lf 0.0 0.0\n", 4 * R );
	fprintf(theLattice, "0.0 %lf 0.0\n", 4 * R );
	fprintf(theLattice, "0.0 0.0 %lf\n", PeriodicLengthB );

	int p = 0;
	for( ; p < npoints; p++ )
	{
		if( pt_map[p] < 0 )
			continue;

		fprintf(theLattice, "%d %lf %lf %lf %d", pt_map[p], r[3*p+0], r[3*p+1], r[3*p+2] , n_new_bonds[p] );

		for( int b = 0; b < n_new_bonds[p]; b++ )
			fprintf(theLattice, " %d", pt_map[new_bonds[max_bonds*p+b]] );

		fprintf(theLattice, "\n"); 
	}

	fprintf(theLattice, "ntri %d\n", ntri );
*/
}
