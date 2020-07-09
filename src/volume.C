#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "interp.h"
#include "mutil.h"

//#define DISABLE_ALL
#define HARD_EPS (1e-10)
#define SOFT_EPS (1e-6)

int ntally_target = 10;
int nbuild_pts_tallied = 0;

//#define DEBUG_VOLUME
#ifdef DEBUG_VOLUME
FILE *edgePts = NULL;
int edge_space = 10;
int nbonds = 0;
int *edge_bonds=NULL;
int nat = 0;
double alt_vol = 0;
double target_area_A[6] = {0,0,0,0,0,0};
double target_area_B[6] = {0,0,0,0,0,0};
double vol_q[27];
static int debug_trigger = 0;
#endif

static FILE *debug_triangles = NULL;
static int ncalls = 0;


/*
 * uses the divergence-theorem to compute an approximate volume either inside or outside of the surface.
 *
 * that's straightforward except for the periodic boundary conditions.
 *
 * */

struct edge_data
{	
	// indexing of the 12 edges:
	// 0-3:
	
	// y == 0, z == 0, x
	// y == 0, z == Lz, x
	// y == Ly, z == 0, x
	// y == Ly, z == Lz, x
	
	// 4-7:
	
	// x == 0, z == 0, y
	// x == 0, z == Lz, y
	// x == Lx, z == 0, y
	// x == Lx, z == Lz, y

	// 8-11:
	
	// x == 0, y == 0, z
	// x == 0, y == Ly, z
	// x == Lx, y == 0, z
	// x == Lx, y == Ly, z

	// a number between zero and one. If it is -1, there is no point there. 
	double cur_edge_pointing_neg[12];
	double cur_edge_pointing_pos[12];	
};

double splitTrianglesOnPBC( double *verts, double *nrm, 
		double Lx, double Ly, double Lz, edge_data *pbc_edge_data, double *dvol, double *dnrm ); // , double *face_factors, double *spec_points, int *npoints);
void drawVMDArrow( FILE *theFile, double *from, double *to );
// fills cross-product derivative
void fillcp( double *der, double *r );
void fillcpder( double *der, double *dr, 
		int wrt,
		double sign
		);
void darea( double *dvol, double *r1, double *r2, double *r3, double fac );
void dcross( double *dvol, double *r1, double *r2, double *der1, double *der2, double factor, double *to_dot );

void surface::new_volume( double *vol_inside, double *vol_outside, double *rsurf, double *alpha, 
	// gradient of volume (inside) wrt mesh verts:
	double *dvol )
{
	if( !theVolumeFormulas )
		generateVolumePlan();
	
	if( dvol )
		memset( dvol, 0, sizeof(double) * 3 * nv );
	
	double Lx = PBC_vec[0][0];
	double Ly = PBC_vec[1][1];
	double Lz = PBC_vec[2][2];

	*vol_inside  = 0;
	*vol_outside = 0;

	nbuild_pts_tallied = 0;

#ifdef DEBUG_VOLUME
	nat = 0;
	nbonds = 0;
	target_area_A[0]=0;
	target_area_A[1]=0;
	target_area_A[2]=0;
	target_area_A[3]=0;
	target_area_A[4]=0;
	target_area_A[5]=0;
	target_area_B[0]=0;
	target_area_B[1]=0;
	target_area_B[2]=0;
	target_area_B[3]=0;
	target_area_B[4]=0;
	target_area_B[5]=0;
	debug_triangles = fopen("volume_triangles.vmd","w");
	edgePts = fopen("edge_pts.xyz","w");
	edge_bonds = (int *)malloc( sizeof(int) * 2 * edge_space );
	memset( vol_q, 0, sizeof(double) * 27 );
#if 0	
	double verts[9];

	verts[0] = 1+189.9;
	verts[1] = 0+189.9;
	verts[2] = 0+189.9;
	
	verts[3] = 0+189.9;
	verts[4] = 1+189.9;
	verts[5] = 0+189.9;
	
	verts[6] = 0+189.9;
	verts[7] = 0+189.9;
	verts[8] = 1+189.9;

	double nrm[3];
	double dr1[3] = { verts[3]-verts[0], verts[4]-verts[1],verts[5]-verts[2] };
	double dr2[3] = { verts[6]-verts[0], verts[7]-verts[1],verts[8]-verts[2] };
	cross( dr1, dr2, nrm);
	normalize(nrm);
	double face_factors[6];
	double spec_points[9];
	int nspec=0;

	splitTrianglesOnPBC( verts, nrm,
						Lx, Ly, Lz );
	exit(1);
#endif

#endif	

	double vol = 0;	

	double rav = 0;
	double nav = 0;
	
	int has_pbc[3] = {0,0,0};

	// each PBC face has an area computed by the 2D divergence theorem.
	double PBC_face_areas[6] = { 0,0,0,0,0,0 };

	double shift[3] = { 0,0,0};
	
	for( int t = 0; t < nv; t++ )
	{
//		shift[0] += -rsurf[3*t+0]; 
//		shift[1] += -rsurf[3*t+1]; 
//		shift[2] += -rsurf[3*t+2]; 
	}
	
	shift[0] /= nv;
	shift[1] /= nv;
	shift[2] /= nv;

//big effect
//	shift[0] += M_PI*2.4*1.20;
//	shift[1] += M_PI*1.5;
//	shift[2] += M_PI*4.5*1.05;


//big effect
//	shift[0] += M_PI*2.4*1.05;
//	shift[1] += M_PI*2.2;
//	shift[2] += M_PI*4.2*1.05;

	shift[0] += PBC_vec[0][0]/2;
	shift[1] += PBC_vec[1][1]/2;
	shift[2] += PBC_vec[2][2]/2;

	// jostle it a bit to get things off edges if the user did that by construction.
	shift[0] += (1e-5) * rand()/(double)RAND_MAX;
	shift[1] += (1e-5) * rand()/(double)RAND_MAX;
	shift[2] += (1e-5) * rand()/(double)RAND_MAX;

	edge_data pbc_edge_data;

	for( int edge = 0; edge < 12; edge++ )
	{
		// > 1 and it will be neglected.
		pbc_edge_data.cur_edge_pointing_neg[edge] = 2;
		pbc_edge_data.cur_edge_pointing_pos[edge] = 2;
	}
	for( int f = 0; f < nfaces; f++ )
	{
		int ni = theVolumeFormulas[f].ni;

		int nv_space = ni*(ni+1)/2;

		double coords[3*ni*ni];

		// based on the control points most of the time we don't even need to check PBCs.
		// the triangle will be completely contained in the convex hull of its control points.
		int activate_PBC[3] = {0,0,0};
		int *cp = theVolumeFormulas[f].cp;
		int np = theVolumeFormulas[f].ncoor;

		for( int p = 0; p < np; p++ )
		{
			if( rsurf[3*cp[p]+0] + shift[0] + theVolumeFormulas[f].r_pbc[3*p+0] < 0 ) activate_PBC[0] = 1;		
			if( rsurf[3*cp[p]+0] + shift[0] + theVolumeFormulas[f].r_pbc[3*p+0] > Lx ) activate_PBC[0] = 1;		
			if( rsurf[3*cp[p]+1] + shift[1] + theVolumeFormulas[f].r_pbc[3*p+1] < 0 ) activate_PBC[1] = 1;		
			if( rsurf[3*cp[p]+1] + shift[1] + theVolumeFormulas[f].r_pbc[3*p+1] > Ly ) activate_PBC[1] = 1;		
			if( rsurf[3*cp[p]+2] + shift[2] + theVolumeFormulas[f].r_pbc[3*p+2] < 0 ) activate_PBC[2] = 1;		
			if( rsurf[3*cp[p]+2] + shift[2] + theVolumeFormulas[f].r_pbc[3*p+2] > Lz ) activate_PBC[2] = 1;		
		}

		int ioff=0;
		for( int iu = 0; iu < ni; iu++ )
		for( int iv = 0; iv < ni-iu; iv++, ioff++ )
		{
			double R[3] = {shift[0], shift[1], shift[2] };


			for( int p = 0; p < np; p++ )
			{
				R[0] += theVolumeFormulas[f].r_w[ioff*np+p] * (rsurf[cp[p]*3+0] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+0]); 
				R[1] += theVolumeFormulas[f].r_w[ioff*np+p] * (rsurf[cp[p]*3+1] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+1]); 
				R[2] += theVolumeFormulas[f].r_w[ioff*np+p] * (rsurf[cp[p]*3+2] + theVolumeFormulas[f].r_pbc[ioff*np*3+3*p+2]); 
			}

			coords[(iu*ni+iv)*3+0] = R[0];
			coords[(iu*ni+iv)*3+1] = R[1];
			coords[(iu*ni+iv)*3+2] = R[2];
		}
		
		double dv_dr[(ni*ni)*3];
		memset( dv_dr, 0, sizeof(double) * ni*ni*3 );

		for( int i = 0; i < ni; i++ )
		for( int j = 0; j < ni-i; j++ )
		{
			
			for( int pass = 0; pass < 2; pass++ )
			{
				double *r1;
				double *r2;
				double *r3;

				int ir1, ir2, ir3;
				if( pass == 0 && (i + j + 1 < ni && i < ni-1 && j < ni-1) )
				{
					ir1 = i*ni+j;
					ir2 = ((i+1)*ni+j);
					ir3 = i*ni+j+1;
					r1 = coords + ir1*3;
					r2 = coords + ir2*3;
					r3 = coords + ir3*3;
				}
				else if( pass == 1 && (i-1 >= 0 && j+1 < ni && i+j+1 < ni) )
				{
					ir1 = i*ni+j;
					ir2 = i*ni+j+1;
					ir3 = (i-1)*ni+j+1;

					r1 = coords + ir1*3;
					r2 = coords + ir2*3;
					r3 = coords + ir3*3;
				}
				else
					continue;

				double dnrm[27]; // dnrmx, d r1-3_{x,y,z}
				double dA[9];

				if( dvol )
				{	
					double dr12[3] = { r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
					double dr32[3] = { r3[0]-r1[0],r3[1]-r1[1],r3[2]-r1[2]};
					double unrm[3];
					cross( dr12, dr32, unrm );
	
					memset( dnrm, 0, sizeof(double) * 27 );

					// derivative of a cross product with a vector wrt the other: 

					double d_cp1_dr[9], d_cp2_dr[9], d_cp3_dr[9];

					fillcp( d_cp1_dr, r1 );
					fillcp( d_cp2_dr, r2 );
					fillcp( d_cp3_dr, r3 );

					// (r2-r1) x ( r3-r1) 
					//  r2 x r3 - r1 x r3 - r2 x r1 

					// first, fill the derivative wrt the numerator.
					// r2 x r3 
					fillcpder( dnrm, d_cp3_dr, 1, +1 );
					fillcpder( dnrm, d_cp2_dr, 2, -1 ); //right side cross product gets neg

					// - r1 x r3
					fillcpder( dnrm, d_cp3_dr, 0, -1 ); 
					fillcpder( dnrm, d_cp1_dr, 2, +1 ); 
					
					// - r2 x r1
					fillcpder( dnrm, d_cp1_dr, 1, -1 );
					fillcpder( dnrm, d_cp2_dr, 0, +1 );

					// that's the numerator (the derivative of unrm)
					double d_nrm_len[9];
					memset( d_nrm_len, 0, sizeof(double) * 9 );
					
					double l = unrm[0]*unrm[0]+unrm[1]*unrm[1]+unrm[2]*unrm[2];

					// d_nrm_len is the derivative of unrm . unrm
					for( int p = 0; p < 3; p++ )
					for( int c = 0; c < 3; c++ )
					{
						d_nrm_len[p*3+c] = 2 * unrm[0] * dnrm[0 +3*p+c] +
							           2 * unrm[1] * dnrm[9 +3*p+c] +	 
							           2 * unrm[2] * dnrm[18+3*p+c];

						dA[p*3+c] = 0.5*d_nrm_len[p*3+c]/(2*sqrt(l));
					}

					// the normalization factor:
					for( int t = 0; t < 27; t++ )
						dnrm[t] /= sqrt(l);
					// that completes the derivative of the numerator.


					// denominator:
					for( int nc = 0; nc < 3; nc++ )
					for( int pc = 0; pc < 3; pc++ )
					for( int p  = 0; p  < 3; p++ )
					{
						// nc: the cartesian comp of the norm.
						//  p: which vector
						// pc: cartesian comp of that vector
						dnrm[nc*9+p*3+pc] +=  unrm[nc] * (-1.0/2.0) * pow(l,-3.0/2.0) * d_nrm_len[p*3+pc];	 
					}
				}
					
				double dr1[3] = { r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};
				double dr2[3] = { r3[0]-r1[0],r3[1]-r1[1],r3[2]-r1[2]};
				
				double nrm[3];
				cross( dr1, dr2, nrm );

				double A = normalize(nrm)/2;

				// this is the normal dotted into an arbitrary point on the triangle.
				double k = (r1[0] * nrm[0] + r1[1] * nrm[1] + r1[2] * nrm[2]);

				//vol += (1.0/3.0) * (bc[0] * nrm[0] + bc[1] * nrm[1] + bc[2] * nrm[2]) * A; 
				
				/*			
					if our triangle intersects the periodic edge we must wrap it inside.
				*/

				int n_sub_tri = 1;

				int sub_act = 0;

				if( r1[0] < 0 || r1[0] >= Lx ) sub_act = 1;
				if( r1[1] < 0 || r1[1] >= Ly ) sub_act = 1;
				if( r1[2] < 0 || r1[2] >= Lz ) sub_act = 1;
				
				if( r2[0] < 0 || r2[0] >= Lx ) sub_act = 1;
				if( r2[1] < 0 || r2[1] >= Ly ) sub_act = 1;
				if( r2[2] < 0 || r2[2] >= Lz ) sub_act = 1;
				
				if( r3[0] < 0 || r3[0] >= Lx ) sub_act = 1;
				if( r3[1] < 0 || r3[1] >= Ly ) sub_act = 1;
				if( r3[2] < 0 || r3[2] >= Lz ) sub_act = 1;
					
				int ir_offs[3] = { ir1, ir2, ir3 };

				if( sub_act ) 
				{	
					double verts[9]={r1[0],r1[1],r1[2],r2[0],r2[1],r2[2],r3[0],r3[1],r3[2]};

					double vol_temp[9];
					memset( vol_temp, 0, sizeof(double) * 9 );


					*vol_inside += splitTrianglesOnPBC(verts , nrm,
						Lx, Ly, Lz, &pbc_edge_data, (dvol ? vol_temp : NULL ), dnrm );

#ifdef DEBUG_VOLUME
					if( debug_trigger )
					{
						splitTrianglesOnPBC(verts , nrm,
							Lx, Ly, Lz, &pbc_edge_data, (dvol ? vol_temp : NULL ), dnrm );
						debug_trigger = 0;
					}
#endif
						if( dvol )
					{
						for( int p = 0; p < 3; p++ )
						for( int c = 0; c < 3; c++ )
							dv_dr[ir_offs[p]*3+c] += vol_temp[p*3+c];
					}
				}
				else
				{
#ifdef DEBUG_VOLUME
					fprintf(debug_triangles, "draw color orange3\n");
					fprintf(debug_triangles, "draw triangle { %lf %lf %lf } { %lf %lf %lf } { %lf %lf %lf }\n",
						r1[0], r1[1], r1[2],
						r2[0], r2[1], r2[2],
						r3[0], r3[1], r3[2] );
					vol_q[13] += k * A / 3.0;
#endif


#ifndef DISABLE_ALL
					*vol_inside  += k * A / 3.0;				
					if( dvol )
					{
						for( int p = 0; p < 3; p++ )
						for( int c = 0; c < 3; c++ )
						{
							// derivative wrt area.
							dv_dr[ir_offs[p]*3+c] += dA[3*p+c] * k / 3;
							// derivative wrt nrm
							for( int nc = 0; nc < 3; nc++ )
								dv_dr[ir_offs[p]*3+c] += A * r1[nc] * dnrm[nc*9+3*p+c] / 3;
						} 

						// derivative wrt point on surface
						for( int  c = 0; c < 3; c++ )
							dv_dr[ir_offs[0]*3+c] += A * nrm[c] / 3;
					}
#endif
				}
			}
		}	 

		if( dvol )
		{
			ioff=0;
			for( int iu = 0; iu < ni; iu++ )
			for( int iv = 0; iv < ni-iu; iv++, ioff++ )
			{
				int *cp = theVolumeFormulas[f].cp;
				int np = theVolumeFormulas[f].ncoor;
	
				for( int p = 0; p < np; p++ )
				{
					dvol[cp[p]*3+0] += dv_dr[(iu*ni+iv)*3+0] * theVolumeFormulas[f].r_w[ioff*np+p] * alpha[0];
					dvol[cp[p]*3+1] += dv_dr[(iu*ni+iv)*3+1] * theVolumeFormulas[f].r_w[ioff*np+p] * alpha[1];
					dvol[cp[p]*3+2] += dv_dr[(iu*ni+iv)*3+2] * theVolumeFormulas[f].r_w[ioff*np+p] * alpha[2];
				}
			}
		}
	} 
	
#ifdef DEBUG_VOLUME
	for( int q = 0; q < 27; q++ )
	{
		printf("%d %d %d vol %.14le\n", -1+(q/9), -1+(q/3)%3, -1+q%3, vol_q[q] );
	}
	printf("vol_inside: %le\n", *vol_inside );
#endif

	for( int edge = 0; edge < 12; edge++ )
	{
		if( edge != 3 && edge != 7 && edge != 11 )
			continue;

		if( pbc_edge_data.cur_edge_pointing_pos[edge] < pbc_edge_data.cur_edge_pointing_neg[edge] && pbc_edge_data.cur_edge_pointing_pos[edge] <= 1 )
		{
			*vol_inside += Lx*Ly*Lz/3; 
		}
	}


//	printf("vol inside: %lf\n", *vol_inside );
	*vol_outside = Lx*Ly*Lz - *vol_inside;

#ifdef DEBUG_VOLUME
	printf("Target area_A: %lf %lf %lf %lf %lf %lf\n", target_area_A[0], target_area_A[1], target_area_A[2], target_area_A[3], target_area_A[4], target_area_A[5] );
	printf("Target area_B: %lf %lf %lf %lf %lf %lf\n", target_area_B[0], target_area_B[1], target_area_B[2], target_area_B[3], target_area_B[4], target_area_B[5] );
	fclose(debug_triangles);
	
	FILE *thePSF = fopen("edge_pts.psf","w");
	writePSF( thePSF, nat, NULL, edge_bonds, nbonds );	
	fclose(thePSF);
	free(edge_bonds);
	fclose(edgePts);

#endif	

//	return vol;	
}


struct point_info 
{	
	double r[3]; // the real-space coords
	double uv[2]; // the u,v coords
	int cells[27]; // this point is a part of which cells.
	int ncells;
	double dp_dr[27]; // the derivative of the point's coordinates wrt the source.
	int intersect_type;
	// here we need enough information to compute derivatives wrt the vertices.
	//double drdc1
};

void swap (double*a,double*b)
{
	double t = *a;
	*a = *b;
	*b = t;
}

int build_out( 
		point_info *points,
		int *nuv,
		double *intersect_m, double *intersect_p, double *sign_move, int traverse_type, int cur_cell_in, double *r1, double *r2, double *r3, double *Ls );

double splitTrianglesOnPBC( double *verts, double *nrm, 
		double Lx, double Ly, double Lz, edge_data *pbc_edge_data, double *dvol, double *dnrm ) //, double *face_factors, double *spec_points, int *nspec)
{
	// this splits up a triangle into sub-triangles that are within the central PBC cell.
	// if a triangle is outside of the central face, we dot the normal into that periodic displacement
	// face factors keeps track of the area of the periodic face that we are intersecting.
	// we also need to know where on the PBC cell edge we intersect so we can figure out the area of the face intersection.
	
	// eq:
	// r1 + t1 (r1-r2) + t2 (r3-r1)
	double Ls[3] = { Lx, Ly, Lz };
	double k = nrm[0] * verts[0] + nrm[1] * verts[1] + nrm[2] * verts[2];	
				


	double r1[3] = { verts[0], verts[1], verts[2] };	
	double r2[3] = { verts[3], verts[4], verts[5] };	
	double r3[3] = { verts[6], verts[7], verts[8] };	


	double dr_1[3] = { r2[0]-r1[0], r2[1]-r1[1], r2[2]-r1[2] };	
	double dr_2[3] = { r3[0]-r2[0], r3[1]-r2[1], r3[2]-r2[2] };
	double dr_3[3] = { r1[0]-r3[0], r1[1]-r3[1], r1[2]-r3[2] };
	
	/*
	for( int c = 0; c < 3; c++ )
	{
		if( fabs(r1[c]/dr_1[c]) < HARD_EPS ) r1[c] += HARD_EPS * dr_1[c] * (r1[c]*dr_1[c])/(HARD_EPS+fabs(r1[c]*dr_1[c]))*2;
		if( fabs((r1[c]-Ls[c])/dr_1[c]) < HARD_EPS ) r1[c] += HARD_EPS * dr_1[c] * ((r1[c]-Ls[c])*dr_1[c])/(HARD_EPS+fabs( (r1[c]-Ls[c])*dr_1[c]))*2;
		
		if( fabs(r2[c]/dr_2[c]) < HARD_EPS ) r2[c] += HARD_EPS * dr_2[c] * (r2[c]*dr_2[c])/(HARD_EPS+fabs(r2[c]*dr_2[c]))*2;
		if( fabs((r2[c]-Ls[c])/dr_2[c]) < HARD_EPS ) r2[c] += HARD_EPS * dr_2[c] * ((r2[c]-Ls[c])*dr_2[c])/(HARD_EPS+fabs((r2[c]-Ls[c])*dr_2[c]))*2;
		
		if( fabs(r3[c]/dr_3[c]) < HARD_EPS ) r3[c] += HARD_EPS * dr_3[c] * (r3[c]*dr_3[c])/(HARD_EPS+fabs(r3[c]*dr_3[c]))*2;
		if( fabs((r3[c]-Ls[c])/dr_3[c]) < HARD_EPS ) r3[c] += HARD_EPS * dr_3[c] * ((r3[c]-Ls[c])*dr_3[c])/(HARD_EPS+fabs((r3[c]-Ls[c])*dr_3[c]))*2;		
	}*/
 
	

	// use a right-hand rule to determine segment traversal
	double vec1[3] = { dr_1[0], dr_1[1], dr_1[2] };
	double vec2[3] = { dr_2[0], dr_2[1], dr_2[2] };
	double check[3];
	cross( vec1, vec2, check );

	double sign = check[0] * nrm[0] + check[1] * nrm[1] + check[2] * nrm[2];
	// dot cross with vec1 into nrm to determine direction.

	// the t values at which this edge intersects an edge
	double t_u_intersect_m[3] = { -r1[0] / dr_1[0], -r1[1] / dr_1[1], -r1[2] / dr_1[2] }; 	
	double t_u_intersect_p[3] = { (Lx-r1[0]) / dr_1[0], (Ly-r1[1]) / dr_1[1], (Lz-r1[2]) / dr_1[2] }; 	
	double t_w_intersect_m[3] = { -r2[0] / dr_2[0], -r2[1] / dr_2[1], -r2[2] / dr_2[2] }; 	
	double t_w_intersect_p[3] = { (Lx-r2[0]) / dr_2[0], (Ly-r2[1]) / dr_2[1], (Lz-r2[2]) / dr_2[2] }; 	
	double t_v_intersect_m[3] = { -r3[0] / dr_3[0], -r3[1] / dr_3[1], -r3[2] / dr_3[2] }; 	
	double t_v_intersect_p[3] = { (Lx-r3[0]) / dr_3[0], (Ly-r3[1]) / dr_3[1], (Lz-r3[2]) / dr_3[2] }; 	

	double sign_u[3] = { 1,1,1};
	double sign_v[3] = { 1,1,1};
	double sign_w[3] = { 1,1,1};
	
	if( dr_1[0] < 0 ) sign_u[0] = -1;
	if( dr_1[1] < 0 ) sign_u[1] = -1; 
	if( dr_1[2] < 0 ) sign_u[2] = -1;
	
	if( dr_2[0] < 0 ) sign_w[0] = -1;
	if( dr_2[1] < 0 ) sign_w[1] = -1; 
	if( dr_2[2] < 0 ) sign_w[2] = -1;
	
	if( dr_3[0] < 0 ) sign_v[0] = -1;
	if( dr_3[1] < 0 ) sign_v[1] = -1; 
	if( dr_3[2] < 0 ) sign_v[2] = -1;


	point_info points[30];	

	int nuv = 0;

	int cur_cell_3[3] = { 1,1,1};

	if( r1[0] < 0 ) cur_cell_3[0] -= 1;
	if( r1[0] > Lx ) cur_cell_3[0] += 1;
	if( r1[1] < 0 ) cur_cell_3[1] -= 1;
	if( r1[1] > Ly ) cur_cell_3[1] += 1;
	if( r1[2] < 0 ) cur_cell_3[2] -= 1;
	if( r1[2] > Lz ) cur_cell_3[2] += 1;
	
	int cur_cell = cur_cell_3[0]*9+cur_cell_3[1]*3+cur_cell_3[2];

	points[nuv].uv[0] = 0;
	points[nuv].uv[1] = 0;
	points[nuv].intersect_type = 0;
	points[nuv].cells[0] = cur_cell; // in the zero-cell.
	points[nuv].ncells = 1;
	memset( points[nuv].dp_dr, 0, sizeof(double) * 27 );
	points[nuv].dp_dr[0   ]  = 1; // dx/dr1x
	points[nuv].dp_dr[9 +1]  = 1; // dy/dr1y
	points[nuv].dp_dr[18+2]  = 1; // dz/dr1z
	
	nuv++;

	cur_cell = build_out( points,
		&nuv, t_u_intersect_m, t_u_intersect_p, sign_u, 0, cur_cell,r1,r2,r3,Ls );
	
	points[nuv].uv[0] = 1;
	points[nuv].uv[1] = 0;
	points[nuv].intersect_type = 0;
	points[nuv].cells[0] = cur_cell; // in the zero-cell.
	points[nuv].ncells = 1;
	
	memset( points[nuv].dp_dr, 0, sizeof(double) * 27 );
	points[nuv].dp_dr[0 +3  ]  = 1; // dx/dr2x
	points[nuv].dp_dr[9 +3+1]  = 1; // dy/dr2y
	points[nuv].dp_dr[18+3+2]  = 1; // dz/dr2z
	nuv++;

	cur_cell = build_out( points,
			&nuv, t_w_intersect_m, t_w_intersect_p, sign_w, 2, cur_cell,r1,r2,r3,Ls );
	
	points[nuv].uv[0] = 0;
	points[nuv].uv[1] = 1;
	points[nuv].intersect_type = 0;
	points[nuv].cells[0] = cur_cell; // in the zero-cell.
	points[nuv].ncells = 1;
	
	memset( points[nuv].dp_dr, 0, sizeof(double) * 27 );
	points[nuv].dp_dr[0 +6  ]  = 1; // dx/dr3x
	points[nuv].dp_dr[9 +6+1]  = 1; // dy/dr3y
	points[nuv].dp_dr[18+6+2]  = 1; // dz/dr3z
	nuv++;
	
	cur_cell = build_out( points,
			&nuv, t_v_intersect_m, t_v_intersect_p, sign_v, 1, cur_cell,r1,r2,r3,Ls );

	// count  the interior points: they cross if 
	
	int prev_type = -1;

	double intersect[3][3] =
	{
		{0,0,0},
		{0,0,0},
		{0,0,0}
	};
	
	double intersect_pts[3][3][2] =
	{
		{ {-1,-1},{-1,-1},{-1,-1} },
		{ {-1,-1},{-1,-1},{-1,-1} },
		{ {-1,-1},{-1,-1},{-1,-1} }
	};

	int checked[3] = {0,0,0};
	int axis_sign[3] = { 0,0,0 };
	double int_axes[3][4];
	double int_axes_pts[3][2];
	
	for( int uv = 0; uv < nuv; uv++ )
	{
		if( points[uv].intersect_type != 0)
		{
			int axis = points[uv].intersect_type;
			if( axis < 0 ) axis *= -1;

			if(checked[axis-1] ) continue;

			int_axes[axis-1][0] = points[uv].uv[0]; 
			int_axes[axis-1][1] = points[uv].uv[1]; 

			for( int uv2 = uv+1; uv2 < nuv; uv2++ )
			{
				if( points[uv2].intersect_type == points[uv].intersect_type )
				{
					int_axes[axis-1][2] = points[uv2].uv[0]; 
					int_axes[axis-1][3] = points[uv2].uv[1]; 

					int_axes_pts[axis-1][0] = uv;
					int_axes_pts[axis-1][1] = uv2;
				}
			}

			checked[axis-1] = 1;
		}
	}
	
	checked[0] = 0;
	checked[1] = 0;
	checked[2] = 0;

	for( int uv = 0; uv < nuv; uv++ )
	{
		if( points[uv].intersect_type != 0)
		{
			int axis = points[uv].intersect_type;
			if( axis < 0 ) axis *= -1;
			if( checked[axis-1] ) continue;

			if( points[uv].intersect_type < 0 )
				axis_sign[axis-1] = -1;
			else
				axis_sign[axis-1] = 1;

			intersect[axis-1][0] = 0;		
			intersect[axis-1][1] = 0;		
			intersect[axis-1][2] = 0;		

			int done = 0;
			int uv2 = uv+1;

			while( ! done )
			{
				if( uv2 >= nuv ) uv2 = 0;
				if( uv2 == uv ) break;

				if( points[uv].intersect_type == points[uv2].intersect_type )
					break;	

				if( points[uv2].intersect_type != 0 )
				{
					int axis2 = points[uv2].intersect_type;
					if( axis2 < 0 ) axis2 *= -1;

					intersect[axis-1][axis2-1] = !intersect[axis-1][axis2-1];
					intersect_pts[axis-1][axis2-1][0] = uv;
					intersect_pts[axis-1][axis2-1][1] = uv2;

				
				}

				uv2++;
			}

			checked[axis] = 1;
		}
	}



//	*nspec = 0;


	for( int axis1 = 0; axis1 < 3; axis1++ )
	for( int axis2 = axis1+1; axis2 < 3; axis2++ )
	{
		if( intersect[axis1][axis2] )
		{
			double a[2] = { int_axes[axis1][0], int_axes[axis1][1] };
			double da[2] = { int_axes[axis1][2]-a[0], int_axes[axis1][3]-a[1] };
			double b[2] = { int_axes[axis2][0], int_axes[axis2][1] };
			double db[2] = { int_axes[axis2][2]-b[0], int_axes[axis2][3]-b[1] };
	
			double denom=da[1]*db[0]-da[0] * db[1];
			double t1 = -(a[1]*db[0]-b[1]*db[0]-a[0]*db[1]+b[0]*db[1])/denom;
			double isect[2] = { a[0] + t1 * da[0],
					    a[1] + t1 * da[1] };
		

			points[nuv].uv[0] = isect[0];
			points[nuv].uv[1] = isect[1];
			points[nuv].intersect_type = -1;
			points[nuv].ncells = 0;
			
			// triangle holds the intersection of two axes. 


			int uv = intersect_pts[axis1][axis2][0];
			int uv2 = intersect_pts[axis1][axis2][1];

			

#if 0	 // this doesn't work. it's wrong.
			for( int axis = 0; axis < 3; axis++ )
			for( int t_p = 0; t_p < 2; t_p ++ )
			{
				if( axis != axis1 && axis != axis2 ) continue;

				int uv = int_axes_pts[axis][t_p];

				for( int t = 0; t < points[uv].ncells; t++ )
				{
					int gotit = 0;

					for( int t2 = 0; t2 < points[nuv].ncells; t2++ )
					{
						if( points[nuv].cells[t2] == points[uv].cells[t] ) { gotit = 1; break; }
					}

					if( !gotit )
					{
						points[nuv].cells[points[nuv].ncells] = points[uv].cells[t];
						points[nuv].ncells++;
					}
				}
			}
#endif

			double tr[3] = { 
r1[0] + isect[0] * (r2[0]-r1[0]) +isect[1] * (r3[0]-r1[0]),
r1[1] + isect[0] * (r2[1]-r1[1]) +isect[1] * (r3[1]-r1[1]),
r1[2] + isect[0] * (r2[2]-r1[2]) +isect[1] * (r3[2]-r1[2])
			};
			
			// BEGIN DERIVATIVE
			{
				double k1 = 0;
				if( fabs(tr[axis1]-Ls[axis1]) < HARD_EPS )
					k1 = Ls[axis1];
				double k2 = 0;
				if( fabs(tr[axis2]-Ls[axis2]) < HARD_EPS )
					k2 = Ls[axis2];
	
				// not x y, axis1 or axis2.
				double x0 =r1[axis1];
				double y0 =r1[axis2];
				double dx1 = r2[axis1]-r1[axis1];
				double dx2 = r3[axis1]-r1[axis1];
				double dy1 = r2[axis2]-r1[axis2];
				double dy2 = r3[axis2]-r1[axis2];
	
				double t1 = -(dy2*k1-dx2*k2-dy2*x0+dx2*y0)/(dx2*dy1-dx1*dy2);
				double t2 = -(-dy1*k1+dx1*k2+dy1*x0-dx1*y0)/(dx2*dy1-dx1*dy2);
	
				double d_t1_x0 = dy2/(dx2*dy1-dx1*dy2);
				double d_t1_dx1 = -dy2*(dy2*k1-dx2*k2-dy2*x0+dx2*y0)/(dx2*dy1-dx1*dy2)/(dx2*dy1-dx1*dy2);
				double d_t1_dx2 = -((-k2 + y0)/(dx2*dy1 - dx1*dy2)) +  (dy1*(dy2*k1 - dx2*k2 - dy2*x0 + dx2*y0))/pow(dx2*dy1 - dx1*dy2,2.);
				
				double d_t1_y0 = -dx2/(dx2*dy1-dx1*dy2); 
				double d_t1_dy1 = (dx2*(dy2*k1 - dx2*k2 - dy2*x0 + dx2*y0))/pow(dx2*dy1 - dx1*dy2,2.);
				double d_t1_dy2 = -((k1 - x0)/(dx2*dy1 - dx1*dy2)) - (dx1*(dy2*k1 - dx2*k2 - dy2*x0 + dx2*y0))/pow(dx2*dy1 - dx1*dy2,2.);
	
				double d_t2_x0 = -(dy1/(dx2*dy1 - dx1*dy2));
				double d_t2_dx1 = -((k2 - y0)/(dx2*dy1 - dx1*dy2)) - (dy2*(-(dy1*k1) + dx1*k2 + dy1*x0 - dx1*y0))/pow(dx2*dy1 - dx1*dy2,2.);
				double d_t2_dx2 = (dy1*(-(dy1*k1) + dx1*k2 + dy1*x0 - dx1*y0))/pow(dx2*dy1 - dx1*dy2,2.);
	
				double d_t2_y0 = dx1/(dx2*dy1 - dx1*dy2);
				double d_t2_dy1 = -((-k1 + x0)/(dx2*dy1 - dx1*dy2)) + (dx2*(-(dy1*k1) + dx1*k2 + dy1*x0 - dx1*y0))/pow(dx2*dy1 - dx1*dy2,2.);
				double d_t2_dy2 = -((dx1*(-(dy1*k1) + dx1*k2 + dy1*x0 - dx1*y0))/pow(dx2*dy1 - dx1*dy2,2.));
	
				memset( points[nuv].dp_dr,0,sizeof(double)*27);
			
				for( int dc = 0; dc < 3; dc++ ) // cartesian component of target vector		
				for( int c  = 0; c < 3; c++ ) // cartesian component of der vector
				{
					if( c == dc )
					{	
						points[nuv].dp_dr[dc*9+0*3+c] = 1-t1-t2;
						points[nuv].dp_dr[dc*9+1*3+c] = t1;
						points[nuv].dp_dr[dc*9+2*3+c] = t2;
					}
					// d t d r1x
			
					if( c == axis1 )
					{
						// derivatives of t1/t2 wrt the base position, first cartesian axis
						points[nuv].dp_dr[dc*9+0*3+c] += d_t1_x0 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] += d_t2_x0 * (r3[dc]-r1[dc]);
	
						// derivatives of t1/t2/ wrt the first dx, first cartesian axis
						
						points[nuv].dp_dr[dc*9+1*3+c] += d_t1_dx1 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+1*3+c] += d_t2_dx1 * (r3[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t1_dx1 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t2_dx1 * (r3[dc]-r1[dc]);
						
						points[nuv].dp_dr[dc*9+2*3+c] += d_t1_dx2 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+2*3+c] += d_t2_dx2 * (r3[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t1_dx2 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t2_dx2 * (r3[dc]-r1[dc]);
					}	
	
					if( c == axis2 )
					{
						// derivatives of t1/t2 wrt the base position, second cartesian axis
						points[nuv].dp_dr[dc*9+0*3+c] += d_t1_y0 * (r2[dc]-r1[dc]); 
						points[nuv].dp_dr[dc*9+0*3+c] += d_t2_y0 * (r3[dc]-r1[dc]); 
						
						points[nuv].dp_dr[dc*9+1*3+c] += d_t1_dy1 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+1*3+c] += d_t2_dy1 * (r3[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t1_dy1 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t2_dy1 * (r3[dc]-r1[dc]);
						
						points[nuv].dp_dr[dc*9+2*3+c] += d_t1_dy2 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+2*3+c] += d_t2_dy2 * (r3[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t1_dy2 * (r2[dc]-r1[dc]);
						points[nuv].dp_dr[dc*9+0*3+c] -= d_t2_dy2 * (r3[dc]-r1[dc]);
					}	
				}
			}	
			// END DERIVATIVE
				
			//spec_points[(*nspec)*3+0] = tr[0];
			//spec_points[(*nspec)*3+1] = tr[1];
			//spec_points[(*nspec)*3+2] = tr[2];

			int pbc_type[3] = { 0,0,0 };

			for( int c = 0; c < 3; c++ )
			{
				if( fabs(tr[c]) < HARD_EPS) 
					pbc_type[c] = -1;
				else if( fabs(tr[c]-Ls[c]) < HARD_EPS )
					pbc_type[c] = 1;
				else
					pbc_type[c] = 0;
			}

			int base_cell[3] = { 1,1,1};

			for( int cx = -1; cx <= 1; cx++ )
			for( int cy = -1; cy <= 1; cy++ )
			for( int cz = -1; cz <= 1; cz++ )
			{
				if( !( (tr[0] < HARD_EPS && cx == -1 ) || (tr[0] > -HARD_EPS && tr[0] <= Lx+HARD_EPS && cx == 0) || (tr[0] >= Lx-HARD_EPS && cx == 1 ) ) ) continue;
				if( !( (tr[1] < HARD_EPS && cy == -1 ) || (tr[1] > -HARD_EPS && tr[1] <= Ly+HARD_EPS && cy == 0) || (tr[1] >= Ly-HARD_EPS && cy == 1 ) ) ) continue;
				if( !( (tr[2] < HARD_EPS && cz == -1 ) || (tr[2] > -HARD_EPS && tr[2] <= Lz+HARD_EPS && cz == 0) || (tr[2] >= Lz-HARD_EPS && cz == 1 ) ) ) continue;

				int the_cell = (cx+1) * 9 + (cy+1) * 3 + cz + 1;

				points[nuv].cells[points[nuv].ncells] = the_cell;
				points[nuv].ncells++;
			}
			nuv++;
//			(*nspec)++;
		}			
	}

#ifdef DEBUG_VOLUME
	const char *colors[27] =
	{
"mauve",
"red",
"gray",
"orange",
"yellow",
"tan",
"silver",
"green",
"white",
"pink",
"cyan",
"purple",
"lime",
"blue",//13
"ochre",
"iceblue",
"magenta",
"yellow2",
"yellow3",
"green2",
"green3",
"cyan2",
"cyan3",
"blue2",
"blue3",
"violet",
"violet2"
	};

#endif

	double the_volume = 0;
	int altc[3][3] = 
	{
		{ -1, 2, 1 },
		{  2, -1, 0 },
		{  1, 0,  -1 }
	};
	
	for( int p = 0; p < nuv; p++ )
	{
		points[p].r[0] = r1[0] + points[p].uv[0] * (r2[0]-r1[0]) + points[p].uv[1] * (r3[0]-r1[0]);
		points[p].r[1] = r1[1] + points[p].uv[0] * (r2[1]-r1[1]) + points[p].uv[1] * (r3[1]-r1[1]);
		points[p].r[2] = r1[2] + points[p].uv[0] * (r2[2]-r1[2]) + points[p].uv[1] * (r3[2]-r1[2]);
	}

// temporary volume debug
#if 0 
	for( int u = 0; u < nuv; u++ )
	{
		if( points[u].ncells > 2)
		{
			if( nbuild_pts_tallied == ntally_target )
			{ 
				printf("Target pt : %le %le %le\n", points[u].r[0], points[u].r[1], points[u].r[2] );
				int ic = 2; // the icth cartesian-component of the build point
				if( dvol )
					memcpy( dvol, points[u].dp_dr+ic*9, sizeof(double) * 9 );
				
				nbuild_pts_tallied++;	
				return points[u].r[ic];	
			}
			nbuild_pts_tallied++;	

		}
	}
	return 0;
#endif
	int drew = 0;
	if(nuv == 3 )
	{
		int q = points[0].cells[0];

		double cell[3] = { q / 9, (q/3)%3, q %3 };

		// the shift to bring the element into the unit cell.
		double shift[3] = { (1-cell[0]) * Lx, (1-cell[1]) * Ly, (1-cell[2]) * Lz }; 
#ifdef DEBUG_VOLUME
		fprintf(debug_triangles, "draw color %s\n", colors[13] );
		fprintf(debug_triangles, "draw triangle { %lf %lf %lf } { %lf %lf %lf } { %lf %lf %lf }\n",
			r1[0]+shift[0], r1[1]+shift[1], r1[2]+shift[2],
			r2[0]+shift[0], r2[1]+shift[1], r2[2]+shift[2],
			r3[0]+shift[0], r3[1]+shift[1], r3[2]+shift[2]
		);
		drew = 1;
#endif
		double ordered_pts[9] = {
			r1[0] + shift[0], r1[1] + shift[1], r1[2] + shift[2],
			r2[0] + shift[0], r2[1] + shift[1], r2[2] + shift[2],
			r3[0] + shift[0], r3[1] + shift[1], r3[2] + shift[2],
		};


#ifdef DEBUG_VOLUME
		for( int p = 0; p < 3; p++ )
		{
			if( ordered_pts[3*p+0] < -HARD_EPS || ordered_pts[3*p+0] > Ls[0] + HARD_EPS ||	
			    ordered_pts[3*p+1] < -HARD_EPS || ordered_pts[3*p+1] > Ls[1] + HARD_EPS ||
			    ordered_pts[3*p+2] < -HARD_EPS || ordered_pts[3*p+2] > Ls[2] + HARD_EPS )
				debug_trigger = 1;
		}
#endif			

		double cp[3];
		cross( dr_1, dr_2, cp );
		double face_area = normalize(cp)/2;
		double k = (nrm[0] * ordered_pts[0] + nrm[1] * ordered_pts[1] + nrm[2] * ordered_pts[2]);
#ifndef DISABLE_ALL	
		the_volume += (1.0/3.0) * face_area * k;
		if( dvol )
		{	// derivative wrt vol, [d ordered_pts]
			dvol[0] += (1.0/3.0) * face_area * nrm[0]; 
			dvol[1] += (1.0/3.0) * face_area * nrm[1]; 
			dvol[2] += (1.0/3.0) * face_area * nrm[2]; 
		
			// wrt nrm	
			for( int nc = 0; nc < 3; nc++ )
			for( int p = 0; p < 3; p++ )
			for( int c = 0; c < 3; c++ )
				dvol[3*p+c] += (1.0/3.0) * face_area * dnrm[nc*9+3*p+c] * ordered_pts[nc]; 
			
			// wrt area
			darea( dvol, r1, r2, r3,  (1.0/3.0) * k ); 
		}
#endif
#ifdef DEBUG_VOLUME
//		vol_q[q] += (1.0/3.0) * face_area * (nrm[0] * ordered_pts[0] + nrm[1] * ordered_pts[1] + nrm[2] * ordered_pts[2]);
		vol_q[q] += (1.0/3.0) * face_area * k;
#endif		
	}
	else if( 1 ) 
//	if( nuv > 5 ) // most of the work will be spent doing nuv == 3 or nuv == 5.
	{
		// more than two triangles, requires giftwrapping
		double quadrant_points[27][12];
		int quad_spec[27][6];
		int quad_src[27][6];
		int nquad[27];
		memset( nquad, 0, sizeof(int) * 27 );

		for( int p = 0; p < nuv; p++ )
		{
			for( int t = 0; t < points[p].ncells; t++ )
			{	
				int q = points[p].cells[t];

				if( q < 0 || q >= 27 )
				{
#ifdef DEBUG_VOLUME
					printf("PBC cell breach.\n");
					debug_trigger = 1;
#endif
					continue;
				}

				if( nquad[q] > 6 )
				{
				printf("Logical error in volume/pbc routines.\n");
				exit(1);
				}
				quadrant_points[q][nquad[q]*2+0] = points[p].uv[0];
				quadrant_points[q][nquad[q]*2+1] = points[p].uv[1];
				if( points[p].ncells > 2 )
					quad_spec[q][nquad[q]] = 2;
				else if( points[p].ncells == 2 )
					quad_spec[q][nquad[q]] = 1;
				else
					quad_spec[q][nquad[q]] = 0;

				quad_src[q][nquad[q]] = p;

				nquad[q] += 1;
			}
		}
		
		int ordered[12];	
		for( int q = 0; q < 27; q++ )
		{
			if( nquad[q] == 0 ) continue;

			if( nquad[q] > 6)
			{
				printf("Logical error in volume/pbc routines.\n");
				exit(1);
			}


/*			if( nquad[q] < 3 )
			{
				printf("Logical error in volume/pbc routines.\n");
				exit(1);
			}
*/		
			if( nquad[q] >= 3 )
			{
				double cell[3] = { q / 9, (q/3)%3, q %3 };

				// the shift to bring the element into the unit cell.
				double shift[3] = { (1-cell[0]) * Lx, (1-cell[1]) * Ly, (1-cell[2]) * Lz }; 
				
				int nconvex=0;

//				double ordered_points_uv[2*nquad[q]];
				int ordered[nquad[q]*2];
				giftwrap( quadrant_points[q], ordered, nquad[q], &nconvex, 1 /*do expand */ );



				double area = 0;

				double ordered_pts[nquad[q]*3];
				for( int t = 0; t < nconvex; t++ )
				{
					double u = quadrant_points[q][ordered[t]*2+0];
					double v = quadrant_points[q][ordered[t]*2+1];

					ordered_pts[t*3+0] = r1[0] + u * (r2[0]-r1[0]) + v * (r3[0]-r1[0]) + shift[0];
					ordered_pts[t*3+1] = r1[1] + u * (r2[1]-r1[1]) + v * (r3[1]-r1[1]) + shift[1];
					ordered_pts[t*3+2] = r1[2] + u * (r2[2]-r1[2]) + v * (r3[2]-r1[2]) + shift[2];

				}

#ifdef DEBUG_VOLUME
				for( int p = 0; p < nconvex; p++ )
				{
					if( ordered_pts[3*p+0] < -HARD_EPS || ordered_pts[3*p+0] > Ls[0] + HARD_EPS ||	
					    ordered_pts[3*p+1] < -HARD_EPS || ordered_pts[3*p+1] > Ls[1] + HARD_EPS ||
					    ordered_pts[3*p+2] < -HARD_EPS || ordered_pts[3*p+2] > Ls[2] + HARD_EPS )
						debug_trigger = 1;
				}
#endif			

						
				
				for( int qx = 0; qx < nconvex; qx++ )
				{
					int qxp1 = qx+1;
					if( qxp1 >= nconvex )
						qxp1 -= nconvex;
	
					if( quad_spec[q][ordered[qxp1]] > 1 )
					{
						// the following point is on the edge of the cube: it contributes to the area of two faces
						double *p = ordered_pts+qxp1*3;
						// this is the vector flowing into the point.
						double dr[3] = { p[0] - ordered_pts[qx*3+0],
								 p[1] - ordered_pts[qx*3+1],
								 p[2] - ordered_pts[qx*3+2] };
						

						int coor_edge = 0;

						for( int c = 0; c < 3; c++ )
						{
							if( fabs(p[c]) > HARD_EPS && fabs(p[c]-Ls[c]) > HARD_EPS )
								coor_edge = c;
						}
						dr[coor_edge] = 0;
						normalize(dr);

						int approach_face=0;
						double abs_dr = HARD_EPS;
						for( int c = 0; c < 3; c++ )
							if( fabs(dr[c]) > abs_dr )
							{
								abs_dr = fabs(dr[c]);
								approach_face = c;
							}
						if( approach_face == -1 )
						{
							printf("Logical error on approach face.\n");
							exit(1);
						}

						int coming_in = 1;
						if( dr[approach_face] > 0 && p[approach_face] < Ls[approach_face]/2 )
							coming_in = 0;
						if( dr[approach_face] < 0 && p[approach_face] > Ls[approach_face]/2 )
						 	coming_in = 0;

						double the_axis[3] = { 0,0,0};
						the_axis[coor_edge] = 1;
	
						double dcoor_edge = p[coor_edge];

						// if the + cell edge is involved, need to do p-L:

						if( nrm[coor_edge] < 0 ) // dotting into edge axis
							dcoor_edge -= Ls[coor_edge];
						dcoor_edge -= Ls[coor_edge];


						// for the vector for which we are "coming into the edge", dr[approach_face] != 0
						// for this we are the *first* point in the edge sequence.
						
						

						for( int c = 0; c < 3; c++ )
						for( int pm = 0; pm < 2; pm++ )
						{
							double face_normal[3] = { 0,0,0};
							if( pm == 0 )
							{
								face_normal[c] = -1;
	
								if( fabs(p[c]) > HARD_EPS ) // it is not on this face. 
									continue;
							}								
							else
							{
								face_normal[c] = 1;

								if( fabs(p[c]-Ls[c]) > HARD_EPS)
									continue;
							}
							double sense = 1.0;

							// does the surface normal point out or in?
					//		if( nrm[0] * the_axis[0] + nrm[1] * the_axis[1] + nrm[2] * the_axis[2] < 0 )
					//			sense *= -1;

							// edge indexing scheme defined in edge_data
							int edge_index = coor_edge * 4;
							int alt_edges[3][2] =
							{
								{1,2},	
								{0,2},
								{0,1}
							};
						
							if( fabs(p[alt_edges[coor_edge][0]]-Ls[alt_edges[coor_edge][0]])<HARD_EPS )
								edge_index += 2;
							if( fabs(p[alt_edges[coor_edge][1]]-Ls[alt_edges[coor_edge][1]])<HARD_EPS )
								edge_index += 1;

							double p1[3] = { p[0], p[1], p[2] };
							double p2[3] = { p[0], p[1], p[2] };

							int der_type = 0;


							// c is the Cartesian coordinate of the face we are computing.
							if( 
									(coming_in && c != approach_face ) // in this case, the vector is in-coming
								      ||(!coming_in && c == approach_face ) // in this case, we are computing for the other direction, and we are in-coming
									)
							{
								if( nrm[coor_edge] < 0 )
								{
									p2[coor_edge] = 0; // Ls[coor_edge];
									
									double x = 1-p[coor_edge]/Ls[coor_edge];
									if( x > pbc_edge_data->cur_edge_pointing_pos[edge_index] )
										pbc_edge_data->cur_edge_pointing_pos[edge_index] = x; 
								}
								else
								{
									p2[coor_edge] = 0; 
									
									double x = 1-p[coor_edge]/Ls[coor_edge];
									if( x > pbc_edge_data->cur_edge_pointing_neg[edge_index] )
										pbc_edge_data->cur_edge_pointing_neg[edge_index] = x; 
								}
								der_type = 0;
							}	
							else
							{ // we are the second point.
								if( nrm[coor_edge] < 0 ) 
								{
									// pointing up: first point is the PBC edge: may need L^2/2 factor
									p1[coor_edge] = 0; // Ls[coor_edge];
								
									double x = 1-p[coor_edge]/Ls[coor_edge];
									if( x < pbc_edge_data->cur_edge_pointing_pos[edge_index] )
										pbc_edge_data->cur_edge_pointing_pos[edge_index] = x; 
								}
								else
								{
									p1[coor_edge] = 0; 
									double x = 1-p[coor_edge]/Ls[coor_edge];
									if( x < pbc_edge_data->cur_edge_pointing_neg[edge_index] )
										pbc_edge_data->cur_edge_pointing_neg[edge_index] = x; 
								}
								der_type = 1;
							}	

							double crossp[3];
							cross(p1,p2,crossp);

#ifdef DEBUG_VOLUME
							if( face_normal[0] == 1 )
								target_area_B[0] += ( face_normal[c] * crossp[c] /2 ); 
							else if( face_normal[0] == -1 )
								target_area_B[1] += ( face_normal[c] * crossp[c] /2 ); 
							else if( face_normal[1] == 1 )
								target_area_B[2] += ( face_normal[c] * crossp[c] /2 ); 
							else if( face_normal[1] == -1 )
								target_area_B[3] += ( face_normal[c] * crossp[c] /2 ); 	
							else if( face_normal[2] == 1 )
								target_area_B[4] += ( face_normal[c] * crossp[c] /2 ); 
							else if( face_normal[2] == -1 )
								target_area_B[5] += ( face_normal[c] * crossp[c] /2 ); 
#endif
							double factor = (1.0/3.0)*(p[c] * face_normal[c]) * ( face_normal[c] * crossp[c] /2 );
							const char *faces="xyz";
//							printf("factor: %le Edge quadrant [%d %d %d] point [%lf %lf %lf] face %c nrm %lf %lf %lf sign: %lf\n", factor,
//							(q/9)-1,(q/3)%3-1,q%3-1, p[0],p[1],p[2], faces[c], nrm[0], nrm[1], nrm[2], sense );

							
#ifdef DEBUG_VOLUME
							vol_q[q] += factor;
#endif
#ifndef DISABLE_ALL	
							the_volume += factor;
							if( dvol )
							{
								
								int srcp = quad_src[q][ordered[qxp1]];	

								// I think this should be zero:

								// derivative wrt p[c]:
//								for( int tp = 0; tp < 3; tp++ )
//								for( int tc = 0; tc < 3; tc++ )
//									dvol[tp*3+tc] += (1.0/3.0)*(points[srcp].dp_dr[c*9+3*tp+tc] * face_normal[c]) * ( face_normal[c] * crossp[c] /2 );

								// the derivative of cross(p1,p2)
								// HOWEVER: only one variable of the cross-product varies, which one?
								double val = 1;
								double sign = 1;
								double the_sign[3][3] = 
								{
									{	0,	-1,	1 },
									{	1,	 0,	-1 },
									{	-1,	 1,	0 }	
								};
								int ider = altc[coor_edge][c];	
								// derivative wrt crossp[c]
								if( der_type == 0 )
									val = p2[ider];
								else
									val = -p1[ider];

								for( int tp = 0; tp < 3; tp++ )
								for( int tc = 0; tc < 3; tc++ )
									dvol[tp*3+tc] += (1.0/3.0)*(p[c] * face_normal[c]) * ( face_normal[c] * val * points[srcp].dp_dr[coor_edge*9+3*tp+tc] * the_sign[coor_edge][c] )/2;
							}
#endif
						}
					}

					
					if( quad_spec[q][ordered[qx]] && quad_spec[q][ordered[qxp1]] )
					{
						// segment falls on periodic boundary.
				
						double crossp[3];
						cross( ordered_pts+qx*3, ordered_pts+qxp1*3, crossp );

						double dr[3] = { 
							ordered_pts[qx*3+0] - ordered_pts[qxp1*3+0],
							ordered_pts[qx*3+1] - ordered_pts[qxp1*3+1],
							ordered_pts[qx*3+2] - ordered_pts[qxp1*3+2] };
				
	
						// this segment better be along an axis.
						int on_axis[3] = {0,0,0};
	
						for( int c = 0; c < 3; c++ )
						{
							if( fabs(dr[c]) < HARD_EPS ) // moving in the x=k plane
							{
								if( fabs(ordered_pts[qx*3+c]) < HARD_EPS ||
								    fabs(ordered_pts[qx*3+c]-Ls[c])<HARD_EPS)
									on_axis[c] =1;
							}
						}
						if( on_axis[0] || on_axis[1] || on_axis[2] )
						{		
							int iface = 0;	
							double face_normal[3] = { 0,0,0};
							if( fabs(dr[1]) < HARD_EPS && on_axis[1] )
								iface = 1;
							if( fabs(dr[2]) < HARD_EPS && on_axis[2] )
								iface = 2;

							if( ordered_pts[qx*3+iface] < Ls[iface]/2 )
								face_normal[iface] = -1;
							else	
								face_normal[iface] = 1;
					
							double x_vec1[3], x_vec2[3];
							cross( vec1, dr, x_vec1 );
							double dp1 = x_vec1[0] * nrm[0] + x_vec1[1] * nrm[1] + x_vec1[2] * nrm[2];
							double the_sign = 1;
	/*						if( dp1 < -EPS )
							{
								the_sign = -1;
							}
							else if( dp1 < EPS )
							{
								cross( vec2, dr, x_vec2 );
								double dp2 = x_vec2[0] * nrm[0] + x_vec2[1] * nrm[1] + x_vec2[2] * nrm[2];
								if( dp2 < 0 )
									the_sign = -1;
							}
	*/
#ifdef DEBUG_VOLUME
							fprintf( debug_triangles, "draw color %s\n", colors[q] );
							if( edge_space == nbonds )
							{
								edge_space *= 2;
								edge_bonds = (int *)realloc( edge_bonds, sizeof(int) * 2 * edge_space );
							}
							if( the_sign < 0 )
							{
								drawVMDArrow( debug_triangles, ordered_pts+3*qxp1, ordered_pts+3*qx );
								fprintf(edgePts, "C %lf %lf %lf\n", ordered_pts[3*qxp1+0], ordered_pts[3*qxp1+1], ordered_pts[3*qxp1+2] );
								fprintf(edgePts, "O %lf %lf %lf\n", ordered_pts[3*qx+0], ordered_pts[3*qx+1], ordered_pts[3*qx+2] );
								edge_bonds[2*nbonds+0] = nat;
								edge_bonds[2*nbonds+1] = nat+1;
								nat+= 2;
								nbonds++;
							}
							else //if( face_normal[2] == 1 )
							{
								drawVMDArrow( debug_triangles, ordered_pts+3*qx, ordered_pts+3*qxp1 );
								fprintf(edgePts, "C %lf %lf %lf\n", ordered_pts[3*qx+0], ordered_pts[3*qx+1], ordered_pts[3*qx+2] );
								fprintf(edgePts, "O %lf %lf %lf\n", ordered_pts[3*qxp1+0], ordered_pts[3*qxp1+1], ordered_pts[3*qxp1+2] );
								edge_bonds[2*nbonds+0] = nat;
								edge_bonds[2*nbonds+1] = nat+1;
								nat+= 2;
								nbonds++;
							}
#endif
	
#ifdef DEBUG_VOLUME
							if( face_normal[0] == 1 )
								target_area_A[0] += (crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2; 
							else if( face_normal[0] == -1 )
								target_area_A[1] += (crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2; 
							else if( face_normal[1] == 1 )
								target_area_A[2] += (crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2; 
							else if( face_normal[1] == -1 )
								target_area_A[3] += (crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2; 
							else if( face_normal[2] == 1 )
							{
								printf("pta: %le %le %le\nptb: %le %le %le contrib %le\n", ordered_pts[3*qx+0], ordered_pts[3*qx+1], ordered_pts[3*qx+2],
									ordered_pts[3*qxp1+0], ordered_pts[3*qxp1+1], ordered_pts[3*qxp1+2], (1.0/3.0) * (face_normal[0]*ordered_pts[3*qx+0] + face_normal[1] * ordered_pts[3*qx+1] + face_normal[2]*ordered_pts[3*qx+2]) *(crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2 ); 
								target_area_A[4] += (crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2; 
							}
							else if( face_normal[2] == -1 )
								target_area_A[5] += (crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2; 
#endif


#ifndef DISABLE_ALL	
							the_volume += (1.0/3.0) * // this determines the ``sense'' of how we traverse to get the PBC face area. 
								// this is the dot product of the normal with a point on the face.
								(face_normal[0]*ordered_pts[3*qx+0] + face_normal[1] * ordered_pts[3*qx+1] + face_normal[2]*ordered_pts[3*qx+2]) * 
								(crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2; // contribution to area; 					
							if( dvol )
							{	
								int srcp1 = quad_src[q][ordered[qx]];
								int srcp2 = quad_src[q][ordered[qxp1]];

								// derivative wrt ordered_pts[3*qx]... this should be zero:
								for( int pc = 0; pc < 3; pc++ )
								for( int s = 0; s < 3; s++ )
								for( int sc = 0; sc < 3; sc++ )
								{
									dvol[3*s+sc] += (1.0/3.0) * face_normal[pc] * points[srcp1].dp_dr[9*pc+3*s+sc] * (crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2;
								}
	
								double prefac = (1.0/3.0) *(face_normal[0]*ordered_pts[3*qx+0] + face_normal[1] * ordered_pts[3*qx+1] + face_normal[2]*ordered_pts[3*qx+2]) / 2.0;
	
								dcross( dvol, ordered_pts+qx*3, ordered_pts+qxp1*3, points[srcp1].dp_dr, points[srcp2].dp_dr, prefac, face_normal );						
							}
#endif

#ifdef DEBUG_VOLUME
							vol_q[q] += (1.0/3.0)* 
								(face_normal[0]*ordered_pts[3*qx+0] + face_normal[1] * ordered_pts[3*qx+1] + face_normal[2]*ordered_pts[3*qx+2])* 
								(crossp[0]*face_normal[0]+crossp[1]*face_normal[1]+crossp[2]*face_normal[2])/2;
#endif
						}
					}
				}

				double face_area = 0;



				for( int t = 0; t < nconvex; t++ )
				{
		
					int t2 = t+1;

					if( t2 >= nconvex ) t2 -= nconvex;

					double crossp[3];
					cross( ordered_pts+3*t, ordered_pts+3*t2, crossp );
					face_area += (crossp[0] * nrm[0] + crossp[1] * nrm[1] + crossp[2] * nrm[2])/2;
					

#ifdef DEBUG_VOLUME
					if( t < nconvex-1 )
					{
						drew = 1;
						fprintf(debug_triangles, "draw color %s\n", colors[q] );
						fprintf(debug_triangles, "draw triangle { %lf %lf %lf } { %lf %lf %lf } { %lf %lf %lf }\n",
							ordered_pts[0], ordered_pts[1], ordered_pts[2],
							ordered_pts[t*3+0], ordered_pts[t*3+1], ordered_pts[t*3+2],
							ordered_pts[t2*3+0], ordered_pts[t2*3+1], ordered_pts[t2*3+2]
							);
					}
#endif

#ifndef DISABLE_ALL
					if( dvol )
					{
						// deriv wrt crossp * nrm						

						int srcp1 = quad_src[q][ordered[t]];
						int srcp2 = quad_src[q][ordered[t2]];
						double prefac =  -(1.0/3.0) * (nrm[0] * ordered_pts[0] + nrm[1] * ordered_pts[1] + nrm[2] * ordered_pts[2]) /2;

						dcross( dvol, ordered_pts+t*3, ordered_pts+t2*3, points[srcp1].dp_dr, points[srcp2].dp_dr, prefac, nrm );						
					}
#endif										
				}
					
			
				face_area = fabs(face_area);

#ifndef DISABLE_ALL		
				the_volume += (1.0/3.0) * face_area * (nrm[0] * ordered_pts[0] + nrm[1] * ordered_pts[1] + nrm[2] * ordered_pts[2]);
				if( dvol )
				{
					int srcp = quad_src[q][ordered[0]];

					for( int c = 0; c < 3; c++ )
					for( int s = 0; s < 3; s++ )
					for( int sc = 0; sc < 3; sc++ )
					{
						// derivative wrt ordered_pts
						dvol[3*s+sc] += (1.0/3.0) * face_area * points[srcp].dp_dr[c*9+s*3+sc] * nrm[c] ;
						// derivative wrt normal
						dvol[3*s+sc] += (1.0/3.0) * face_area * ordered_pts[c] * dnrm[c*9+3*s+sc];
					}
	
				}
#endif

#ifdef DEBUG_VOLUME
				vol_q[q] += (1.0/3.0) * face_area * (nrm[0] * ordered_pts[0] + nrm[1] * ordered_pts[1] + nrm[2] * ordered_pts[2]);
#endif
			}
		}
	}
/*	else
	{
		// one has three, one has four.

	}	 */


#ifdef DEBUG_VOLUME	
	if( !drew )
	{
		printf("Didn't draw anything on call %d.\n", ncalls );
	}
#endif
	ncalls++;
	return the_volume;
} 

int build_out( 
		point_info *points,
		int *nuv,
		double *intersect_m, double *intersect_p, double *sign_move, int traverse_type, int cur_cell_in, double *r1, double *r2, double *r3, double *Ls)
{
	int the_cell[3] = { cur_cell_in/9, (cur_cell_in/3)%3, cur_cell_in%3 }; 
	double to_sort[6];
	int itype[6];
	int n_to_sort = 0;

	double fr_eps = -HARD_EPS;

	for( int c = 0; c < 3; c++ )
	{
		if( intersect_m[c] >= fr_eps && intersect_m[c] <= 1.0-fr_eps )
		{
			to_sort[n_to_sort] = intersect_m[c];
			itype[n_to_sort] = -(1+c);
			n_to_sort++;		
		}	
		if( intersect_p[c] >= fr_eps && intersect_p[c] <= 1.0-fr_eps )
		{
			to_sort[n_to_sort] = intersect_p[c];
			itype[n_to_sort] = (1+c);
			n_to_sort++;		
		}	
	}

	int done = 0;
	while(!done)
	{
		done = 1;
	
		for( int x = 0; x < n_to_sort-1; x++ )
		{
			if( to_sort[x+1] < to_sort[x] )
			{
				int iswap = itype[x];
				double val = to_sort[x];

				to_sort[x] = to_sort[x+1];
				itype[x] = itype[x+1];
				
				to_sort[x+1] = val;
				itype[x+1] = iswap;

				done = 0;
			}
		}	
	}

	for( int i = 0; i < n_to_sort; i++ )
	{
		points[*nuv].cells[0] = the_cell[0]*9+the_cell[1]*3+the_cell[2];
		int pm = 1;
		int axis = itype[i]-1;
		if( itype[i] < 0 )
		{
			axis = -itype[i]-1;
			pm = 0;
		}
//		if( itype[i] < 0 )
//			the_cell[axis] -= sign_move[axis];	
//		else
			the_cell[axis] += sign_move[axis];	
			
		points[*nuv].cells[1] = the_cell[0]*9+the_cell[1]*3+the_cell[2];
		points[*nuv].ncells = 2;
		double t=to_sort[i];
		int ip1 = 0;
		int ip2 = 0;	
		double *p1,*p2;
		if( traverse_type == 0 )
		{
			p1 = r1;
			p2 = r2;	
			ip1 = 0;
			ip2 = 1;
			points[*nuv].uv[0] = to_sort[i];	
			points[*nuv].uv[1] = 0;

//			uv_points[2*(*nuv)+0] = to_sort[i];
//			uv_points[2*(*nuv)+1] = 0;
		}
		else if ( traverse_type == 1 )
		{
			p1 = r3;
			p2 = r1;	
			ip1 = 2;
			ip2 = 0;
			points[*nuv].uv[0] = 0;	
			points[*nuv].uv[1] = 1-to_sort[i];
			
//			uv_points[2*(*nuv)+0] = 0;
//			uv_points[2*(*nuv)+1] = 1-to_sort[i];
		}
		else
		{
			p1 = r2;
			p2 = r3;	
			ip1 = 1;
			ip2 = 2;
			points[*nuv].uv[0] = 1-to_sort[i];	
			points[*nuv].uv[1] = to_sort[i];
//			uv_points[2*(*nuv)+0] = 1-to_sort[i];
//			uv_points[2*(*nuv)+1] = to_sort[i];
		}
		
		double d_t_d_p1[3] = { 
			0,0,0			
		};
		double d_t_d_p2[3] = { 
			0,0,0			
		};

		if( pm == 0 )
		{
			d_t_d_p1[axis] = -p1[axis] / (p2[axis]-p1[axis]) / (p2[axis]-p1[axis]) - 1.0 / (p2[axis]-p1[axis]);
			d_t_d_p2[axis] =  p1[axis] / (p2[axis]-p1[axis]) / (p2[axis]-p1[axis]);
		}
		else
		{
			d_t_d_p1[axis] =  (Ls[axis]-p1[axis]) / (p2[axis]-p1[axis]) / (p2[axis]-p1[axis]) - 1.0 / (p2[axis]-p1[axis]);
			d_t_d_p2[axis] = -(Ls[axis]-p1[axis]) / (p2[axis]-p1[axis]) / (p2[axis]-p1[axis]);
		}

		memset(points[*nuv].dp_dr,0,sizeof(double)*27);
		for( int dc = 0; dc < 3; dc++ ) // cartesian component of target vector		
		for( int c  = 0; c < 3; c++ ) // cartesian component of der vector
		{
			if( c == dc )
			{	
				points[*nuv].dp_dr[dc*9+ip1*3+c] = 1-t;
				points[*nuv].dp_dr[dc*9+ip2*3+c] = t;
			}
			// d t d r1x
		
			points[*nuv].dp_dr[dc*9+ip1*3+c] += d_t_d_p1[c] * (p2[dc]-p1[dc]);
			points[*nuv].dp_dr[dc*9+ip2*3+c] += d_t_d_p2[c] * (p2[dc]-p1[dc]);
		}
		points[*nuv].intersect_type = itype[i];	
		//intersect_type[*nuv] = itype[i];	
		(*nuv)++;
	}

	return the_cell[0]*9+the_cell[1]*3+the_cell[2];
}


void drawVMDArrow( FILE *theFile, double *from, double *to )
{
	double dr[3] = { to[0]-from[0], to[1]-from[1], to[2]-from[2] };
	double len = normalize(dr);

	double min_radius = 0.25;



	double interm[3] = { from[0] +0.8*len*dr[0],
			     from[1] +0.8*len*dr[1],
			     from[2] +0.8*len*dr[2] };
	fprintf( debug_triangles, "draw cylinder { %lf %lf %lf} {%lf %lf %lf} radius %lf\n",
		from[0], from[1], from[2],
		interm[0], interm[1], interm[2], (len * 0.1 > min_radius ? len*0.1 : min_radius) );
	
	fprintf( debug_triangles, "draw cone { %lf %lf %lf} {%lf %lf %lf} radius %lf\n",
		interm[0], interm[1], interm[2], 
		to[0], to[1], to[2], (len * 0.15 > min_radius ? len * 0.15 : min_radius) );

	//0.000000 8.654184 38.778522} {0.000000 8.653986 38.777034
	double pt_select_1[3] = { 0, 8.654184, 38.778522};
	double pt_select_2[3] = { 0, 8.653986, 38.777034};

	double dr1[3] = { 
		pt_select_1[0] - from[0],
		pt_select_1[1] - from[1],
		pt_select_1[2] - from[2] };
	
	double dr2[3] = { 
		pt_select_2[0] - interm[0],
		pt_select_2[1] - interm[1],
		pt_select_2[2] - interm[2] };

	double l1 =sqrt(dr1[0]*dr1[0]+dr1[1]*dr1[1]+dr1[2]*dr1[2]);
	double l2 =sqrt(dr2[0]*dr2[0]+dr2[1]*dr2[1]+dr2[2]*dr2[2]);

	if( l1 < 1e-3 && l2 < 1e-3 )
	{
		printf("Debug.\n");
	}
}	

/*					
void fillcp( double *der, double *r )
{ // fast is the vector component that we are differentiating
  // slow is the vector component that we are differentiating wrt
	der[0] = 0;
	der[1] = -r[2];
	der[2] =  r[1];
	der[3] =  r[2];
	der[4] =  0;
	der[5] = -r[0];
	der[6] = -r[1];
	der[7] =  r[0];
	der[8] = 0;
}

void fillcpder( double *der, double *dr, 
		int wrt,
		double sign
		)
{
	// what we are differentiating is slow (n_x, n_y, n_z)
	// what we are differentiating wrt is fast.
	// dnrmx, dr_x

	der[3*wrt+0] += sign * dr[0]; 	// x_comp
	der[3*wrt+1] += sign * dr[3]; 	
	der[3*wrt+2] += sign * dr[6]; 	
	der[9+3*wrt+0] += sign * dr[1+0];   // y_comp	
	der[9+3*wrt+1] += sign * dr[1+3]; 	
	der[9+3*wrt+2] += sign * dr[1+6]; 	
	der[18+3*wrt+0] += sign * dr[2+0];  // z_comp	
	der[18+3*wrt+1] += sign * dr[2+3]; 	
	der[18+3*wrt+2] += sign * dr[2+6]; 	
}
*/

/* derivative of the area of a triangle defined by r1,r2,r3, dnrm is the nrm deriv, fac is leading fac */
				
void darea( double *dvol, double *r1, double *r2, double *r3, double fac )
{
	double dnrm[27];
	memset(dnrm,0,sizeof(double)*27);
	double dr1[3] = {r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2]};	
	double dr2[3] = {r3[0]-r1[0],r3[1]-r1[1],r3[2]-r1[2]};	
	double unrm[3];
	cross( dr1, dr2, unrm );
	double d_cp1_dr[9], d_cp2_dr[9], d_cp3_dr[9];

	fillcp( d_cp1_dr, r1 );
	fillcp( d_cp2_dr, r2 );
	fillcp( d_cp3_dr, r3 );

	// (r2-r1) x ( r3-r1) 
	//  r2 x r3 - r1 x r3 - r2 x r1 

	// first, fill the derivative wrt the numerator.
	// r2 x r3 
	fillcpder( dnrm, d_cp3_dr, 1, +1 );
	fillcpder( dnrm, d_cp2_dr, 2, -1 ); //right side cross product gets neg

	// - r1 x r3
	fillcpder( dnrm, d_cp3_dr, 0, -1 ); 
	fillcpder( dnrm, d_cp1_dr, 2, +1 ); 
	
	// - r2 x r1
	fillcpder( dnrm, d_cp1_dr, 1, -1 );
	fillcpder( dnrm, d_cp2_dr, 0, +1 );

	// that's the numerator (the derivative of unrm)
	double d_nrm_len[9];
	memset( d_nrm_len, 0, sizeof(double) * 9 );
	
	double l = unrm[0]*unrm[0]+unrm[1]*unrm[1]+unrm[2]*unrm[2];

	// d_nrm_len is the derivative of unrm . unrm
	for( int p = 0; p < 3; p++ )
	for( int c = 0; c < 3; c++ )
	{
		double val = 2 * unrm[0] * dnrm[0 +3*p+c] +
			           2 * unrm[1] * dnrm[9 +3*p+c] +	 
			           2 * unrm[2] * dnrm[18+3*p+c];

		dvol[p*3+c] += 0.5*fac * val/(2*sqrt(l));
	}
}

void dcross( double *dvol, double *r1, double *r2, double *der1, double *der2, double factor, double *to_dot )
{
	// der1 is the derivative of the vector wrt the three source vectors.

	// factor * [d cross(r1,r2), src_{x,y,z}] . to_dot
	 
	double dcross_dr1[9];
	double dcross_dr2[9];

	// fast is what we're differentiating.
	fillcp( dcross_dr1, r2 );
	fillcp( dcross_dr2, r1 );

	for( int s = 0; s < 3; s++ )
	for( int sc = 0; sc < 3; sc++ )
	for( int c1 = 0; c1 < 3; c1++ ) // the component of the cross-product 
	for( int c2 = 0; c2 < 3; c2++ ) // the component of r1
	{
		// dcross is the derivative of the cross product with respect to r1
		// quel disastre. =^(
		dvol[3*s+sc] += factor * to_dot[c1] * dcross_dr1[c1+3*c2] * der1[c2*9+3*s+sc];
		dvol[3*s+sc] -= factor * to_dot[c1] * dcross_dr2[c1+3*c2] * der2[c2*9+3*s+sc];
	}
}
