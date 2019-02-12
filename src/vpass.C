#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mutil.h"
#include "vpass.h"

#define DEBUG_1
#define DEBUG_2
#define ALL_EDGES

static int max_edges = 30;
double compute_F( int ntarget, int ntri, int *tri, double *gamma_i, double *area_i, int *states, double PBC_vec[3][3], double *load_cen, double *aux_i );
void label_lattices( int *bonds, int nbonds, int KMED, int *labels );


void vpass( double *verts, int *tris, int *eft_tris, int ntri, int *edges, int nedges, int edge_dim, int ntarget, double PBC_vec[3][3] ) 
//		double *out_verts, int **out_tris, int **out_eft_tris, int *out_ntri, int **out_edges, int *out_nedges, double PBC_vec[3][3] )
{
	int *states = (int *)malloc( sizeof(int) * ntri );

	vpass_inner( states, verts, tris, eft_tris, ntri, edges, nedges, edge_dim, ntarget, PBC_vec  );
	{
		double *sgamma = (double *)malloc( sizeof(double) * ntarget * 3 );
		double *srho = (double *)malloc( sizeof(double) * ntarget );
		memset( sgamma, 0, sizeof(double) * 3 * ntarget );
		memset( srho, 0, sizeof(double) * ntarget );

		double *area_i = (double *)malloc( sizeof(double) * ntri );
		double *gamma_i = (double *)malloc( sizeof(double) * 3 * ntri );
		double putv[3];

		for( int t = 0; t < ntri; t++ )
		{
			int p1 = tris[3*t+0];
			int p2 = tris[3*t+1];
			int p3 = tris[3*t+2];
	
			double r1[3] = { verts[3*p1+0], verts[3*p1+1], verts[3*p1+2] };
			double dr2[3] = { verts[3*p2+0]-r1[0], verts[3*p2+1]-r1[1], verts[3*p2+2]-r1[2] };
			double dr3[3] = { verts[3*p3+0]-r1[0], verts[3*p3+1]-r1[1], verts[3*p3+2]-r1[2] };
	
			MinImage3D( dr2, PBC_vec, putv );
			MinImage3D( dr3, PBC_vec, putv );
	
			gamma_i[3*t+0] = r1[0]+(dr2[0]+dr3[0])/3.0;
			gamma_i[3*t+1] = r1[1]+(dr2[1]+dr3[1])/3.0;
			gamma_i[3*t+2] = r1[2]+(dr2[2]+dr3[2])/3.0;
	
			double o[3] = {0,0,0};
	
			area_i[t] = triangle_area( o, dr2, dr3 );
		}

		for( int t = 0; t < ntri; t++ )
		{
			int s = states[t];
			sgamma[3*s+0] += area_i[t] * gamma_i[3*t+0];
			sgamma[3*s+1] += area_i[t] * gamma_i[3*t+1];
			sgamma[3*s+2] += area_i[t] * gamma_i[3*t+2];
			srho[s] += area_i[t];
		}
		free(area_i);
		free(gamma_i);

		for( int s = 0; s < ntarget; s++ )
		{
			double put[3];
			double r1[3] = { sgamma[3*s+0] / srho[s], sgamma[3*s+1] / srho[s], sgamma[3*s+2] / srho[s] };
			MinImage3D(r1,PBC_vec,put);
			sgamma[3*s+0] = r1[0] * srho[s];
			sgamma[3*s+1] = r1[1] * srho[s];
			sgamma[3*s+2] = r1[2] * srho[s];
		}
		FILE *vFile = fopen("v.lattice", "w" );
		fprintf(vFile, "3D Lattice\n" );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[0][0], PBC_vec[0][1], PBC_vec[0][2] );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[1][0], PBC_vec[1][1], PBC_vec[1][2] );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[2][0], PBC_vec[2][1], PBC_vec[2][2] );
	
		int *lbonds = (int *)malloc( sizeof(int) * max_edges * ntarget );
		int *nlbonds = (int *)malloc( sizeof(int) * ntarget );
		memset( nlbonds, 0, sizeof(int) * ntarget );
	
		int *all_bonds = (int*)malloc( sizeof(int) * ntarget * max_edges);
		int xbonds = 0;
	
		for( int e = 0; e < nedges; e++ )
		{
			int t1 = edges[(edge_dim)*e+2];
			int t2 = edges[(edge_dim)*e+3];

			int s1 = states[t1];
			int s2 = states[t2];

			if( s1 == s2 ) continue; 

			int gotit =0;
			for( int x = 0; x < nlbonds[s1]; x++ )
			{
				if( lbonds[max_edges*s1+x] == s2 )
					gotit =1;
			}
			if( !gotit )
			{
				lbonds[s1*max_edges+nlbonds[s1]] = s2;
				nlbonds[s1]++;
			}
			gotit =0;
			for( int x = 0; x < nlbonds[s2]; x++ )
			{
				if( lbonds[max_edges*s2+x] == s1 )
					gotit =1;
			}
			if( !gotit )
			{
				lbonds[s2*max_edges+nlbonds[s2]] = s1;
				nlbonds[s2]++;
			}
		}
		
		for( int x = 0; x < ntarget; x++ )
		{
			fprintf(vFile, "%d %lf %lf %lf %d", x, sgamma[3*x+0] / srho[x], sgamma[3*x+1] / srho[x], sgamma[3*x+2] / srho[x], nlbonds[x] );
	
			for( int b = 0; b < nlbonds[x]; b++ )
				fprintf( vFile, " %d", lbonds[x*max_edges+b] );
			fprintf(vFile, "\n");
		} 
		fclose(vFile);
		free(sgamma);
		free(srho);
	}
/*
	for( int x = 0; x < ntri; x++ )
		states[x] = -1;

	double *area_i = (double *)malloc( sizeof(double) * ntri );
	double *gamma_i = (double *)malloc( sizeof(double) * 3 * ntri );
	double putv[3];

	for( int t = 0; t < ntri; t++ )
	{
		int p1 = tris[3*t+0];
		int p2 = tris[3*t+1];
		int p3 = tris[3*t+2];

		double r1[3] = { verts[3*p1+0], verts[3*p1+1], verts[3*p1+2] };
		double dr2[3] = { verts[3*p2+0]-r1[0], verts[3*p2+1]-r1[1], verts[3*p2+2]-r1[2] };
		double dr3[3] = { verts[3*p3+0]-r1[0], verts[3*p3+1]-r1[1], verts[3*p3+2]-r1[2] };

		MinImage3D( dr2, PBC_vec, putv );
		MinImage3D( dr3, PBC_vec, putv );

		gamma_i[3*t+0] = r1[0]+(dr2[0]+dr3[0])/3.0;
		gamma_i[3*t+1] = r1[1]+(dr2[1]+dr3[1])/3.0;
		gamma_i[3*t+2] = r1[2]+(dr2[2]+dr3[2])/3.0;

		double o[3] = {0,0,0};

		area_i[t] = triangle_area( o, dr2, dr3 );
	}
	
	double *sgamma = (double *)malloc( sizeof(double) * ntarget*3 );
	double *srho   = (double *)malloc( sizeof(double) * ntarget );
	int t = 0;
	
	while( t != ntarget )
	{
		int the_tri = rand() % ntri;
	
		if( states[the_tri] == -1 )	
		{
			sgamma[3*t+0] = gamma_i[3*the_tri+0] * area_i[the_tri];  
			sgamma[3*t+1] = gamma_i[3*the_tri+1] * area_i[the_tri];  
			sgamma[3*t+2] = gamma_i[3*the_tri+2] * area_i[the_tri];
			srho[t] = area_i[the_tri]; 
			states[the_tri] = t;

			t++;
		}
	}

	// assign the other states.
	
	int assignment_done = 0;
	int bad_pass = 0;
	while( !assignment_done )
	{
		assignment_done = 1;
		bad_pass = 1;

		for( int t = 0; t < ntri; t++ )
		{
			if( states[t] >= 0 ) continue;
	
			int best_target = 0;
			double best_dF = 1e10;
	
			for( int s = 0; s < ntarget; s++ )
			{
				double val1 = (sgamma[3*s+0]*sgamma[3*s+0] + sgamma[3*s+1]*sgamma[3*s+1] +sgamma[3*s+2]*sgamma[3*s+2])/srho[s];
				double dr[3] = { gamma_i[3*t+0] - sgamma[3*s+0]/srho[s], gamma_i[3*t+1] - sgamma[3*s+1]/srho[s], gamma_i[3*t+2] - sgamma[3*s+2]/srho[s] };
				MinImage3D( dr, PBC_vec, putv );
	
				dr[0] += sgamma[3*s+0]/srho[s];
				dr[1] += sgamma[3*s+1]/srho[s];
				dr[2] += sgamma[3*s+2]/srho[s];
	
				// with
				double val2 = 
					((sgamma[3*s+0]+area_i[t]*dr[0])*(sgamma[3*s+0]+area_i[t]*dr[0]) + 
					(sgamma[3*s+1]+area_i[t]*dr[1])*(sgamma[3*s+1]+area_i[t]*dr[1]) + 
					(sgamma[3*s+2]+area_i[t]*dr[2])*(sgamma[3*s+2]+area_i[t]*dr[2]))/(srho[s]+area_i[t]); 
				
				double dF = area_i[t] * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]) - val2 + val1;
	
				if( dF < best_dF )
				{
					best_target = s;
					best_dF = dF;
				}
			}
			//if( border_best || bad_pass )
			{
				int s = best_target;
				
				double dr[3] = { gamma_i[3*t+0] - sgamma[3*s+0]/srho[s], gamma_i[3*t+1] - sgamma[3*s+1]/srho[s], gamma_i[3*t+2] - sgamma[3*s+2]/srho[s] };
				MinImage3D( dr, PBC_vec, putv );
	
				dr[0] += sgamma[3*s+0]/srho[s];
				dr[1] += sgamma[3*s+1]/srho[s];
				dr[2] += sgamma[3*s+2]/srho[s];
	
				states[t] = s;
				sgamma[3*s+0] += area_i[t] * dr[0]; 
				sgamma[3*s+1] += area_i[t] * dr[1]; 
				sgamma[3*s+2] += area_i[t] * dr[2]; 
				srho[s] += area_i[t];
				bad_pass = 0;
			}
//			else
//				assignment_done = 0;
		}
	}


	int *border_edges = (int *)malloc( sizeof(int) * nedges );
	int *is_border = (int *)malloc( sizeof(int) * nedges );
	memset( is_border, 0, sizeof(int) * nedges );
	int nborder = 0;

	for( int e = 0; e < nedges; e++ )
	{
		int t1 = edges[edge_dim*e+2];
		int t2 = edges[edge_dim*e+3];

		if( t1 < 0 || t2 <0 ) continue;

		if( states[t1] != states[t2] )
		{
			is_border[e] = 1;
			border_edges[nborder] = e;
			nborder++;
		}
	}


	for( int pass = 0; pass < 2; pass ++ )
	{
		int v_done = 0;
		while( !v_done )
		{
			printf("Border iteration.\n");
			v_done = 1;
	
#ifdef ALL_EDGES
			for( int e = 0; e < nedges; e++ )
			{
#else
			for( int b = 0; b < nborder; b++ )
			{
				int e = border_edges[b];
#endif		
				int t1 = edges[edge_dim*e+2];
				int t2 = edges[edge_dim*e+3];
	

				double dF12 = 0; // one to two
				double dF21 = 0; // two to one
	
				int s1 = states[t1];
				int s2 = states[t2];
	
				if( s1 == s2 ) continue;
	
				double dr11[3] = { 
					gamma_i[3*t1+0] - sgamma[3*s1+0]/srho[s1],
					gamma_i[3*t1+1] - sgamma[3*s1+1]/srho[s1],
					gamma_i[3*t1+2] - sgamma[3*s1+2]/srho[s1] };
				double dr12[3] = { 
					gamma_i[3*t1+0] - sgamma[3*s2+0]/srho[s2],
					gamma_i[3*t1+1] - sgamma[3*s2+1]/srho[s2],
					gamma_i[3*t1+2] - sgamma[3*s2+2]/srho[s2] };
				
				double dr21[3] = { 
					gamma_i[3*t2+0] - sgamma[3*s1+0]/srho[s1],
					gamma_i[3*t2+1] - sgamma[3*s1+1]/srho[s1],
					gamma_i[3*t2+2] - sgamma[3*s1+2]/srho[s1] };
				double dr22[3] = { 
					gamma_i[3*t2+0] - sgamma[3*s2+0]/srho[s2],
					gamma_i[3*t2+1] - sgamma[3*s2+1]/srho[s2],
					gamma_i[3*t2+2] - sgamma[3*s2+2]/srho[s2] };
						
				MinImage3D( dr11, PBC_vec, putv );
				MinImage3D( dr12, PBC_vec, putv );
				MinImage3D( dr21, PBC_vec, putv );
				MinImage3D( dr22, PBC_vec, putv );
		
				dr11[0] += sgamma[3*s1+0]/srho[s1];
				dr11[1] += sgamma[3*s1+1]/srho[s1];
				dr11[2] += sgamma[3*s1+2]/srho[s1];
				
				dr12[0] += sgamma[3*s2+0]/srho[s2];
				dr12[1] += sgamma[3*s2+1]/srho[s2];
				dr12[2] += sgamma[3*s2+2]/srho[s2];
				
				dr21[0] += sgamma[3*s1+0]/srho[s1];
				dr21[1] += sgamma[3*s1+1]/srho[s1];
				dr21[2] += sgamma[3*s1+2]/srho[s1];
				
				dr22[0] += sgamma[3*s2+0]/srho[s2];
				dr22[1] += sgamma[3*s2+1]/srho[s2];
				dr22[2] += sgamma[3*s2+2]/srho[s2];
		
				double state_1_with_2 = 0;
				double state_2_with_1 = 0;
	
				double F11 = area_i[t1] * (dr11[0]*dr11[0]+dr11[1]*dr11[1]+dr11[2]*dr11[2]);
				double F12 = area_i[t1] * (dr12[0]*dr12[0]+dr12[1]*dr12[1]+dr12[2]*dr12[2]);
				double F21 = area_i[t2] * (dr21[0]*dr21[0]+dr21[1]*dr21[1]+dr21[2]*dr21[2]);
				double F22 = area_i[t2] * (dr22[0]*dr22[0]+dr22[1]*dr22[1]+dr22[2]*dr22[2]);
	
				double v1_1[3] = { sgamma[3*s1+0], sgamma[3*s1+1], sgamma[3*s1+2] };
				double v1_0[3] = { 
							v1_1[0] - dr11[0] * area_i[t1],
							v1_1[1] - dr11[1] * area_i[t1],
							v1_1[2] - dr11[2] * area_i[t1] };
	
	
				double v1_12[3] = { 
							v1_1[0] + dr21[0] * area_i[t2],
							v1_1[1] + dr21[1] * area_i[t2],
							v1_1[2] + dr21[2] * area_i[t2] };
				
				double v2_2[3] = { sgamma[3*s2+0], sgamma[3*s2+1], sgamma[3*s2+2] };
				double v2_0[3] = { 
							v2_2[0] - dr22[0] * area_i[t2],
							v2_2[1] - dr22[1] * area_i[t2],
							v2_2[2] - dr22[2] * area_i[t2] };
				double v2_12[3] = { 
							v2_2[0] + dr12[0] * area_i[t1],
							v2_2[1] + dr12[1] * area_i[t1],
							v2_2[2] + dr12[2] * area_i[t1] };
	
				double t_F1_00 = (v1_0[0]*v1_0[0] + v1_0[1]*v1_0[1] + v1_0[2]*v1_0[2])/(srho[s1]-area_i[t1]+1e-15);
				double t_F1_X0 = (v1_1[0]*v1_1[0] + v1_1[1]*v1_1[1] + v1_1[2]*v1_1[2])/(srho[s1]);
				double t_F1_XX = (v1_12[0]*v1_12[0] + v1_12[1]*v1_12[1] + v1_12[2]*v1_12[2])/(srho[s1]+area_i[t2]);
				double t_F2_00 = (v2_0[0]*v2_0[0] + v2_0[1]*v2_0[1] + v2_0[2]*v2_0[2])/(srho[s2]-area_i[t2]+1e-15);
				double t_F2_0X = (v2_2[0]*v2_2[0] + v2_2[1]*v2_2[1] + v2_2[2]*v2_2[2])/(srho[s2]);
				double t_F2_XX = (v2_12[0]*v2_12[0] + v2_12[1]*v2_12[1] + v2_12[2]*v2_12[2])/(srho[s2]+area_i[t1]);
	
				double F0 = F11 + F22 - t_F1_X0 - t_F2_0X;
				double F1 = F11 + F21 - t_F1_XX - t_F2_00;
				double F2 = F12 + F22 - t_F1_00 - t_F2_XX;
			
				if( F0 <= F1 && F0 <= F2 )
				{
	//				printf("REJECTING F0 %le F1 %le F2 %le\n", F0, F1, F2 );	
				}
				else 
				{
					
					v_done = 0;
#ifdef ALL_EDGES
	
#else
					is_border[e] = 0;
					border_edges[b] = border_edges[nborder-1];
					nborder--;
#endif
					double dF_expected = 0;
					int *es;
	
#ifdef DEBUG_1				
					double FCUR = 0;
					{
						for( int t = 0; t < ntri; t++ )
						{
							int st = states[t];
							double dr[3] = { gamma_i[3*t+0] - sgamma[3*st+0]/srho[st],
									 gamma_i[3*t+1] - sgamma[3*st+1]/srho[st],
									 gamma_i[3*t+2] - sgamma[3*st+2]/srho[st] };
							double put[3];
							MinImage3D( dr, PBC_vec, put );		
							dr[0] += sgamma[3*st+0]/srho[st];
							dr[1] += sgamma[3*st+1]/srho[st];
							dr[2] += sgamma[3*st+2]/srho[st];
							FCUR += area_i[t] *(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);					
						}	
						for( int s = 0; s < ntarget; s++ )
							FCUR -= (sgamma[3*s+0] * sgamma[3*s+0] + sgamma[3*s+1] * sgamma[3*s+1] + sgamma[3*s+2] * sgamma[3*s+2])/srho[s];
					}
#endif
					if( F1 <= F0 && F1 <= F2 )
					{
						states[t2] = s1;
						es = eft_tris+3*t2;

	
						sgamma[3*s2+0] -= dr22[0] * area_i[t2];	
						sgamma[3*s2+1] -= dr22[1] * area_i[t2];	
						sgamma[3*s2+2] -= dr22[2] * area_i[t2];	
						                         
						sgamma[3*s1+0] += dr21[0] * area_i[t2];	
						sgamma[3*s1+1] += dr21[1] * area_i[t2];	
						sgamma[3*s1+2] += dr21[2] * area_i[t2];	
						
						srho[s2] -= area_i[t2];
						srho[s1] += area_i[t2];
						dF_expected = F1 - F0;
					}
					else
					{
						states[t1] = s2;	
						es = eft_tris+3*t1;	
						
						sgamma[3*s1+0] -= dr11[0] * area_i[t1];	
						sgamma[3*s1+1] -= dr11[1] * area_i[t1];	
						sgamma[3*s1+2] -= dr11[2] * area_i[t1];	
						                         
						sgamma[3*s2+0] += dr12[0] * area_i[t1];	
						sgamma[3*s2+1] += dr12[1] * area_i[t1];	
						sgamma[3*s2+2] += dr12[2] * area_i[t1];	
	
						srho[s1] -= area_i[t1];
						srho[s2] += area_i[t1];
						dF_expected = F2 - F0;
					}
	
#ifndef ALL_EDGES
					for( int ex = 0; ex < 3; ex++ )
					{
						int ne = es[ex];
						if( is_border[ne] ) continue;
	
						int t1 = edges[edge_dim*ne+2];
						int t2 = edges[edge_dim*ne+3];
			
						if( t1 < 0 || t2 <0 ) continue;
	
						if( states[t1] != states[t2] )
						{
							is_border[ne] = 1;
							border_edges[nborder] = ne;
							nborder++;
						}
					}
#endif
	
#ifdef DEBUG_1				
					double FNEW = 0;
					{
						for( int t = 0; t < ntri; t++ )
						{
							int st = states[t];
							double dr[3] = { gamma_i[3*t+0] - sgamma[3*st+0]/srho[st],
									 gamma_i[3*t+1] - sgamma[3*st+1]/srho[st],
									 gamma_i[3*t+2] - sgamma[3*st+2]/srho[st] };
							double put[3];
							MinImage3D( dr, PBC_vec, put );		
							dr[0] += sgamma[3*st+0]/srho[st];
							dr[1] += sgamma[3*st+1]/srho[st];
							dr[2] += sgamma[3*st+2]/srho[st];
							FNEW += area_i[t] *(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);					
						}	
						for( int s = 0; s < ntarget; s++ )
							FNEW -= (sgamma[3*s+0] * sgamma[3*s+0] + sgamma[3*s+1] * sgamma[3*s+1] + sgamma[3*s+2] * sgamma[3*s+2])/srho[s];
					}
#ifdef DEBUG_2
					double abs_F = compute_F( ntarget, ntri, tris, gamma_i, area_i, states, PBC_vec);
#else
					double abs_F = -1;
#endif
	
					if( fabs(dF_expected-(FNEW-FCUR)) > 1e-4 )
					{
						printf("PROBLEM_dF_expected: %lf got: %lf FCUR %lf FNEW %lf ABSF %lf F0: %le F1: %le F2: %le\n", dF_expected, FNEW-FCUR, FCUR, FNEW, abs_F, F0, F1, F2  ); 
					}
					else
						printf("dF_expected: %lf got: %lf FCUR %lf FNEW %lf ABSF %lf F0: %le F1: %le F2: %le\n", dF_expected, FNEW-FCUR, FCUR, FNEW, abs_F, F0, F1, F2  ); 
#endif
				}
			}
#ifdef ALL_EDGES
		}
#else
		}
#endif

		for( int t = 0; t < ntri; t++ )
		{
			int e_t[3];
			int b = 0;
	
			for( int xt = 0; xt < 3; xt++ )
			{
				int be = eft_tris[3*t+xt];

				if( edges[edge_dim*be+2] == t )
					e_t[xt] = edges[edge_dim*be+3];
				else
					e_t[xt] = edges[edge_dim*be+2];
	
				if( states[t] == states[e_t[xt]] )
					b = 1;
			}

			if( ! b )
			{
				printf("triangle %d has no common borders.\n", t);
			}

		}

	}
	{
		FILE *vFile = fopen("v.lattice", "w" );
		fprintf(vFile, "3D Lattice\n" );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[0][0], PBC_vec[0][1], PBC_vec[0][2] );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[1][0], PBC_vec[1][1], PBC_vec[1][2] );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[2][0], PBC_vec[2][1], PBC_vec[2][2] );
	
		int *lbonds = (int *)malloc( sizeof(int) * max_edges * ntarget );
		int *nlbonds = (int *)malloc( sizeof(int) * ntarget );
		memset( nlbonds, 0, sizeof(int) * ntarget );
	
		int *all_bonds = (int*)malloc( sizeof(int) * ntarget * max_edges);
		int xbonds = 0;
	
		for( int s = 0; s < ntarget; s++ )
		{
			double put[3];
			double r1[3] = { sgamma[3*s+0] / srho[s], sgamma[3*s+1] / srho[s], sgamma[3*s+2] / srho[s] };
			MinImage3D(r1,PBC_vec,put);
			sgamma[3*s+0] = r1[0] * srho[s];
			sgamma[3*s+1] = r1[1] * srho[s];
			sgamma[3*s+2] = r1[2] * srho[s];
		}

		for( int e = 0; e < nedges; e++ )
		{
			int t1 = edges[(edge_dim)*e+2];
			int t2 = edges[(edge_dim)*e+3];

			int s1 = states[t1];
			int s2 = states[t2];

			if( s1 == s2 ) continue; 

			int gotit =0;
			for( int x = 0; x < nlbonds[s1]; x++ )
			{
				if( lbonds[max_edges*s1+x] == s2 )
					gotit =1;
			}
			if( !gotit )
			{
				lbonds[s1*max_edges+nlbonds[s1]] = s2;
				nlbonds[s1]++;
			}
			gotit =0;
			for( int x = 0; x < nlbonds[s2]; x++ )
			{
				if( lbonds[max_edges*s2+x] == s1 )
					gotit =1;
			}
			if( !gotit )
			{
				lbonds[s2*max_edges+nlbonds[s2]] = s1;
				nlbonds[s2]++;
			}
		}
		

		FILE *vxyz = fopen("c.xyz","w");
		fprintf(vxyz, "%d\n", ntarget );
		fprintf(vxyz, "c.xyz\n");
		for( int x = 0; x < ntarget;x++)
			fprintf(vxyz, "C %lf %lf %lf\n", sgamma[3*x+0]/srho[x], sgamma[3*x+1]/srho[x], sgamma[3*x+2]/srho[x] );
		fclose(vxyz);
	
		FILE *vpsf = fopen("c.psf","w");
		writePSF( vpsf, ntarget, NULL, all_bonds, xbonds );
		fclose(vpsf);
	
		for( int x = 0; x < ntarget; x++ )
		{
			fprintf(vFile, "%d %lf %lf %lf %d", x, sgamma[3*x+0] / srho[x], sgamma[3*x+1] / srho[x], sgamma[3*x+2] / srho[x], nlbonds[x] );
	
			for( int b = 0; b < nlbonds[x]; b++ )
				fprintf( vFile, " %d", lbonds[x*max_edges+b] );
			fprintf(vFile, "\n");
		} 
		fclose(vFile);
	
	FILE *triFile = fopen("v.vmd","w");
	
	for( int t = 0; t < ntri; t++ )
	{
		int curp = tris[3*t+0];
		int curp1 = tris[3*t+1];
		int curp2 = tris[3*t+2];
	
		const char *colors[] = { "red", "green", "blue", "white", "orange", "yellow", "pink", "purple", "cyan", "lime", "mauve", "ochre", "iceblue", "yellow2", "yellow3", "green2", "green3", "cyan2", "cyan3" };
		int ncol = sizeof(colors)/sizeof(const char *);

		int st = states[t];
		fprintf(triFile, "draw color %s\n", colors[st%ncol] );
		int x1 = tris[3*t+1];
		int x2 = tris[3*t+2];
		
		double p1[3] = { verts[3*curp+0], verts[3*curp+1], verts[3*curp+2] };
		double p2[3] = { verts[3*x1+0], verts[3*x1+1], verts[3*x1+2] };
		double p3[3] = { verts[3*x2+0], verts[3*x2+1], verts[3*x2+2] };

		double dr1[3] = { p2[0]-p1[0],p2[1]-p1[1],p2[2]-p1[2] };
		double dr2[3] = { p3[0]-p1[0],p3[1]-p1[1],p3[2]-p1[2] };
		double put[3];
	
		MinImage3D( dr1, PBC_vec, put );
		MinImage3D( dr2, PBC_vec, put );
	
		p2[0] = p1[0] + dr1[0];
		p2[1] = p1[1] + dr1[1];
		p2[2] = p1[2] + dr1[2];
		
		p3[0] = p1[0] + dr2[0];
		p3[1] = p1[1] + dr2[1];
		p3[2] = p1[2] + dr2[2];


		int dx = 0, dy = 0;
//		for( int dx = -1; dx <= 1; dx++ )
//		for( int dy = -1; dy <= 1; dy++ )
		{
			fprintf(triFile, "draw triangle { %lf %lf %lf } { %lf %lf %lf } { %lf %lf %lf }\n",
				verts[3*curp+0] + dx * PBC_vec[0][0], verts[3*curp+1]+dy*PBC_vec[1][1], verts[3*curp+2],
				p2[0]+dx*PBC_vec[0][0],p2[1]+dy*PBC_vec[1][1],p2[2],
				p3[0]+dx*PBC_vec[0][0],p3[1]+dy*PBC_vec[1][1],p3[2] 
				 );
		}
	}
	}
	free(sgamma);
	free(srho);
	free(gamma_i);
	free(area_i);
*/
	free(states);
}	

void cpass( double *gamma_i, double *gamma_nrm, int nv, int *bonds, int nbonds_per_site, int ntarget, double PBC_vec[3][3])
{
	int *states = (int *)malloc( sizeof(int) * nv );

	for( int x = 0; x < nv; x++ )
		states[x] = -1;

	double *area_i = (double *)malloc( sizeof(double) * nv );
	double putv[3];

	for( int t = 0; t < nv; t++ )
		area_i[t] = 1;

	double *nrm_init = (double* )malloc( sizeof(double) * ntarget * 3 );	
	double *sgamma = (double *)malloc( sizeof(double) * ntarget*3 );
	double *srho   = (double *)malloc( sizeof(double) * ntarget );
	int t = 0;
	
	while( t != ntarget )
	{
		int the_tri = rand() % nv;
	
		if( states[the_tri] == -1 )	
		{
			nrm_init[3*t+0] = gamma_nrm[3*the_tri+0];
			nrm_init[3*t+1] = gamma_nrm[3*the_tri+1];
			nrm_init[3*t+2] = gamma_nrm[3*the_tri+2];
			sgamma[3*t+0] = gamma_i[3*the_tri+0] * area_i[the_tri];  
			sgamma[3*t+1] = gamma_i[3*the_tri+1] * area_i[the_tri];  
			sgamma[3*t+2] = gamma_i[3*the_tri+2] * area_i[the_tri];
			srho[t] = area_i[the_tri]; 
			states[the_tri] = t;

			t++;
		}
	}

	// assign the other states.
	
	int assignment_done = 0;
	int bad_pass = 0;
	while( !assignment_done )
	{
		assignment_done = 1;
		bad_pass = 1;

		for( int t = 0; t < nv; t++ )
		{
			if( states[t] >= 0 ) continue;
	
			int best_target = 0;
			double best_dF = 1e10;
	
			for( int s = 0; s < ntarget; s++ )
			{
				double val1 = (sgamma[3*s+0]*sgamma[3*s+0] + sgamma[3*s+1]*sgamma[3*s+1] +sgamma[3*s+2]*sgamma[3*s+2])/srho[s];
				double dr[3] = { gamma_i[3*t+0] - sgamma[3*s+0]/srho[s], gamma_i[3*t+1] - sgamma[3*s+1]/srho[s], gamma_i[3*t+2] - sgamma[3*s+2]/srho[s] };
				MinImage3D( dr, PBC_vec, putv );

				double dp =gamma_nrm[3*t+0] * nrm_init[3*s+0] +  
					gamma_nrm[3*t+1] * nrm_init[3*s+1] +
					gamma_nrm[3*t+2] * nrm_init[3*s+2];

				if( dp < -0.5 ) continue;
	
				dr[0] += sgamma[3*s+0]/srho[s];
				dr[1] += sgamma[3*s+1]/srho[s];
				dr[2] += sgamma[3*s+2]/srho[s];
	
				// with
				double val2 = 
					((sgamma[3*s+0]+area_i[t]*dr[0])*(sgamma[3*s+0]+area_i[t]*dr[0]) + 
					(sgamma[3*s+1]+area_i[t]*dr[1])*(sgamma[3*s+1]+area_i[t]*dr[1]) + 
					(sgamma[3*s+2]+area_i[t]*dr[2])*(sgamma[3*s+2]+area_i[t]*dr[2]))/(srho[s]+area_i[t]); 
				
				double dF = area_i[t] * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]) - val2 + val1;
	
				if( dF < best_dF )
				{
					best_target = s;
					best_dF = dF;
				}
			}

			int border_best = 0;

			// is it connected to this state.
			for( int xt = 0; xt < nbonds_per_site; xt ++)
			{
				int t2 = bonds[t*nbonds_per_site+xt];

				if( states[t2] == best_target || states[t] == best_target )
					border_best = 1;
			}
	
			if( border_best || bad_pass )
			{
				int s = best_target;
				
				double dr[3] = { gamma_i[3*t+0] - sgamma[3*s+0]/srho[s], gamma_i[3*t+1] - sgamma[3*s+1]/srho[s], gamma_i[3*t+2] - sgamma[3*s+2]/srho[s] };
				MinImage3D( dr, PBC_vec, putv );
	
				dr[0] += sgamma[3*s+0]/srho[s];
				dr[1] += sgamma[3*s+1]/srho[s];
				dr[2] += sgamma[3*s+2]/srho[s];
	
				states[t] = s;
				sgamma[3*s+0] += area_i[t] * dr[0]; 
				sgamma[3*s+1] += area_i[t] * dr[1]; 
				sgamma[3*s+2] += area_i[t] * dr[2]; 
				srho[s] += area_i[t];
				bad_pass = 0;
			}
			else
				assignment_done = 0;
		}
	}

	for( int t = 0; t < ntarget; t++ )
	{
//		printf("state %d n %lf\n", t, srho[t] );
	}

	int *border_edges = (int *)malloc( sizeof(int) * 2 * nv * nbonds_per_site );
	int *is_border = (int *)malloc( sizeof(int) * nv * nbonds_per_site );
	memset( is_border, 0, sizeof(int) * nv * nbonds_per_site);
	int nborder = 0;

	for( int v = 0; v < nv; v++ )
	{
		for( int bx = 0; bx < nbonds_per_site; bx++ )
		{
			int v2 = bonds[v*nbonds_per_site+bx];

			if( v2 < v )
				continue;
			
			if( states[v] != states[v2] )
			{
				is_border[v*nbonds_per_site+bx] = 1;

				border_edges[nborder*2+0] = v;
				border_edges[nborder*2+1] = bx;
				nborder++;
			}
		}
	}

	int v_done = 0;

	while( !v_done )
	{
		printf("Border iteration.\n");
		v_done = 1;
		printf("nborder: %d\n", nborder);
		for( int b = 0; b < nborder; b++ )
		{
			int t1 = border_edges[b*2+0];
			int t2 = bonds[t1*nbonds_per_site+border_edges[b*2+1]];

			double dF12 = 0; // one to two
			double dF21 = 0; // two to one

			int s1 = states[t1];
			int s2 = states[t2];
		
//			printf("t1: %d t2: %d s1: %d s2: %d\n", t1, t2, s1, s2 );
			if( s1 == s2 ) 
				continue;
	
	
			double dr11[3] = { 
				gamma_i[3*t1+0] - sgamma[3*s1+0]/srho[s1],
				gamma_i[3*t1+1] - sgamma[3*s1+1]/srho[s1],
				gamma_i[3*t1+2] - sgamma[3*s1+2]/srho[s1] };
			double dr12[3] = { 
				gamma_i[3*t1+0] - sgamma[3*s2+0]/srho[s2],
				gamma_i[3*t1+1] - sgamma[3*s2+1]/srho[s2],
				gamma_i[3*t1+2] - sgamma[3*s2+2]/srho[s2] };
			
			double dr21[3] = { 
				gamma_i[3*t2+0] - sgamma[3*s1+0]/srho[s1],
				gamma_i[3*t2+1] - sgamma[3*s1+1]/srho[s1],
				gamma_i[3*t2+2] - sgamma[3*s1+2]/srho[s1] };
			double dr22[3] = { 
				gamma_i[3*t2+0] - sgamma[3*s2+0]/srho[s2],
				gamma_i[3*t2+1] - sgamma[3*s2+1]/srho[s2],
				gamma_i[3*t2+2] - sgamma[3*s2+2]/srho[s2] };
					
			MinImage3D( dr11, PBC_vec, putv );
			MinImage3D( dr12, PBC_vec, putv );
			MinImage3D( dr21, PBC_vec, putv );
			MinImage3D( dr22, PBC_vec, putv );
	
			dr11[0] += sgamma[3*s1+0]/srho[s1];
			dr11[1] += sgamma[3*s1+1]/srho[s1];
			dr11[2] += sgamma[3*s1+2]/srho[s1];
			
			dr12[0] += sgamma[3*s2+0]/srho[s2];
			dr12[1] += sgamma[3*s2+1]/srho[s2];
			dr12[2] += sgamma[3*s2+2]/srho[s2];
			
			dr21[0] += sgamma[3*s1+0]/srho[s1];
			dr21[1] += sgamma[3*s1+1]/srho[s1];
			dr21[2] += sgamma[3*s1+2]/srho[s1];
			
			dr22[0] += sgamma[3*s2+0]/srho[s2];
			dr22[1] += sgamma[3*s2+1]/srho[s2];
			dr22[2] += sgamma[3*s2+2]/srho[s2];
	
			double state_1_with_2 = 0;
			double state_2_with_1 = 0;

			double F11 = area_i[t1] * (dr11[0]*dr11[0]+dr11[1]*dr11[1]+dr11[2]*dr11[2]);
			double F12 = area_i[t1] * (dr12[0]*dr12[0]+dr12[1]*dr12[1]+dr12[2]*dr12[2]);
			double F21 = area_i[t2] * (dr21[0]*dr21[0]+dr21[1]*dr21[1]+dr21[2]*dr21[2]);
			double F22 = area_i[t2] * (dr22[0]*dr22[0]+dr22[1]*dr22[1]+dr22[2]*dr22[2]);

			double v1_1[3] = { sgamma[3*s1+0], sgamma[3*s1+1], sgamma[3*s1+2] };
			double v1_0[3] = { 
						v1_1[0] - dr11[0] * area_i[t1],
						v1_1[1] - dr11[1] * area_i[t1],
						v1_1[2] - dr11[2] * area_i[t1] };


			double v1_12[3] = { 
						v1_1[0] + dr21[0] * area_i[t2],
						v1_1[1] + dr21[1] * area_i[t2],
						v1_1[2] + dr21[2] * area_i[t2] };
			
			double v2_2[3] = { sgamma[3*s2+0], sgamma[3*s2+1], sgamma[3*s2+2] };
			double v2_0[3] = { 
						v2_2[0] - dr22[0] * area_i[t2],
						v2_2[1] - dr22[1] * area_i[t2],
						v2_2[2] - dr22[2] * area_i[t2] };
			double v2_12[3] = { 
						v2_2[0] + dr12[0] * area_i[t1],
						v2_2[1] + dr12[1] * area_i[t1],
						v2_2[2] + dr12[2] * area_i[t1] };

			double t_F1_00 = (v1_0[0]*v1_0[0] + v1_0[1]*v1_0[1] + v1_0[2]*v1_0[2])/(srho[s1]-area_i[t1]+1e-15);
			double t_F1_X0 = (v1_1[0]*v1_1[0] + v1_1[1]*v1_1[1] + v1_1[2]*v1_1[2])/(srho[s1]);
			double t_F1_XX = (v1_12[0]*v1_12[0] + v1_12[1]*v1_12[1] + v1_12[2]*v1_12[2])/(srho[s1]+area_i[t2]);
			double t_F2_00 = (v2_0[0]*v2_0[0] + v2_0[1]*v2_0[1] + v2_0[2]*v2_0[2])/(srho[s2]-area_i[t2]+1e-15);
			double t_F2_0X = (v2_2[0]*v2_2[0] + v2_2[1]*v2_2[1] + v2_2[2]*v2_2[2])/(srho[s2]);
			double t_F2_XX = (v2_12[0]*v2_12[0] + v2_12[1]*v2_12[1] + v2_12[2]*v2_12[2])/(srho[s2]+area_i[t1]);

			double F0 = F11 + F22 - t_F1_X0 - t_F2_0X;
			double F1 = F11 + F21 - t_F1_XX - t_F2_00;
			double F2 = F12 + F22 - t_F1_00 - t_F2_XX;

			if( F1 <= F0 && F1 <= F2 && srho[s2] < 1.5 )
			{	// can't remove single site
			}
			else if( F2 <= F0 && F2 <= F1 && srho[s1] < 1.5 )
			{	// can't remove single site
			}
			else if( F0 <= F1 && F0 <= F2 )
			{
				
			}
			else 
			{
				double FCUR = 0;
				{
					for( int t = 0; t < nv; t++ )
					{
						int st = states[t];
						double dr[3] = { gamma_i[3*t+0] - sgamma[3*st+0]/srho[st],
								 gamma_i[3*t+1] - sgamma[3*st+1]/srho[st],
								 gamma_i[3*t+2] - sgamma[3*st+2]/srho[st] };
						double put[3];
						MinImage3D( dr, PBC_vec, put );		
						dr[0] += sgamma[3*st+0]/srho[st];
						dr[1] += sgamma[3*st+1]/srho[st];
						dr[2] += sgamma[3*st+2]/srho[st];
						FCUR += area_i[t] *(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);					
					}	
					for( int s = 0; s < ntarget; s++ )
						FCUR -= (sgamma[3*s+0] * sgamma[3*s+0] + sgamma[3*s+1] * sgamma[3*s+1] + sgamma[3*s+2] * sgamma[3*s+2])/srho[s];
				}
				v_done = 0;
				is_border[t1*nbonds_per_site+border_edges[b*2+1]] = 0;
				border_edges[2*b+0] = border_edges[2*(nborder-1)+0];
				border_edges[2*b+1] = border_edges[2*(nborder-1)+1];
				nborder--;

				int site;

				double dF_expected = 0;
				if( F1 <= F0 && F1 <= F2 )
				{
					states[t2] = s1;
					site = t2;
					dF_expected = F1 - F0;
						if( srho[s2] < 1.5 )
						{
							printf("we are removing a single site, that can't be right\n");
							exit(1);
						}
					
					sgamma[3*s2+0] -= dr22[0] * area_i[t2];	
					sgamma[3*s2+1] -= dr22[1] * area_i[t2];	
					sgamma[3*s2+2] -= dr22[2] * area_i[t2];	
					
					sgamma[3*s1+0] += dr21[0] * area_i[t2];	
					sgamma[3*s1+1] += dr21[1] * area_i[t2];	
					sgamma[3*s1+2] += dr21[2] * area_i[t2];	

					srho[s2] -= area_i[t2];
					srho[s1] += area_i[t2];
				}
				else
				{
					dF_expected = F2 - F0;
					states[t1] = s2;	
					site = t1;
						
						if( srho[s1] < 1.5 )
						{
							printf("we are removing a single site, that can't be right\n");
							exit(1);
						}
					
					sgamma[3*s1+0] -= dr11[0] * area_i[t1];	
					sgamma[3*s1+1] -= dr11[1] * area_i[t1];	
					sgamma[3*s1+2] -= dr11[2] * area_i[t1];	
					                         
					sgamma[3*s2+0] += dr12[0] * area_i[t1];	
					sgamma[3*s2+1] += dr12[1] * area_i[t1];	
					sgamma[3*s2+2] += dr12[2] * area_i[t1];	

					srho[s1] -= area_i[t1];
					srho[s2] += area_i[t1];
				}

				for( int x = 0; x < nbonds_per_site; x++ )
				{
					int p1 = site;
					int p2 = bonds[p1*nbonds_per_site+x];
				
					if( states[p1] != states[p2] && !is_border[p1*nbonds_per_site+x] )
					{
						if( p1 < p2 )
						{
							is_border[p1*nbonds_per_site+x] = 1;
							border_edges[nborder*2+0] = p1;
							border_edges[nborder*2+1] = x;
							nborder++;
						}
						else
						{
							for( int x2 = 0;x2 < nbonds_per_site; x2++ )
							{
								if( bonds[p2*nbonds_per_site+x2] == p1 && !is_border[p2*nbonds_per_site+x2] )
								{
									is_border[p2*nbonds_per_site+x2] = 1;
									border_edges[nborder*2+0] = p2;
									border_edges[nborder*2+1] = x2;
									nborder++;
								}	
							}
						}
					}
				}
				
				double FNEW = 0;
				{
					
					for( int t = 0; t < nv; t++ )
					{
						int st = states[t];
						double dr[3] = { gamma_i[3*t+0] - sgamma[3*st+0]/srho[st],
								 gamma_i[3*t+1] - sgamma[3*st+1]/srho[st],
								 gamma_i[3*t+2] - sgamma[3*st+2]/srho[st] };
						double put[3];
						MinImage3D( dr, PBC_vec, put );		
						dr[0] += sgamma[3*st+0]/srho[st];
						dr[1] += sgamma[3*st+1]/srho[st];
						dr[2] += sgamma[3*st+2]/srho[st];
						FNEW += area_i[t] *(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);					
					}	
					for( int s = 0; s < ntarget; s++ )
						FNEW -= (sgamma[3*s+0] * sgamma[3*s+0] + sgamma[3*s+1] * sgamma[3*s+1] + sgamma[3*s+2] * sgamma[3*s+2])/srho[s];
				}

			//	printf("dF_expected: %lf got: %lf FCUR %lf FNEW %lf\n", dF_expected, FNEW-FCUR, FCUR, FNEW  ); 
			}
		}
	}

	// remove tri-valent states in which all three neighbors are connected to each other.




	const char *colors[] = { "red", "green", "blue", "white", "orange", "yellow", "pink", "purple", "cyan", "lime", "mauve", "ochre", "iceblue" };
	int ncol = sizeof(colors)/sizeof(const char *);
	

	{
		FILE *vFile = fopen("c.lattice", "w" );
		fprintf(vFile, "3D Lattice\n" );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[0][0], PBC_vec[0][1], PBC_vec[0][2] );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[1][0], PBC_vec[1][1], PBC_vec[1][2] );
		fprintf(vFile, "%lf %lf %lf\n", PBC_vec[2][0], PBC_vec[2][1], PBC_vec[2][2] );
	
		int *lbonds = (int *)malloc( sizeof(int) * max_edges * ntarget );
		int *nlbonds = (int *)malloc( sizeof(int) * ntarget );
		memset( nlbonds, 0, sizeof(int) * ntarget );
	
		int *all_bonds = (int*)malloc( sizeof(int) * ntarget * max_edges);
		int xbonds = 0;
		
		int *conn_bonds = (int*)malloc( sizeof(int) * ntarget * max_edges);
		int cbonds = 0;
	
		for( int s = 0; s < ntarget; s++ )
		{
			double put[3];
			double r1[3] = { sgamma[3*s+0] / srho[s], sgamma[3*s+1] / srho[s], sgamma[3*s+2] / srho[s] };
			MinImage3D(r1,PBC_vec,put);
			sgamma[3*s+0] = r1[0] * srho[s];
			sgamma[3*s+1] = r1[1] * srho[s];
			sgamma[3*s+2] = r1[2] * srho[s];
		}

		for( int t1 = 0; t1 < nv; t1++ )
		{
			for( int ex = 0; ex < nbonds_per_site; ex++ )
			{
				int t2 = bonds[t1*nbonds_per_site+ex];
	
				if( t2 < t1 ) continue;
		
				if ( t1 < 0 || t2 < 0 ) continue;
		
				int s1 = states[t1];
				int s2 = states[t2];
		
				if( s1 == s2 ) continue;
				if( s1 < 0 || s2 < 0 ) continue;

				int gotit = 0;
	
				for( int bx = 0; bx < nlbonds[s1]; bx++ )
				{
					if( lbonds[s1*max_edges+bx] == s2 )
						gotit = 1;
				}
	
				if( nlbonds[s1] >= max_edges || nlbonds[s2] >= max_edges )
				{
					printf("EDGE ERROR.\n");
					exit(1);
				}
		
				if( !gotit )
				{
					lbonds[s1*max_edges+nlbonds[s1]] = s2;
					lbonds[s2*max_edges+nlbonds[s2]] = s1;
		
					nlbonds[s1]++;
					nlbonds[s2]++;
		
					double r1[3] = { sgamma[3*s1+0]/srho[s1], sgamma[3*s1+1]/srho[s1], sgamma[3*s1+2]/srho[s1]  };
					double r2[3] = { sgamma[3*s2+0]/srho[s2], sgamma[3*s2+1]/srho[s2], sgamma[3*s2+2]/srho[s2]  };
					double dr[3] = { r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2] };
					double tput[3];
					MinImage3D(dr,PBC_vec,tput);
				
					conn_bonds[2*cbonds+0] = s1;
					conn_bonds[2*cbonds+1] = s2;
					cbonds++;
		
					if( tput[0]*tput[0]+tput[1]*tput[1]+tput[2]*tput[2] < 1e-5 )
					{
						all_bonds[2*xbonds+0] = s1;
						all_bonds[2*xbonds+1] = s2;
						xbonds++;
					}
				}
					
			}
		}
		
		{
			FILE *bxyz = fopen ("bc.xyz","w");

			const char * atom_types[] = { "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Na", "Cl", "I" };
			int nat = sizeof(atom_types) / sizeof(const char *);
	
			fprintf(bxyz, "%d\n", nv );
			fprintf(bxyz, "bc.xyz\n");
			for( int x = 0; x < nv;x++)
			{	
				int s1 = states[x];

				fprintf(bxyz, "%s %lf %lf %lf\n", atom_types[s1%nat], gamma_i[3*x+0], gamma_i[3*x+1], gamma_i[3*x+2] );
			}
			fclose(bxyz);
	
			int *all_bonds = (int *)malloc( sizeof(int) * nv * max_edges );
			int xbonds = 0;
	
			for( int x = 0; x < nv; x++ )
			{
				for( int bx = 0; bx < nbonds_per_site; bx++ )
				{
					int x2 = bonds[x*nbonds_per_site+bx];

					if( x2 > x )
					{
						double r1[3] = { gamma_i[3*x+0], gamma_i[3*x+1], gamma_i[3*x+2] };
						double r2[3] = { gamma_i[3*x2+0], gamma_i[3*x2+1], gamma_i[3*x2+2]  };
						double dr[3] = { r1[0]-r2[0],r1[1]-r2[1],r1[2]-r2[2] };
						double tput[3];
						MinImage3D(dr,PBC_vec,tput);
			
						if( tput[0]*tput[0]+tput[1]*tput[1]+tput[2]*tput[2] < 1e-5 )
						{
							all_bonds[xbonds*2+0] = x;
							all_bonds[xbonds*2+1] = x2;
							xbonds++;
						}
					}
				}
			}

			FILE *bpsf = fopen("bc.psf","w");
			writePSF( bpsf, nv, NULL, all_bonds, xbonds );
			fclose(bpsf);
			
		}		

		int elim[ntarget];
		int map[ntarget];
		memset( elim, 0, sizeof(int) * ntarget );

		for( int x = 0; x < ntarget; x++ )
		{
			if( nlbonds[x] == 3 )
			{	
				int *ps = lbonds+x*max_edges;	

				int bonds[9];
				memset( bonds, 0, sizeof(int) * 9 );

				int ok = 0;

				for( int i = 0; i < 3; i++ )
				for( int j = 0; j < 3; j++ )
				{
					if( i == j ) continue;
						
					for( int b = 0; b < nlbonds[ps[i]]; b++ )
					{
						if( lbonds[ps[i]*max_edges+b] == ps[j] )
							bonds[i*3+j] = 1;
					}

					if( bonds[i*3+j] != 1 ) ok = 1;
				}

				if( !ok  )
					elim[x] = 1;
			}
		}

		int m =0;
		for( int x = 0; x < ntarget; x++ )
		{
			if( !elim[x] )
			{
				map[x] = m;
				m++;
			}
		}


		int labels[ntarget];
		label_lattices( conn_bonds, cbonds, ntarget, labels );
	
		FILE *vxyz = fopen("c.xyz","w");
		fprintf(vxyz, "%d\n", ntarget );
		fprintf(vxyz, "c.xyz\n");
		for( int x = 0; x < ntarget;x++)
			fprintf(vxyz, "C %lf %lf %lf\n", sgamma[3*x+0]/srho[x], sgamma[3*x+1]/srho[x], sgamma[3*x+2]/srho[x] );
		fclose(vxyz);
	
		FILE *vpsf = fopen("c.psf","w");
		writePSF( vpsf, ntarget, NULL, all_bonds, xbonds );
		fclose(vpsf);
	
		for( int x = 0; x < ntarget; x++ )
		{
			int nb = 0;

			if( elim[x] ) continue;

			for( int b = 0; b < nlbonds[x]; b++ )
			{
				if( !elim[lbonds[x*max_edges+b]] ) nb++;
			}

			fprintf(vFile, "%d %lf %lf %lf %d", map[x], sgamma[3*x+0] / srho[x], sgamma[3*x+1] / srho[x], sgamma[3*x+2] / srho[x], nb );
	
			for( int b = 0; b < nlbonds[x]; b++ )
			{
				if( !elim[lbonds[x*max_edges+b]] )
					fprintf( vFile, " %d", map[lbonds[x*max_edges+b]] );
			}
			fprintf(vFile, "\n");
		} 
		fclose(vFile);
	}
	free(sgamma);
	free(srho);
	free(area_i);
	free(states);
	
}

double compute_F( int ntarget, int ntri, int *tri, double *gamma_i, double *area_i, int *states, double PBC_vec[3][3], double *load_cens, double *aux_i, double aux_w )
{
	double *sgamma = (double *)malloc( sizeof(double) * ntarget * 3 );
	memset( sgamma, 0, sizeof(double) * 3 * ntarget );
	double *srho = (double *)malloc( sizeof(double) * ntarget  );
	memset( srho, 0, sizeof(double) *  ntarget );
	double *saux = (double *)malloc( sizeof(double) * ntarget  );
	memset( saux, 0, sizeof(double) *  ntarget );

	double putv[3];
	double F = 0;

	for( int t = 0; t < ntri; t++ )
	{
		int s = states[t];

		double dr[3] = { 
			gamma_i[3*t+0] - load_cens[3*s+0],
			gamma_i[3*t+1] - load_cens[3*s+1],
			gamma_i[3*t+2] - load_cens[3*s+2] };
				
		MinImage3D( dr, PBC_vec, putv );
	
		dr[0] += load_cens[3*s+0];
		dr[1] += load_cens[3*s+1];
		dr[2] += load_cens[3*s+2];

		F += area_i[t] * ( dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
		if( aux_i ) F += aux_w * area_i[t] * aux_i[t] * aux_i[t]; 				
		sgamma[3*s+0] += dr[0] * area_i[t];
		sgamma[3*s+1] += dr[1] * area_i[t];
		sgamma[3*s+2] += dr[2] * area_i[t];
							
		srho[s] += area_i[t];
	}
	
	for( int s = 0; s < ntarget; s++ )
	{
//		printf("load_cen %le %le %le sgamma %le %le %le\n", load_cens[3*s+0], load_cens[3*s+1], load_cens[3*s+2], sgamma[3*s+0]/srho[s], sgamma[3*s+1]/srho[s], sgamma[3*s+2]/srho[s] );
		F -= (sgamma[3*s+0]*sgamma[3*s+0] + sgamma[3*s+1]*sgamma[3*s+1]+sgamma[3*s+2]*sgamma[3*s+2])/srho[s];
		if( aux_i ) F -= aux_w * saux[s] * saux[s] / srho[s]; 				
	}

	free(sgamma);
	free(srho);
	free(saux);

	return F;
}

void label_lattices( int *bonds, int nbonds, int KMED, int *labels )
{ 
	for( int i = 0; i < KMED; i++ )
		labels[i] = i;
	
	int start_pt = 0;
	
	int done = 0;

	while( !done )
	{
		done = 1;

		int touched[KMED];
		memset( touched, 0, sizeof(int) * KMED );
	
		for( int b = 0; b < nbonds; b++ )
		{
			int i1 = bonds[2*b+0];	
			int i2 = bonds[2*b+1];

			if( touched[i1] || touched[i2] ) continue;

			if( labels[i1] == labels[i2] ) continue;

			if( labels[i1] < labels[i2] )
			{
				if( labels[i1] == 0 )
					printf("index %d or ", i2 );
				labels[i2] = labels[i1];
				touched[i2] = 1;
			}
			else
			{
				labels[i1] = labels[i2];	
				touched[i1] = 1;
			}
			done = 0;
		}	
		printf("\n");
	}	

	int npts[KMED];
	memset( npts, 0, sizeof(int) * KMED );

	for( int p = 0; p < KMED; p++ )
	{
		npts[labels[p]] += 1;
	}
	for( int p = 0; p < KMED; p++ )
	{
		if( npts[p] > 0 ) 
			printf("Connected group of size %d.\n", npts[p] );
	}
}

// the real vpass loop
void vpass_inner(
		int *states, // the output. 
		double *verts, int *tris, int *eft_tris, int ntri, int *edges, int nedges, int edge_dim, int ntarget, double PBC_vec[3][3] ) 
{
	int n_fix_regions = 0;
	int *fix_regions = NULL;
	double *aux_i = NULL;
	double aux_w = 0;

	int fixed[ntri];
	memset( fixed, 0, sizeof(int) * ntri );
	for( int x = 0; x < ntri; x++ )
		states[x] = -1;


	double *area_i = (double *)malloc( sizeof(double) * ntri );
	double *gamma_i = (double *)malloc( sizeof(double) * 3 * ntri );
	double putv[3];

	for( int t = 0; t < ntri; t++ )
	{
		int p1 = tris[3*t+0];
		int p2 = tris[3*t+1];
		int p3 = tris[3*t+2];

		double r1[3] = { verts[3*p1+0], verts[3*p1+1], verts[3*p1+2] };
		double dr2[3] = { verts[3*p2+0]-r1[0], verts[3*p2+1]-r1[1], verts[3*p2+2]-r1[2] };
		double dr3[3] = { verts[3*p3+0]-r1[0], verts[3*p3+1]-r1[1], verts[3*p3+2]-r1[2] };

		MinImage3D( dr2, PBC_vec, putv );
		MinImage3D( dr3, PBC_vec, putv );

		gamma_i[3*t+0] = r1[0]+(dr2[0]+dr3[0])/3.0;
		gamma_i[3*t+1] = r1[1]+(dr2[1]+dr3[1])/3.0;
		gamma_i[3*t+2] = r1[2]+(dr2[2]+dr3[2])/3.0;

		double o[3] = {0,0,0};

		area_i[t] = triangle_area( o, dr2, dr3 );
	}
	
	double *sgamma = (double *)malloc( sizeof(double) * ntarget*3 );
	memset( sgamma, 0, sizeof(double) * ntarget*3);
	double *srho   = (double *)malloc( sizeof(double) * ntarget );
	memset( srho, 0, sizeof(double) * ntarget);
	double *saux   = (double *)malloc( sizeof(double) * ntarget );
	memset( saux, 0, sizeof(double) * ntarget);

	// assign any fixed states right here.
	int assigned[ntarget];
	memset( assigned, 0, sizeof(int) * ntarget );
	
	if( n_fix_regions > 0 )
	{	
		for( int i = 0; i < n_fix_regions; i++ )
		{
			int s = fix_regions[2*i+1];
			int tri = fix_regions[2*i+0];

			if( s >= ntarget )
			{
				printf("Error: attempting to fix a region outside the bounds.\n");
				exit(1);
			}
			if( tri >= ntri )
			{
				printf("Error: attempting to fix a triangle outside the bounds.\n");
				exit(1);
			}

			states[tri] = s;
			assigned[s] = 1;
			sgamma[3*s+0] += gamma_i[3*tri+0] * area_i[tri];
			sgamma[3*s+1] += gamma_i[3*tri+1] * area_i[tri];
			sgamma[3*s+2] += gamma_i[3*tri+2] * area_i[tri];
			srho[s] += area_i[tri];
			if( aux_i ) saux[s] += area_i[tri] * aux_i[tri];

			fixed[tri] = 1;
		}
	}
	int t = 0;
	
	while( t != ntarget )
	{
		if( assigned[t] )
		{
			t++;
			continue;
		} 

		int the_tri = rand() % ntri;
	
		if( states[the_tri] == -1 )	
		{
			sgamma[3*t+0] = gamma_i[3*the_tri+0] * area_i[the_tri];  
			sgamma[3*t+1] = gamma_i[3*the_tri+1] * area_i[the_tri];  
			sgamma[3*t+2] = gamma_i[3*the_tri+2] * area_i[the_tri];
			srho[t] = area_i[the_tri]; 

			if( aux_i ) saux[t] = area_i[the_tri] * aux_i[the_tri];

			states[the_tri] = t;

			t++;
		}
	}

	// assign the other states.
	
	int assignment_done = 0;
	int bad_pass = 0;
	while( !assignment_done )
	{
		assignment_done = 1;
		bad_pass = 1;

		for( int t = 0; t < ntri; t++ )
		{
			if( states[t] >= 0 ) continue;
	
			int best_target = 0;
			double best_dF = 1e10;
	
			for( int s = 0; s < ntarget; s++ )
			{
				double val1 = (sgamma[3*s+0]*sgamma[3*s+0] + sgamma[3*s+1]*sgamma[3*s+1] +sgamma[3*s+2]*sgamma[3*s+2])/srho[s];
				double dr[3] = { gamma_i[3*t+0] - sgamma[3*s+0]/srho[s], gamma_i[3*t+1] - sgamma[3*s+1]/srho[s], gamma_i[3*t+2] - sgamma[3*s+2]/srho[s] };
				MinImage3D( dr, PBC_vec, putv );
	
				dr[0] += sgamma[3*s+0]/srho[s];
				dr[1] += sgamma[3*s+1]/srho[s];
				dr[2] += sgamma[3*s+2]/srho[s];
	
				// with
				double val2 = 
					((sgamma[3*s+0]+area_i[t]*dr[0])*(sgamma[3*s+0]+area_i[t]*dr[0]) + 
					(sgamma[3*s+1]+area_i[t]*dr[1])*(sgamma[3*s+1]+area_i[t]*dr[1]) + 
					(sgamma[3*s+2]+area_i[t]*dr[2])*(sgamma[3*s+2]+area_i[t]*dr[2]))/(srho[s]+area_i[t]); 
				
				double dF = area_i[t] * (dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]) - val2 + val1;
	
				if( dF < best_dF )
				{
					best_target = s;
					best_dF = dF;
				}
			}
/*
			int border_best = 0;
			int *es = eft_tris+3*t;

			for( int xt = 0; xt < 3; xt ++)
			{
				int e = es[xt];
				
				int t1 = edges[edge_dim*e+2];
				int t2 = edges[edge_dim*e+3];

				if( states[t1] == best_target || states[t2] == best_target )
					border_best = 1;
			}
	*/
			//if( border_best || bad_pass )
			{
				int s = best_target;
				
				double dr[3] = { gamma_i[3*t+0] - sgamma[3*s+0]/srho[s], gamma_i[3*t+1] - sgamma[3*s+1]/srho[s], gamma_i[3*t+2] - sgamma[3*s+2]/srho[s] };
				MinImage3D( dr, PBC_vec, putv );
	
				dr[0] += sgamma[3*s+0]/srho[s];
				dr[1] += sgamma[3*s+1]/srho[s];
				dr[2] += sgamma[3*s+2]/srho[s];
	
				states[t] = s;
				sgamma[3*s+0] += area_i[t] * dr[0]; 
				sgamma[3*s+1] += area_i[t] * dr[1]; 
				sgamma[3*s+2] += area_i[t] * dr[2];

				if( aux_i ) saux[s] += area_i[t] * aux_i[t]; 
				srho[s] += area_i[t];

				bad_pass = 0;
			}
//			else
//				assignment_done = 0;
		}
	}


	int *border_edges = (int *)malloc( sizeof(int) * nedges );
	int *is_border = (int *)malloc( sizeof(int) * nedges );
	memset( is_border, 0, sizeof(int) * nedges );
	int nborder = 0;

	for( int e = 0; e < nedges; e++ )
	{
		int t1 = edges[edge_dim*e+2];
		int t2 = edges[edge_dim*e+3];

		if( t1 < 0 || t2 <0 ) continue;

		if( states[t1] != states[t2] )
		{
			is_border[e] = 1;
			border_edges[nborder] = e;
			nborder++;
		}
	}

	double load_cens[3*ntarget];
	memcpy( load_cens, sgamma, 3 * ntarget * sizeof(double) );
	for( int s = 0; s < ntarget; s++ )
	{
		load_cens[3*s+0] /= srho[s];
		load_cens[3*s+1] /= srho[s];
		load_cens[3*s+2] /= srho[s];
	}

	aux_w = 0;

	double running_F = compute_F( ntarget, ntri, tris, gamma_i, area_i, states, PBC_vec, load_cens, aux_i, aux_w);


	for( int pass = 0; pass < 2; pass ++ )
	{
		int v_done = 0;
		int niter = 0;
		while( !v_done )
		{
			printf("Border iteration.\n");
			v_done = 1;
			niter++;
			if( niter > 50 )
				break;
#ifdef ALL_EDGES
			for( int e = 0; e < nedges; e++ )
			{
#else
			for( int b = 0; b < nborder; b++ )
			{
				int e = border_edges[b];
#endif		
	memcpy( load_cens, sgamma, 3 * ntarget*sizeof(double) );
	for( int s = 0; s < ntarget; s++ )
	{
		load_cens[3*s+0] /= srho[s];
		load_cens[3*s+1] /= srho[s];
		load_cens[3*s+2] /= srho[s];
	}
///				printf("CheckF running %le fn %le\n", running_F, compute_F( ntarget, ntri, tris, gamma_i, area_i, states, PBC_vec, load_cens, aux_i, aux_w ) );

				int t1 = edges[edge_dim*e+2];
				int t2 = edges[edge_dim*e+3];
	
				if( fixed[t1] || fixed[t2] ) continue;

				double dF12 = 0; // one to two
				double dF21 = 0; // two to one
	
				int s1 = states[t1];
				int s2 = states[t2];
	
				if( s1 == s2 ) continue;
	

				double dr11[3] = { 
					gamma_i[3*t1+0] - sgamma[3*s1+0]/srho[s1],
					gamma_i[3*t1+1] - sgamma[3*s1+1]/srho[s1],
					gamma_i[3*t1+2] - sgamma[3*s1+2]/srho[s1] };
				double dr12[3] = { 
					gamma_i[3*t1+0] - sgamma[3*s2+0]/srho[s2],
					gamma_i[3*t1+1] - sgamma[3*s2+1]/srho[s2],
					gamma_i[3*t1+2] - sgamma[3*s2+2]/srho[s2] };
				
				double dr21[3] = { 
					gamma_i[3*t2+0] - sgamma[3*s1+0]/srho[s1],
					gamma_i[3*t2+1] - sgamma[3*s1+1]/srho[s1],
					gamma_i[3*t2+2] - sgamma[3*s1+2]/srho[s1] };
				double dr22[3] = { 
					gamma_i[3*t2+0] - sgamma[3*s2+0]/srho[s2],
					gamma_i[3*t2+1] - sgamma[3*s2+1]/srho[s2],
					gamma_i[3*t2+2] - sgamma[3*s2+2]/srho[s2] };

				MinImage3D( dr11, PBC_vec, putv );
				MinImage3D( dr12, PBC_vec, putv );
				MinImage3D( dr21, PBC_vec, putv );
				MinImage3D( dr22, PBC_vec, putv );
		
				dr11[0] += sgamma[3*s1+0]/srho[s1];
				dr11[1] += sgamma[3*s1+1]/srho[s1];
				dr11[2] += sgamma[3*s1+2]/srho[s1];
				
				dr12[0] += sgamma[3*s2+0]/srho[s2];
				dr12[1] += sgamma[3*s2+1]/srho[s2];
				dr12[2] += sgamma[3*s2+2]/srho[s2];
				
				dr21[0] += sgamma[3*s1+0]/srho[s1];
				dr21[1] += sgamma[3*s1+1]/srho[s1];
				dr21[2] += sgamma[3*s1+2]/srho[s1];
				
				dr22[0] += sgamma[3*s2+0]/srho[s2];
				dr22[1] += sgamma[3*s2+1]/srho[s2];
				dr22[2] += sgamma[3*s2+2]/srho[s2];
		

				double state_1_with_2 = 0;
				double state_2_with_1 = 0;
	
				double F11 = area_i[t1] * (dr11[0]*dr11[0]+dr11[1]*dr11[1]+dr11[2]*dr11[2]);
				double F12 = area_i[t1] * (dr12[0]*dr12[0]+dr12[1]*dr12[1]+dr12[2]*dr12[2]);
				double F21 = area_i[t2] * (dr21[0]*dr21[0]+dr21[1]*dr21[1]+dr21[2]*dr21[2]);
				double F22 = area_i[t2] * (dr22[0]*dr22[0]+dr22[1]*dr22[1]+dr22[2]*dr22[2]);
	
				double v1_1[3] = { sgamma[3*s1+0], sgamma[3*s1+1], sgamma[3*s1+2] };
				double v1_0[3] = { 
							v1_1[0] - dr11[0] * area_i[t1],
							v1_1[1] - dr11[1] * area_i[t1],
							v1_1[2] - dr11[2] * area_i[t1] };
	
	
				double v1_12[3] = { 
							v1_1[0] + dr21[0] * area_i[t2],
							v1_1[1] + dr21[1] * area_i[t2],
							v1_1[2] + dr21[2] * area_i[t2] };
				
				double v2_2[3] = { sgamma[3*s2+0], sgamma[3*s2+1], sgamma[3*s2+2] };
				double v2_0[3] = { 
							v2_2[0] - dr22[0] * area_i[t2],
							v2_2[1] - dr22[1] * area_i[t2],
							v2_2[2] - dr22[2] * area_i[t2] };
				double v2_12[3] = { 
							v2_2[0] + dr12[0] * area_i[t1],
							v2_2[1] + dr12[1] * area_i[t1],
							v2_2[2] + dr12[2] * area_i[t1] };
	
				double t_F1_00 = (v1_0[0]*v1_0[0] + v1_0[1]*v1_0[1] + v1_0[2]*v1_0[2])/(srho[s1]-area_i[t1]+1e-15);
				double t_F1_X0 = (v1_1[0]*v1_1[0] + v1_1[1]*v1_1[1] + v1_1[2]*v1_1[2])/(srho[s1]);
				double t_F1_XX = (v1_12[0]*v1_12[0] + v1_12[1]*v1_12[1] + v1_12[2]*v1_12[2])/(srho[s1]+area_i[t2]);
				double t_F2_00 = (v2_0[0]*v2_0[0] + v2_0[1]*v2_0[1] + v2_0[2]*v2_0[2])/(srho[s2]-area_i[t2]+1e-15);
				double t_F2_0X = (v2_2[0]*v2_2[0] + v2_2[1]*v2_2[1] + v2_2[2]*v2_2[2])/(srho[s2]);
				double t_F2_XX = (v2_12[0]*v2_12[0] + v2_12[1]*v2_12[1] + v2_12[2]*v2_12[2])/(srho[s2]+area_i[t1]);
	
				double F0 = F11 + F22 - t_F1_X0 - t_F2_0X;
				double F1 = F11 + F21 - t_F1_XX - t_F2_00;
				double F2 = F12 + F22 - t_F1_00 - t_F2_XX;
		
				if( aux_i )
				{
					// moving the second into one.
					
					double A_10 = aux_w * (saux[s1]) * (saux[s1]) / srho[s1];
					double A_20 = aux_w * (saux[s2]) * (saux[s2]) / srho[s2];

					double A_11 = aux_w * (saux[s1] + area_i[t2] * aux_i[t2]) * (saux[s1] + area_i[t2] * aux_i[t2]) / (srho[s1]+area_i[t2]);
					double A_21 = aux_w * (saux[s2] - area_i[t2] * aux_i[t2]) * (saux[s2] - area_i[t2] * aux_i[t2]) / (srho[s2]-area_i[t2]);
					
					double A_12 = aux_w * (saux[s1] - area_i[t1] * aux_i[t1]) * (saux[s1] - area_i[t1] * aux_i[t1]) / (srho[s1]-area_i[t1]);
					double A_22 = aux_w * (saux[s2] + area_i[t1] * aux_i[t1]) * (saux[s2] + area_i[t1] * aux_i[t1]) / (srho[s2]+area_i[t1]);
					
					F0 += -A_10 - A_20;
					F1 += -A_11 - A_21;
					F2 += -A_12 - A_22;
				}
				if( F0 <= F1 && F0 <= F2 )
				{
				}
				else 
				{
#ifdef DEBUG_2
					double abs_F_p = compute_F( ntarget, ntri, tris, gamma_i, area_i, states, PBC_vec, load_cens, aux_i, aux_w);
	memcpy( load_cens, sgamma, 3 * ntarget * sizeof(double) );
	for( int s = 0; s < ntarget; s++ )
	{
		load_cens[3*s+0] /= srho[s];
		load_cens[3*s+1] /= srho[s];
		load_cens[3*s+2] /= srho[s];
	}
#else
					double abs_F_p = -1;
#endif
					
					v_done = 0;
#ifdef ALL_EDGES
	
#else
					is_border[e] = 0;
					border_edges[b] = border_edges[nborder-1];
					nborder--;
#endif
					double dF_expected = 0;
					int *es;
	
#ifdef DEBUG_1				
					double FCUR = 0;
					{
						double check_sgamma[3*ntarget];
						memset( check_sgamma, 0, sizeof(double) * 3 * ntarget );
						for( int t = 0; t < ntri; t++ )
						{
							int st = states[t];
							double dr[3] = { gamma_i[3*t+0] - sgamma[3*st+0]/srho[st],
									 gamma_i[3*t+1] - sgamma[3*st+1]/srho[st],
									 gamma_i[3*t+2] - sgamma[3*st+2]/srho[st] };
							double put[3];
							MinImage3D( dr, PBC_vec, put );		

							dr[0] += sgamma[3*st+0]/srho[st];
							dr[1] += sgamma[3*st+1]/srho[st];
							dr[2] += sgamma[3*st+2]/srho[st];
							FCUR += area_i[t] *(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
							if( aux_i ) FCUR += aux_w * aux_i[t] * aux_i[t] * area_i[t];					
		
							check_sgamma[3*st+0] += dr[0] * area_i[t];
							check_sgamma[3*st+1] += dr[1] * area_i[t];
							check_sgamma[3*st+2] += dr[2] * area_i[t];
						}	
						for( int s = 0; s < ntarget; s++ )
						{
							double dr[3] = { sgamma[3*s+0] - check_sgamma[3*s+0],
										sgamma[3*s+1] - check_sgamma[3*s+1],
										sgamma[3*s+2] - check_sgamma[3*s+2] };
							if( dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2] > 1e-4 )
//							printf("%le %le %le versus %le %le %le\n", check_sgamma[3*s+0], check_sgamma[3*s+1], check_sgamma[3*s+2],
//								sgamma[3*s+0], sgamma[3*s+1], sgamma[3*s+2] );
							FCUR -= (sgamma[3*s+0] * sgamma[3*s+0] + sgamma[3*s+1] * sgamma[3*s+1] + sgamma[3*s+2] * sgamma[3*s+2])/srho[s];
							if( aux_i ) FCUR -= aux_w * saux[s] * saux[s] / srho[s]; 				
						}	
					}
#endif
					if( F1 <= F0 && F1 <= F2 )
					{
						states[t2] = s1;
						es = eft_tris+3*t2;

	
						sgamma[3*s2+0] -= dr22[0] * area_i[t2];	
						sgamma[3*s2+1] -= dr22[1] * area_i[t2];	
						sgamma[3*s2+2] -= dr22[2] * area_i[t2];	
						                         
						sgamma[3*s1+0] += dr21[0] * area_i[t2];	
						sgamma[3*s1+1] += dr21[1] * area_i[t2];	
						sgamma[3*s1+2] += dr21[2] * area_i[t2];	
						
						if( aux_i ) saux[s1] += area_i[t2] * aux_i[t2];
						if( aux_i ) saux[s2] -= area_i[t2] * aux_i[t2];
						srho[s2] -= area_i[t2];
						srho[s1] += area_i[t2];
						dF_expected = F1 - F0;
					}
					else
					{
						states[t1] = s2;	
						es = eft_tris+3*t1;	
						
						sgamma[3*s1+0] -= dr11[0] * area_i[t1];	
						sgamma[3*s1+1] -= dr11[1] * area_i[t1];	
						sgamma[3*s1+2] -= dr11[2] * area_i[t1];	
						                         
						sgamma[3*s2+0] += dr12[0] * area_i[t1];	
						sgamma[3*s2+1] += dr12[1] * area_i[t1];	
						sgamma[3*s2+2] += dr12[2] * area_i[t1];	
						
						if( aux_i ) saux[s1] -= area_i[t1] * aux_i[t1];
						if( aux_i ) saux[s2] += area_i[t1] * aux_i[t1];
	
						srho[s1] -= area_i[t1];
						srho[s2] += area_i[t1];
						dF_expected = F2 - F0;
					}
	
#ifndef ALL_EDGES
					for( int ex = 0; ex < 3; ex++ )
					{
						int ne = es[ex];
						if( is_border[ne] ) continue;
	
						int t1 = edges[edge_dim*ne+2];
						int t2 = edges[edge_dim*ne+3];
			
						if( t1 < 0 || t2 <0 ) continue;
	
						if( states[t1] != states[t2] )
						{
							is_border[ne] = 1;
							border_edges[nborder] = ne;
							nborder++;
						}
					}
#endif
	
#ifdef DEBUG_1				
					double FNEW = 0;
					{
						for( int t = 0; t < ntri; t++ )
						{
							int st = states[t];
							double dr[3] = { gamma_i[3*t+0] - sgamma[3*st+0]/srho[st],
									 gamma_i[3*t+1] - sgamma[3*st+1]/srho[st],
									 gamma_i[3*t+2] - sgamma[3*st+2]/srho[st] };
							double put[3];
							MinImage3D( dr, PBC_vec, put );		
							dr[0] += sgamma[3*st+0]/srho[st];
							dr[1] += sgamma[3*st+1]/srho[st];
							dr[2] += sgamma[3*st+2]/srho[st];
							FNEW += area_i[t] *(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);					
							if( aux_i ) FNEW += aux_w * aux_i[t] * aux_i[t] * area_i[t];					
						}	
						for( int s = 0; s < ntarget; s++ )
						{
							FNEW -= (sgamma[3*s+0] * sgamma[3*s+0] + sgamma[3*s+1] * sgamma[3*s+1] + sgamma[3*s+2] * sgamma[3*s+2])/srho[s];
							if( aux_i ) FNEW -= aux_w * saux[s] * saux[s] / srho[s]; 				
						}
					}
#ifdef DEBUG_2
	memcpy( load_cens, sgamma, 3 * ntarget * sizeof(double) );
	for( int s = 0; s < ntarget; s++ )
	{
		load_cens[3*s+0] /= srho[s];
		load_cens[3*s+1] /= srho[s];
		load_cens[3*s+2] /= srho[s];
	}
					double abs_F = compute_F( ntarget, ntri, tris, gamma_i, area_i, states, PBC_vec, load_cens, aux_i, aux_w);
#else
					double abs_F = -1;
#endif
	
					running_F += dF_expected;

					if( fabs(dF_expected-(FNEW-FCUR)) > 1e-4 )
					{
//						printf("PROBLEM_dF_expected: %lf got: %lf FCUR %lf FNEW %lf ABSF %lf dABS: %le F0: %le F1: %le F2: %le\n", dF_expected, FNEW-FCUR, FCUR, FNEW, abs_F, abs_F-abs_F_p,F0, F1, F2  ); 
					}
//					else
//						printf("dF_expected: %lf got: %lf FCUR %lf FNEW %lf ABSF %lf dABS: %le F0: %le F1: %le F2: %le\n", dF_expected, FNEW-FCUR, FCUR, FNEW, abs_F, abs_F-abs_F_p, F0, F1, F2  ); 
#endif
				}
			}
#ifdef ALL_EDGES
		}
#else
		}
#endif

		for( int t = 0; t < ntri; t++ )
		{
			int e_t[3];
			int b = 0;
	
			for( int xt = 0; xt < 3; xt++ )
			{
				int be = eft_tris[3*t+xt];

				if( edges[edge_dim*be+2] == t )
					e_t[xt] = edges[edge_dim*be+3];
				else
					e_t[xt] = edges[edge_dim*be+2];
	
				if( states[t] == states[e_t[xt]] )
					b = 1;
			}

			if( ! b )
			{
				printf("triangle %d has no common borders.\n", t);
			}

		}

	}
	free(sgamma);
	free(srho);
	free(gamma_i);
	free(area_i);
}	

