
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include <sys/time.h>
#include "interp.h"
#include "mutil.h"
#include "l-bfgs.h"
#include "m_triangles.h"
#include "config.h"
#include "gsl/gsl_sf.h"
#include "gauss.h"
#include "lapack_we_use.h"
#include "Bin.h"
#include "Point.h"
void MinImage( double *dr1, double PBC_vec[3][3], double *put_vec );
/*
int surface::subdivideSurface3( surface *copy_surface, int refine_triangle )
{
	// build a path.
	
	int path_len = 2 * copy_surface->theVertices[refine_at].valence;

	int *path = (int *)malloc( sizeof(int) * path_len );

	int val_cen = copy_surface->theVertices[refine_at].valence;

	vertex * overt = copy_surface->theVertices;

	int err_x = 0;

	int i = refine_at;

	for( int ex = 0; ex < val_cen; ex++ )
	{
		int mid_j = overt[i].edges[ex];

		if( overt[mid_j].valence != 6 )
		{
			err_x = 1;

			return err_x;
		}

		int erev = overt[i].edge_rev[ex];

		erev += 3;
		if( erev >= 6 )
			erev -= 6;

		int j = overt[mid_j].edges[erev]; 	
		int rev_j = overt[mid_j].edge_rev[erev]-1;
		if( rev_j < 0 )
			rev_j += overt[j].valence;
		int next_j = overt[j].edges[rev_j];

		path[2*ex+0] = j;
		path[2*ex+1] = next_j;
	}

	for( int p = 0; p < path_len/2; p++ )
		printf("path %d %d\n", path[2*p+0], path[2*p+1] );

	return subdividePath3( copy_surface, path, path_len );
}*/


int skipvert( vertex *overt, int i, int e, int *the_vert, int *mid, int *erev )
{
	int err_x = 0;

	*mid = overt[i].edges[e];

	if( overt[*mid].valence != 6 )
	{
		err_x = 1;
		return err_x;
	}
	int terev = overt[i].edge_rev[e];
	terev += 3;
	if( terev >= 6 )
		terev -= 6;

	*the_vert = overt[*mid].edges[terev]; 		
	*erev = overt[*mid].edge_rev[terev]; // edge from target vertex that hops us back.

	return err_x;
} 


void surface::subdivideTriangle3( surface *copy_surface, int center_i, int center_e, int *new_center, int *new_edge  )
{
	disable_PBC = 0;
	Sq_cache = NULL;
	cumulative_area = NULL;
	max_valence = 15; // points used.
	opencl_init = 0;
#ifdef FFTW
	h_in = NULL;
	h_out = NULL;
#endif
	bcs = NULL;
	nbc = 0;

	edgeFormulas = NULL;
	theFormulas = NULL;
	theVolumeFormulas = NULL;
	ptBoxes = NULL;
	prBoxes = NULL;
	nf_faces = 0;
	nf_g_q_p = 0;
	nf_irr_faces = 0;
	nf_irr_pts = 0;

	printf("Subdividing.\n");
	c0 = 0;
	PBC_vec[0][0] = copy_surface->PBC_vec[0][0];
	PBC_vec[0][1] =  copy_surface->PBC_vec[0][1];
	PBC_vec[0][2] =  copy_surface->PBC_vec[0][2];
	PBC_vec[1][0] =  copy_surface->PBC_vec[1][0];
	PBC_vec[1][1] =  copy_surface->PBC_vec[1][1];
	PBC_vec[1][2] =  copy_surface->PBC_vec[1][2];
	PBC_vec[2][0] =  copy_surface->PBC_vec[2][0];
	PBC_vec[2][1] =  copy_surface->PBC_vec[2][1];
	PBC_vec[2][2] =  copy_surface->PBC_vec[2][2];

	PBC[0] = copy_surface->PBC[0];
	PBC[1] = copy_surface->PBC[1];
	PBC[2] = copy_surface->PBC[2];
	int path_len = 12;
	int *path = (int *)malloc( sizeof(int) * path_len );
	vertex * overt = copy_surface->theVertices;
	int err_x = 0;

	int center_j, mid_j, j_erev;
	err_x += skipvert( overt, center_i, center_e, &center_j, &mid_j, &j_erev );
	
	int cep1 = center_e + 1;
	if( cep1 >= overt[center_i].valence ) cep1 -= overt[center_i].valence;
	
	int cep2 = cep1 + 1;
	if( cep2 >= overt[center_i].valence ) cep2 -= overt[center_i].valence;

	int cem1 = center_e - 1;
	if( cem1 < 0 ) cem1 += overt[center_i].valence;
	
	int center_k, mid_k, k_erev;
	err_x += skipvert( overt, center_i, cep1, &center_k, &mid_k, &k_erev );
	
	int check_k, mid_jk, check_erev;
	int ejm1 = j_erev-1;
	if( ejm1 < 0 )
		ejm1 += overt[center_j].valence;

	err_x += skipvert( overt, center_j, ejm1, &check_k, &mid_jk, &check_erev );
	
	if( check_k != center_k )
	{
		printf("Logical error 1.\n");
		exit(1);
	}

	int vert_ij, mid_ij, e_ij;
	err_x += skipvert( overt, center_i, cem1, &vert_ij, &mid_ij, &e_ij );
	int e_ijm1 = e_ij-1;
	if( e_ijm1 < 0 ) e_ijm1 += overt[vert_ij].valence;	
	int check_j, mid_ijj;
	err_x += skipvert( overt, vert_ij, e_ijm1, &check_j, &mid_ijj, &check_erev);
	
	if( check_j != center_j )
	{
		printf("Logical error 2.\n");
		exit(1);
	}

	int vert_ik, mid_ik, e_ik;
	err_x += skipvert( overt, center_i, cep2, &vert_ik, &mid_ik, &e_ik );
	int e_ikp1 = e_ik+1;
	if( e_ikp1 > overt[vert_ij].valence ) e_ikp1 -= overt[vert_ik].valence;	
	int mid_ikk, rev_ikk;
	err_x += skipvert( overt, vert_ik, e_ikp1, &check_k, &mid_ikk, &rev_ikk);
	if( check_k != center_k )
	{
		printf("Logical error 3.\n");
		exit(1);
	}
	int eikk_next = rev_ikk + 3;
	if( eikk_next >= overt[center_k].valence )
		eikk_next -= overt[center_k].valence;
	int vert_jk, mid_kkj, rev_jk;
	err_x += skipvert( overt, center_k, eikk_next, &vert_jk, &mid_kkj, &rev_jk );
	rev_jk += 1;
	if( rev_jk >= overt[vert_jk].valence )
		rev_jk -= overt[vert_jk].valence;
	int mid_jjk;
	err_x += skipvert( overt, vert_jk, rev_jk, &check_j, &mid_jjk, &check_erev ); 
	// center_i, center_j, center_k are the unsubdivided triangle we will add resolution to.
	if( check_j != center_j )
	{
		printf("Logical error 4.\n");
		exit(1);
	}

	int mid[3] = { mid_j, mid_jk, mid_k };
	int set[12] = { center_i, mid_ij, vert_ij, mid_ijj, center_j, mid_jjk, vert_jk, mid_kkj, center_k, mid_ikk, vert_ik, mid_ik };
	int newv[12];

	int nvo = copy_surface->nv;

	for( int p = 0; p < 12; p++ )
		newv[p] = nvo+p;

	// manually edit the vertices.

	int new_edges[12][3] = 
	{
		{ newv[1], newv[2], -1 },
		{ newv[2], newv[3], -1 },
		{ newv[3], -1 },
		{ newv[3], newv[8], -1 },
		{ newv[8], newv[10], -1 },
		{ newv[10], newv[11], -1 },
		{ newv[11], -1 },
		{ newv[11], newv[9], -1 },
		{ newv[9], newv[5], -1 },
		{ newv[5], newv[0], -1 },
		{ newv[0], -1 },
		{ newv[0], newv[1], -1 }	
	};
	
	int delete_edges[12][3] = 
	{
		{ mid[0], mid[2], -1 },
		{ mid[0], set[3], -1 },
		{ -1 },
		{ set[1], mid[0], -1 },
		{ mid[0], mid[1], -1 },
		{ mid[1], set[7], -1 },
		{ -1 },
		{ set[5], mid[1], -1 },
		{ mid[1], mid[2], -1 },
		{ mid[2], set[11], -1 },
		{ -1 },
		{ set[9], mid[2], -1 } 
	};

	int new_mid_edges[3][5] =
	{
		{ newv[2], newv[3], newv[8], newv[7], newv[4] },
		{ newv[6], newv[7], newv[10], newv[11], newv[9] },
		{ newv[0], newv[1], newv[4], newv[6], newv[5] }
	};

	int subdivision_origins[12][4] =
	{ // first are 3/8 weights, second are 1/8 weights.
		{ set[10], mid[2], set[0], set[8]	},
		{ set[0], mid[2], mid[0], set[10]	},
		{ set[0], mid[0], set[2], mid[2]	},
		{ mid[0], set[2], set[0], set[4] 	},
		{ mid[0], mid[2], set[0], mid[1]	},
		{ mid[2], set[8], set[10], mid[1]	},
		{ mid[1], mid[2], mid[0], set[8] 	},
		{ mid[0], mid[1], mid[2], set[4]	},
		{ mid[0], set[4], mid[1], set[2] 	},
		{ mid[1], set[8], mid[2], set[6]	},
		{ mid[1], set[4], mid[0], set[6]	},
		{ mid[1], set[6], set[4], set[8]	}
	};

	int all_new_edges[12][6] = 
	{
		{ newv[1], mid[2], newv[5], set[9], set[10], set[11] },
		{ newv[2], newv[4], mid[2], newv[0], set[11], set[0] },
		{ newv[3], mid[0], newv[4], newv[1], set[0], set[1] },
		{ set[1], set[2], set[3], newv[8], mid[0], newv[2] },
		{ newv[1], newv[2], mid[0], newv[7], newv[6], mid[2] },
		{ newv[6], newv[9], set[8], set[9], newv[0], mid[2] },
		{ newv[4], newv[7], mid[1], newv[9], newv[5], mid[2] },
		{ newv[8], newv[10], mid[1], newv[6], newv[4], mid[0] },
		{ newv[3], set[3], set[4], newv[10], newv[7], mid[0] },
		{ newv[5], newv[6], mid[1], newv[11], set[7], set[8] },
		{ mid[1], newv[7], newv[8], set[4], set[5], newv[11] },
		{ mid[1], newv[10], set[5], set[6], set[7], newv[9] }
	};
	
	nv = nvo + 12;	

	theVertices = (vertex *)malloc( sizeof(vertex) * nv );

	/*
		This code block computes the positions of the original indices in the next iteration of the subdivision.
		It uses the formula from Cirak.
	*/

	for( int v = 0; v < copy_surface->nv; v++ )
	{
		int val = copy_surface->theVertices[v].valence;

		theVertices[v].r[0] = copy_surface->theVertices[v].r[0];
		theVertices[v].r[1] = copy_surface->theVertices[v].r[1];
		theVertices[v].r[2] = copy_surface->theVertices[v].r[2];

		theVertices[v].valence = copy_surface->theVertices[v].valence;

		for( int e = 0; e < copy_surface->theVertices[v].valence; e++)
			theVertices[v].edges[e] = copy_surface->theVertices[v].edges[e];
	}
	
	// for the old vertices, delete their stuff and add new edges
	
	for( int sx = 0; sx < 12; sx++ )
	{
		int v = set[sx];

		for( int vx = 0; vx < 3; vx++ )
		{
			if( delete_edges[sx][vx] == -1 ) break;

			for( int ex = 0; ex < theVertices[v].valence; ex++ )
			{
				if( theVertices[v].edges[ex] == delete_edges[sx][vx] )
					theVertices[v].edges[ex] = -1;				
			}
		}
		
		for( int vx = 0; vx < 3; vx++ )
		{
			if( new_edges[sx][vx] == -1 ) break;

			int gotit = 0;
			for( int ex = 0; ex < theVertices[v].valence; ex++ )
			{
				if( theVertices[v].edges[ex] == -1 )
				{
					theVertices[v].edges[ex] = new_edges[sx][vx];
					gotit=1;
					break;
				}				
			}
	
			if( !gotit )
			{
				theVertices[v].edges[theVertices[v].valence] = new_edges[sx][vx];
				theVertices[v].valence += 1;
			}
		}
	}

	for( int mx = 0; mx < 3; mx++ )
	{
		int v= mid[mx];

		theVertices[v].valence=5;

		for( int vx = 0; vx < 5; vx++ )
			theVertices[v].edges[vx] = new_mid_edges[mx][vx]; 
	}
	
	for( int nx = 0; nx < 12; nx++ )
	{
		int v= newv[nx];

		theVertices[v].valence = 6;

		int mix_list[6];
		int nmix=0;

		for( int vx = 0; vx < 6; vx++ )
		{
			theVertices[v].edges[vx] = all_new_edges[nx][vx]; 
			
			if( theVertices[v].edges[vx] < nvo )
			{
				mix_list[nmix] = theVertices[v].edges[vx];
				nmix++;
			}	
		}

#define SUBDIVISION_RULE

#ifdef SUBDIVISION_RULE
		int o = subdivision_origins[nx][0];
		
		double weights[4] = { 3.0/8.0, 3.0/8.0, 1.0/8.0, 1.0/8.0 };

		theVertices[v].r[0] = weights[0] * overt[o].r[0];			
		theVertices[v].r[1] = weights[0] * overt[o].r[1];			
		theVertices[v].r[2] = weights[0] * overt[o].r[2];			
		

		for( int px = 1; px < 4; px++ )
		{
			int n = subdivision_origins[nx][px];
			double dr[3] = { 
				overt[n].r[0] - overt[o].r[0],
				overt[n].r[1] - overt[o].r[1],
				overt[n].r[2] - overt[o].r[2] };
			
			double add[3];
			MinImage(dr, PBC_vec, add );

			theVertices[v].r[0] += (overt[o].r[0]+dr[0])*weights[px];
			theVertices[v].r[1] += (overt[o].r[1]+dr[1])*weights[px];
			theVertices[v].r[2] += (overt[o].r[2]+dr[2])*weights[px];
		}
	
#else
		theVertices[v].r[0] = overt[mix_list[0]].r[0];			
		theVertices[v].r[1] = overt[mix_list[0]].r[1];			
		theVertices[v].r[2] = overt[mix_list[0]].r[2];			

		for( int x = 1; x < nmix; x++ )
		{
			double dr[3] = { 
				overt[mix_list[x]].r[0] - overt[mix_list[0]].r[0],
				overt[mix_list[x]].r[1] - overt[mix_list[0]].r[1],
				overt[mix_list[x]].r[2] - overt[mix_list[0]].r[2] };
			
			double add[3];
			MinImage(dr, PBC_vec, add );

			theVertices[v].r[0] += dr[0]/nmix;
			theVertices[v].r[1] += dr[1]/nmix;
			theVertices[v].r[2] += dr[2]/nmix;
		}
#endif
	}
	
	theTriangles = NULL;
	theEdges = NULL;

	// here I fill out the indices that let us figure out edge indices of our partners	
	for( int i = 0; i < nv; i++ )
	{
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];

			for( int e2 = 0; e2 < theVertices[j].valence; e2++ )
			{
				if( theVertices[j].edges[e2] == i )
					theVertices[i].edge_rev[e] = e2;
			}
		}
	}

	for( int i = 0; i < nv; i++ )
		theVertices[i].c0 = 0;

	orientSurface();
//	constructTriangles();
	assignEdgePBC();	
	sortFaces();	

	*new_center = mid[1];
	*new_edge = 0;

	for( int e = 0; e < theVertices[mid[1]].valence; e++ )
	{
		if( theVertices[mid[1]].edges[e] == newv[7] ) 
			*new_edge = e;
	}
}

int surface::subdividePath3( surface *copy_surface, int *path, int npath )
{
	int complete = 0;
	for( int e = 0; e < copy_surface->theVertices[path[0]].valence; e++ )
	{
		if( copy_surface->theVertices[path[0]].edges[e] == path[npath-1] )
			complete = 1;
	}

	if( !complete )
	{
		printf("Incomplete path.\n");
		return 2;
	}
	else
		printf("Complete path.\n");

	int ninterior = 0;
	int *interior_list;
	
	
	for( int px = 0; px < npath; px++ )
	{
		
	}


	int refine_at = path[0];

	exit(1);
	disable_PBC = 0;
	Sq_cache = NULL;
	cumulative_area = NULL;
	max_valence = 15; // points used.
	opencl_init = 0;
#ifdef FFTW
	h_in = NULL;
	h_out = NULL;
#endif
	bcs = NULL;
	nbc = 0;

	edgeFormulas = NULL;
	theFormulas = NULL;
	theVolumeFormulas = NULL;
	ptBoxes = NULL;
	prBoxes = NULL;
	nf_faces = 0;
	nf_g_q_p = 0;
	nf_irr_faces = 0;
	nf_irr_pts = 0;

	printf("Subdividing.\n");
	c0 = 0;
	PBC_vec[0][0] = copy_surface->PBC_vec[0][0];
	PBC_vec[0][1] =  copy_surface->PBC_vec[0][1];
	PBC_vec[0][2] =  copy_surface->PBC_vec[0][2];
	PBC_vec[1][0] =  copy_surface->PBC_vec[1][0];
	PBC_vec[1][1] =  copy_surface->PBC_vec[1][1];
	PBC_vec[1][2] =  copy_surface->PBC_vec[1][2];
	PBC_vec[2][0] =  copy_surface->PBC_vec[2][0];
	PBC_vec[2][1] =  copy_surface->PBC_vec[2][1];
	PBC_vec[2][2] =  copy_surface->PBC_vec[2][2];

	PBC[0] = copy_surface->PBC[0];
	PBC[1] = copy_surface->PBC[1];
	PBC[2] = copy_surface->PBC[2];

/*
	nvo is the Original Number of Vertices: the number from copy_surface.
*/

	int nvo = copy_surface->nv;

	nv = copy_surface->nv;
	
/*
	we calculate the new number of vertices so I can allocate arrays.
*/

	// we create a new vertex for every valence of the refined location.
	
	nv = nvo + copy_surface->theVertices[refine_at].valence;	

	theVertices = (vertex *)malloc( sizeof(vertex) * nv );

	/*
		This code block computes the positions of the original indices in the next iteration of the subdivision.
		It uses the formula from Cirak.
	*/

	for( int v = 0; v < copy_surface->nv; v++ )
	{
		int val = copy_surface->theVertices[v].valence;

		theVertices[v].r[0] = copy_surface->theVertices[v].r[0];
		theVertices[v].r[1] = copy_surface->theVertices[v].r[1];
		theVertices[v].r[2] = copy_surface->theVertices[v].r[2];

		theVertices[v].valence = copy_surface->theVertices[v].valence;

		for( int e = 0; e < copy_surface->theVertices[v].valence; e++)
			theVertices[v].edges[e] = copy_surface->theVertices[v].edges[e];
	}

	// tnv is the index of the created indices. It starts at nvo, the index of the first new vertex.
	int tnv = nvo;	
	
	// this code block fills out some initial values for the new vertices.
	// *all* of the edge indexes for the original vertices change : they point to new subdivision points.
	// we fill these in as we create the vertices.
	
	int old_val = theVertices[refine_at].valence;	

	int list[old_val];
	int newv[old_val];
	for( int ex = 0; ex < old_val; ex++ )
	{
		int i = refine_at;
		int exp1 = ex+1;
		if( exp1 >= theVertices[i].valence )
			exp1 -= theVertices[i].valence;

		int j = copy_surface->theVertices[i].edges[ex];
		int k = copy_surface->theVertices[i].edges[exp1];

		list[ex] = j; // list of the old edges.
		newv[ex] = tnv;

		theVertices[tnv].r[0] = theVertices[i].r[0]; 
		theVertices[tnv].r[1] = theVertices[i].r[1]; 
		theVertices[tnv].r[2] = theVertices[i].r[2]; 
		
		int indos[2] = { j, k };
	
		for( int ix = 0; ix < 2; ix++ )
		{	
			double dr[3] = { theVertices[indos[ix]].r[0] - theVertices[i].r[0],
					 theVertices[indos[ix]].r[1] - theVertices[i].r[1],
					 theVertices[indos[ix]].r[2] - theVertices[i].r[2] };
			double add[3];
			MinImage(dr, PBC_vec, add );

			theVertices[tnv].r[0] += dr[0]/3;
			theVertices[tnv].r[1] += dr[1]/3;
			theVertices[tnv].r[2] += dr[2]/3;
		}

		theVertices[i].edges[ex] = tnv;

		tnv++;
	}

	// redo the edges of the outer (old) triangles.

	for( int ex = 0; ex < old_val; ex++ )
	{
	 	int ind_rep = -1;

		int j = list[ex];

		for( int e = 0; e < theVertices[j].valence; e++ )
		{
			if( theVertices[j].edges[e] == refine_at )
				ind_rep = e;
		}

		int copy[theVertices[j].valence+1];

		memcpy( copy, theVertices[j].edges, sizeof(int) * theVertices[j].valence );

		int put = 0;

		for( int e = 0 ; e < theVertices[j].valence; e++ )
		{	
			if( e == ind_rep )
			{
				int exp1 = ex+1;
				if( exp1 >= old_val )
					exp1 -= old_val;
				int exm1 = ex-1;
				if( exm1< 0 )
					exm1 += old_val;

				theVertices[j].edges[put++] = newv[ex]; 
				theVertices[j].edges[put++] = newv[exm1];
			}
			else
			{
				theVertices[j].edges[put++] = copy[e]; 
			}
		}

		theVertices[j].valence++;
	}

	// create the vertices of the interior triangles.

	for( int ex = 0; ex < old_val; ex++ )
	{
		int m = newv[ex];
		int exp1 = ex+1; if( exp1 >= old_val ) exp1 -= old_val;
		int exm1 = ex-1; if( exm1 < 0 ) exm1 += old_val;

		theVertices[m].edges[0] = refine_at;
		theVertices[m].edges[1] = newv[exm1];
		theVertices[m].edges[2] = list[ex];
		theVertices[m].edges[3] = list[exp1];
		theVertices[m].edges[4] = newv[exp1];

		theVertices[m].valence = 5;
	}
	
	theTriangles = NULL;
	theEdges = NULL;



	// here I fill out the indices that let us figure out edge indices of our partners	
	for( int i = 0; i < nv; i++ )
	{
		for( int e = 0; e < theVertices[i].valence; e++ )
		{
			int j = theVertices[i].edges[e];

			for( int e2 = 0; e2 < theVertices[j].valence; e2++ )
			{
				if( theVertices[j].edges[e2] == i )
					theVertices[i].edge_rev[e] = e2;
			}
		}
	}

	for( int i = 0; i < nv; i++ )
		theVertices[i].c0 = 0;
	orientSurface();
//	constructTriangles();
	assignEdgePBC();	
	sortFaces();	
}
