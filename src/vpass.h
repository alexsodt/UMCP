#ifndef __vpassh__
#define __vpassh__
void vpass( double *verts, int *tris, int *eft_tris, int ntri, int *edges, int nedges, int edge_dim, int ntarget, 
		double *out_verts, int **out_tris, int **out_eft_tris, int *out_ntri, int **out_edges, int *out_nedges, double PBC_vec[3][3] );
void cpass( double *gamma_i, double *gamma_nrm, int nv, int *bonds, int nbonds_per_site, int ntarget, double PBC_vec[3][3]);
void vpass_inner(
		int *states, // the output. 
		double *verts, int *tris, int *eft_tris, int ntri, int *edges, int nedges, int edge_dim, int ntarget, double PBC_vec[3][3] ); 
#endif
