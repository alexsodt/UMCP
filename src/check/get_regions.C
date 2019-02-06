#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "util.h"
#include "interp.h"
#include "mutil.h"
#include "l-bfgs.h"
#include "alignSet.h"
#include "m_triangles.h"
#include "gauss.h"
#include "config.h"
#include "vpass.h"

#ifdef PARALLEL
void surface::getRegions( int *region_for_tri, int nregions )
{
	int edge_dim = 2;
	
	double *verts;
	int nvert;
	int *tris;
	int *eft_tris;
	int ntri;
	int *edges;
	int nedges;

	getVPassData( &verts, &nvert, &tris, &eft_tris, &ntri, &edges, &nedges, edge_dim );	

	double *r = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	get(r);
	r[3*nv+0] = 1.0;
	r[3*nv+1] = 1.0;
	r[3*nv+2] = 1.0;

	int *states = (int *)malloc( sizeof(int) * ntri );

	int nreg = nregions;

	vpass_inner(
		states,	
		verts, tris, eft_tris, ntri, edges, nedges, 2+edge_dim, nreg, PBC_vec ); 

	memcpy( region_for_tri, states, sizeof(int) * ntri );

}
#endif
