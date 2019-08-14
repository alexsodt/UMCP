#define __parallelc__

//#define USE_DGEMV
#define PRINT_REGIONS

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "parallel.h"
#include "interp.h"
#include "pcomplex.h"
#include "lapack_we_use.h"
#include <math.h>
#include "simulation.h"

// global info
parallel_info par_info; 
int setup_for_parallel = 0;

void quietParallel( void )
{
	int ierr;
	int nprocs, taskid;
#ifdef PARALLEL
         /* Change the random number if this is part of a parallel run */
         if ( (ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs)) ) {
           exit(1);
         }
         if ( (ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid)) ) {
           exit(1);
         }

	if( taskid != BASE_TASK )
	{
		freopen( "/dev/null","w", stdout );
	}
#endif
}

void balanceWork( surface *theSurface, int *regions_for_tri )
{
	double cur_work[par_info.nprocs];
	memset( cur_work, 0, sizeof(double) * par_info.nprocs );
	int nt = theSurface->nt;
	for( int t = 0; t < nt; t++ )
	{
		if( theSurface->theTriangles[t].f < theSurface->nf_faces )
			cur_work[regions_for_tri[t]] += theSurface->nf_g_q_p;
		else
			cur_work[regions_for_tri[t]] += theSurface->nf_irr_pts;
	}
	
	int done = 0;

	triangle *tris = theSurface->theTriangles;

	while( !done )
	{
		done = 1;

		for( int t = 0; t < nt; t++ )
		{	
			int face = tris[t].f;

			int my_proc = regions_for_tri[t];
			double the_work = theSurface->nf_g_q_p;
			if( face >= theSurface->nf_faces )
				the_work = theSurface->nf_irr_pts;

			int procs[3] = { 
				regions_for_tri[tris[t].border_tri[0]], 
				regions_for_tri[tris[t].border_tri[1]], 
				regions_for_tri[tris[t].border_tri[2]] };

			int alt[par_info.nprocs];
			memset( alt, 0, sizeof(int) * par_info.nprocs );
			 
			if( procs[0] != my_proc )
				alt[procs[0]]++;
			if( procs[1] != my_proc )
				alt[procs[1]]++;
			if( procs[2] != my_proc )
				alt[procs[2]]++;

//			for( int x = 0; x < par_info.nprocs; x++ )
//				printf(" %d ", alt[x] );
//			printf("\n");

			for( int x = 0; x < par_info.nprocs; x++ )
			{
			//	if( alt[x] == 2 )
				{
//					printf("Checking %lf and %lf\n", fabs((cur_work[x]+the_work)/(cur_work[my_proc]-the_work+1e-10)-1),
//								 fabs(cur_work[x] / (cur_work[my_proc]+1e-10)-1) );
				}

				if(  alt[x] > 0 && fabs((cur_work[x]+the_work)/(cur_work[my_proc]-the_work+1e-10)-1) < fabs(cur_work[x] / (cur_work[my_proc]+1e-10)-1) )
				{
//					printf("Transferring %lf from %d (%lf) to %d (%lf)\n",
//						the_work, my_proc, cur_work[my_proc], x, cur_work[x] ); 
					done=0;
					cur_work[x]       += the_work;
					cur_work[my_proc] -= the_work;
					regions_for_tri[t] = x;
					break;
				}
			}
		}
	}
	
	double check_work[par_info.nprocs];
	memset( check_work, 0, sizeof(double) * par_info.nprocs );

	for( int t = 0; t < nt; t++ )
	{
		if( theSurface->theTriangles[t].f < theSurface->nf_faces )
			check_work[regions_for_tri[t]] += theSurface->nf_g_q_p;
		else
			check_work[regions_for_tri[t]] += theSurface->nf_irr_pts;
	}


#if 0
	printf("\"Balanced\" work:\n");

	for( int x = 0; x < par_info.nprocs; x++ )
		printf(" %lf", cur_work[x] );
	printf("\n");
	printf("Check\n");
	for( int x = 0; x < par_info.nprocs; x++ )
		printf(" %lf", check_work[x] );
	printf("\n");
#endif
//	exit(1);
	
}

void setupParallel( Simulation *theSimulation )
{
	int nTotFaces = 0;
	int nSurfaces = 0;

	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		nTotFaces += sRec->theSurface->nt;
		nSurfaces++;
	}
	setup_for_parallel = 1;

	par_info.nSurfaces = nSurfaces;	
	par_info.faces         = (int **)malloc( sizeof(int *) * nSurfaces );
	par_info.nf            = (int *)malloc( sizeof(int) * nSurfaces );
	par_info.verts         = (int **)malloc( sizeof(int *) * nSurfaces );
	par_info.nv            = (int *)malloc( sizeof(int) * nSurfaces );
	par_info.proc_for_vert = (int **)malloc( sizeof(int *) * nSurfaces );
	par_info.genQ          = (int **)malloc( sizeof(int *) * nSurfaces );
	par_info.NQ            = (int *)malloc( sizeof(int) * nSurfaces );
 
	int *ntri = (int *)malloc( sizeof(int) * nSurfaces);
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		ntri[sRec->id] = sRec->theSurface->nt;
		par_info.faces[sRec->id] = (int *)malloc( sizeof(int) * sRec->theSurface->nt ); 
		par_info.nf[sRec->id] = 0;	
		
		par_info.verts[sRec->id] = (int *)malloc( sizeof(int) * sRec->theSurface->nv ); 
		par_info.nv[sRec->id] = 0;	

		par_info.proc_for_vert[sRec->id] = (int *)malloc( sizeof(int) * sRec->theSurface->nv ); 
	}

	if( theSimulation->ncomplex > 0 )		
		par_info.complexes = (int *)malloc( sizeof(int) * theSimulation->ncomplex );
	else	
		par_info.complexes = NULL;
	par_info.nc = 0;

	// real parallel partition

	int ierr=0;
	// defaults for not parallel.
	int nprocs=1, taskid=0;

#ifdef PARALLEL
         /* Change the random number if this is part of a parallel run */
         if ( (ierr = MPI_Comm_size(MPI_COMM_WORLD, &nprocs)) ) {
           exit(1);
         }
         if ( (ierr = MPI_Comm_rank(MPI_COMM_WORLD, &taskid)) ) {
           exit(1);
         }
#endif
	
	par_info.nprocs = nprocs;
	par_info.my_id = taskid;

	int ncomplex = theSimulation->ncomplex;
	pcomplex **allComplexes = theSimulation->allComplexes;

	if( par_info.nprocs == 1 && 0 )
	{
		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		{
			int id = sRec->id;
			surface *theSurface = sRec->theSurface;
			int nt = theSurface->nt;
			int nv = theSurface->nv;
	
			par_info.nf[id] = theSurface->nt;
			for( int f = 0; f < nt; f++ )
				par_info.faces[id][f] = f;
			par_info.nf[id] = nt;
			
			par_info.nv[id]     = nv;
		
			for( int v = 0; v < nv; v++ )
				par_info.verts[id][v] = v;
		}
			
		for( int c = 0; c < ncomplex; c++ )
			par_info.complexes[c] = c;
		par_info.nc = ncomplex;
	}
	else
	{	
		// partition processors between the surfaces.
		// we need to judge the approximate workload.

		int **regions_for_tri = (int **)malloc( sizeof(int*) * nSurfaces );
		for( int s = 0; s < nSurfaces; s++ )
			regions_for_tri[s] = (int *)malloc( sizeof(int) * ntri[s] );

		if( par_info.my_id == BASE_TASK )
		{ // RUN IN SERIAL NOW.
			for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
			{ // balance each surface individually, that's the work unit.
				sRec->theSurface->getRegions( regions_for_tri[sRec->id], par_info.nprocs );
				balanceWork( sRec->theSurface, regions_for_tri[sRec->id] );
			}
		}
#ifdef PARALLEL
		for( int s = 0; s < nSurfaces; s++ )
			MPI_Bcast( regions_for_tri[s], ntri[s], MPI_INT, BASE_TASK, MPI_COMM_WORLD );
#endif

		for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		{
			surface *theSurface = sRec->theSurface;
			int nv = theSurface->nv;
			int nt = theSurface->nt;
			int id = sRec->id;
			int NQ = sRec->NQ;

			int *vclaim = (int *)malloc( sizeof(int) * MAX_VALENCE * nv );
			int *nclaim = (int *)malloc( sizeof(int) * nv );
			memset( nclaim, 0, sizeof(int) * nv );

			for( int t = 0; t < nt; t++ )
			{
				for( int ix = 0; ix < 3; ix++ )
				{
					int v = theSurface->theTriangles[t].ids[ix];

					vclaim[v*MAX_VALENCE+nclaim[v]] = regions_for_tri[id][t];
					nclaim[v]++;
				}
			}	
			
			int *proc_for_vertex = (int *)malloc( sizeof(int) * nv );
				
			for( int v = 0; v < nv; v++ )
			{
				int val = theSurface->theVertices[v].valence;

				int nclaims[val];
				int pclaims[val];

				memset( nclaims, 0, sizeof(int) * val );
				for( int x = 0; x < val; x++ )
					pclaims[x] = -1;

				int claimed_by = -1;

				for( int x = 0; x < nclaim[v] && claimed_by == -1; x++ )
				{
					int cl = vclaim[v*MAX_VALENCE+x];

					int ind = -1;

					for( int p = 0; p < val;p ++ )
					{
						if( pclaims[p] == cl )
							ind = p;
						else if( pclaims[p] == -1 )
							ind = p;
						if( ind != -1 ) break;
					}

					pclaims[ind] = cl;
					nclaims[ind]++;

					if( nclaims[ind] == 2 )
						claimed_by = pclaims[ind];
				}

				if( claimed_by == -1 )
					claimed_by = vclaim[v*MAX_VALENCE+0];
				
				proc_for_vertex[v] = claimed_by;
			}	

#ifdef PARALLEL
			MPI_Bcast( proc_for_vertex, nv, MPI_INT, BASE_TASK, MPI_COMM_WORLD );
#endif
			par_info.proc_for_vert[id] = proc_for_vertex;

			free(vclaim);
			free(nclaim);

			for( int t = 0; t < nt; t++ )
			{
				if( regions_for_tri[id][t] == par_info.my_id )
				{
					par_info.faces[id][par_info.nf[id]] = theSurface->theTriangles[t].f;
					par_info.nf[id]++;
				}
			}
			
			for( int v = 0; v < nv; v++ )
			{
				if( proc_for_vertex[v] == par_info.my_id )
				{
					par_info.verts[id][par_info.nv[id]] = v;
					par_info.nv[id]++;
				}
			}

			if( NQ == -1 ) //vertices are the general coordinates.
			{
				par_info.genQ[id] = (int *)malloc( sizeof(int) * 3 * par_info.nv[id] );

				for( int vx = 0; vx < par_info.nv[id]; vx++ )
				{
					int v = par_info.verts[id][vx];

					par_info.genQ[id][3*vx+0] = 3*v+0;
					par_info.genQ[id][3*vx+1] = 3*v+1;
					par_info.genQ[id][3*vx+2] = 3*v+2;
				}
			}
			else
			{
				for( int pass = 0; pass < 2; pass++ )
				{
					if( pass == 1 ) par_info.genQ[id] = (int *)malloc( sizeof(int) * par_info.NQ[id] );

					par_info.NQ[id] = 0;

					for( int vq = 0; vq < NQ; vq++ )
					{
						//if( vq % par_info.nprocs == par_info.my_id )
						if( par_info.my_id == BASE_TASK )
						{
							if( pass == 1 ) par_info.genQ[id][par_info.NQ[id]] = vq;
							par_info.NQ[id] += 1;
						}
					}
				}
			}


		}

		for( int c = 0; c < ncomplex; c++ )
		{
			if( c % par_info.nprocs == par_info.my_id )
			{
				par_info.complexes[par_info.nc] = c;
				par_info.nc++;
			}
		}
	}

	free(ntri);

}	


void setupSparseVertexPassing( Simulation *theSimulation )
//SparseMatrix *EFFM, int nv, int do_gen_q )
{
#ifdef PARALLEL
	int nSurfaces = 0;
	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
		nSurfaces++;

	par_info.nget_from = (int **)malloc( sizeof(int*) * nSurfaces );
	par_info.nsend_to = (int **)malloc( sizeof(int*) * nSurfaces );
	par_info.send_to = (int ***)malloc( sizeof(int**) * nSurfaces );
	par_info.get_from = (int ***)malloc( sizeof(int**) * nSurfaces );
	par_info.psum_nget_from = (int **)malloc( sizeof(int*) * nSurfaces );
	par_info.psum_nsend_to = (int **)malloc( sizeof(int*) * nSurfaces );
	par_info.psum_send_to = (int ***)malloc( sizeof(int**) * nSurfaces );
	par_info.psum_get_from = (int ***)malloc( sizeof(int**) * nSurfaces );

	for( surface_record *sRec = theSimulation->allSurfaces; sRec; sRec = sRec->next )
	{
		surface *theSurface = sRec->theSurface;	
		int nv = theSurface->nv;
		int do_gen_q = sRec->do_gen_q;
		int id = sRec->id;
		SparseMatrix *EFFM = sRec->EFFM;

		MPI_Status status;	
		int np = par_info.nprocs;
			
		int *need_list = (int *)malloc( sizeof(int) * nv );

		// if generalized coordinates are being used, the base task needs information from all processors.

		int *all_list = (int *)malloc( sizeof(int) * nv );
		for( int v = 0; v < nv; v++ )
			all_list[v] = v;
			

		for( int pass = 0; pass < 2; pass++ )
		{
			int nsend_mat[np*np];
			int nget_mat[np*np];
		
			// fetch what is needed from each processor.
			int *source_list = EFFM->source_list;
			int nsrc = EFFM->nsource;
 

			if( do_gen_q && par_info.my_id == BASE_TASK )
			{
				nsrc = nv;
				source_list = all_list;
			}	
			else if( do_gen_q ) // supporting processors do not require any vertices.
				nsrc = 0;

			// in "parallel sum" mode, we pass whatever is in our source_list, figuring we are likely to compute this as part of our gradient.

			if( pass == 0 )
			{	
				par_info.nget_from[id] = (int *)malloc( sizeof(int) * par_info.nprocs );
				par_info.nsend_to[id] = (int *)malloc( sizeof(int) * par_info.nprocs );
				par_info.send_to[id] = (int **)malloc( sizeof(int*) * par_info.nprocs );
				par_info.get_from[id] = (int **)malloc( sizeof(int*) * par_info.nprocs );
			}
			else if( pass == 1 )
			{
				par_info.psum_nget_from[id] = (int *)malloc( sizeof(int) * par_info.nprocs );
				par_info.psum_nsend_to[id] = (int *)malloc( sizeof(int) * par_info.nprocs );
				par_info.psum_send_to[id] = (int **)malloc( sizeof(int*) * par_info.nprocs );
				par_info.psum_get_from[id] = (int **)malloc( sizeof(int*) * par_info.nprocs );
			}
			
			for( int p1 = 0; p1 < par_info.nprocs; p1++ )
			{
				for( int p2 = 0; p2 < par_info.nprocs; p2++ )
				{
					if( p1 == p2 ) continue;
		
					if( par_info.my_id == p1 )
					{	// this is what I need from you.
		
						int nneed = 0;
						for( int x = 0; x < nsrc; x++ )
						{
							if( pass == 0 && par_info.proc_for_vert[id][source_list[x]] == p2 )
							{ 	// pass only what we know it provides as its main vertex.
								need_list[nneed] = source_list[x];
								nneed++;
							}
							else if( pass == 1 )
							{	// pass everything as a tentative request.
								need_list[nneed] = source_list[x];
								nneed++;
							}
						}					

						// send my need list to the other processor...			

						MPI_Send( &nneed, 1, MPI_INT, p2, (p1 * np + p2)*6+(pass==0?0:2),  MPI_COMM_WORLD );
		
						if( pass == 0 )
							par_info.nget_from[id][p2] = nneed;
						else
							par_info.psum_nget_from[id][p2] = nneed;
		
						if( nneed > 0 )
						{
							MPI_Send( need_list, nneed, MPI_INT, p2, (p1*np+p2)*6+(pass==0?1:3), MPI_COMM_WORLD );
							if( pass == 0 )
							{
								par_info.get_from[id][p2] = (int *)malloc( sizeof(int) * nneed );
								memcpy( par_info.get_from[id][p2], need_list, sizeof(int) * nneed );
								printf("partial %d getting %d from %d.\n", p1, nneed, p2 );
							}
							else
							{
								// We receive the truncated list.
								MPI_Recv( &nneed, 1, MPI_INT, p2, (p1 * np + p2)*6+4,  MPI_COMM_WORLD, &status );
								MPI_Recv( need_list, nneed, MPI_INT, p2, (p1*np+p2)*6+5, MPI_COMM_WORLD, &status );
								par_info.psum_nget_from[id][p2] = nneed;
								par_info.psum_get_from[id][p2] = (int *)malloc( sizeof(int) * nneed );
								memcpy( par_info.psum_get_from[id][p2], need_list, sizeof(int) * nneed );
							}
						}
					}
		
					if( par_info.my_id == p2 )
					{
						// this is what you need from me.
						int nsend = 0;
		
						MPI_Recv( &nsend, 1, MPI_INT, p1, (p1*np+p2)*6+(pass==0?0:2), MPI_COMM_WORLD, &status );
						if( pass == 0 )
						{
							par_info.nsend_to[id][p1] = nsend;
			
							if( nsend > 0 )
							{
								par_info.send_to[id][p1] = (int *)malloc( sizeof(int) * nsend );
								printf("%d RECEIVING %d ints from %d.\n", p2, nsend, p1 );
								MPI_Recv( par_info.send_to[id][p1], nsend, MPI_INT, p1, (p1*np+p2)*6+1, MPI_COMM_WORLD, &status );
									
							}
						}
						else if( pass == 1 )
						{
							int full_request = nsend;
							int n_provided = 0;

							par_info.psum_nsend_to[id][p1] = nsend;
			
							if( nsend > 0 )
							{
								int *could_pass = (int *)malloc( sizeof(int) * nsend );
								par_info.psum_send_to[id][p1] = (int *)malloc( sizeof(int) * nsend );
									
								MPI_Recv(could_pass, nsend, MPI_INT, p1, (p1*np+p2)*6+3, MPI_COMM_WORLD, &status );
								// loop over my source list, which is everything I could reasonably provide a gradient for.
								// if it matches something the other processor needs a gradient for, we must pass it.

								int npass = 0;

								for( int x = 0; x < nsrc; x++ )
								{
									int pass_it = 0;

									for( int y = 0; y < nsend; y++ )
									{
										if( source_list[x] == could_pass[y] )
										{
											pass_it = 1;
											break;
										}
									}

									if( pass_it )
									{
										par_info.psum_send_to[id][p1][npass] = source_list[x];
										npass++;
									}
								}

								par_info.psum_nsend_to[id][p1] = npass;
								free(could_pass);
						
								MPI_Send( &npass, 1, MPI_INT, p1, (p1 * np + p2)*6+4,  MPI_COMM_WORLD );
								MPI_Send( par_info.psum_send_to[id][p1], npass, MPI_INT, p1, (p1 * np + p2)*6+5,  MPI_COMM_WORLD );
							}
						}
					}
				}
			}
		}
		
		free(all_list);
		free(need_list);
	}
#endif
}

void surface::getVPassData( double **verts, int *nvert, int **tris, int ** eft_tris, int *ntri, int **edges, int *nedges_in, int edge_dim )
{
	*nvert = nv;

	(*verts) = (double *)malloc( sizeof(double) * 3 * nv );

	for( int v = 0; v < nv; v++ )
	{
		(*verts)[3*v+0] = theVertices[v].r[0];
		(*verts)[3*v+1] = theVertices[v].r[1];
		(*verts)[3*v+2] = theVertices[v].r[2];
	}

	(*tris) = (int *)malloc( sizeof(int) * 3 * nt );
	(*eft_tris) = (int *)malloc( sizeof(int) * 3 * nt );

	*ntri = nt;

	for( int t = 0; t < nt; t++ )
	{
		(*tris)[3*t+0] = theTriangles[t].ids[0];
		(*tris)[3*t+1] = theTriangles[t].ids[1];
		(*tris)[3*t+2] = theTriangles[t].ids[2];

		(*eft_tris)[3*t+0] = theTriangles[t].edges[0];
		(*eft_tris)[3*t+1] = theTriangles[t].edges[1];
		(*eft_tris)[3*t+2] = theTriangles[t].edges[2];
	}

	*(nedges_in) = nedges;
	
	(*edges) = (int *)malloc( sizeof(int) * (2+edge_dim) * nedges );	

	for( int e = 0; e < nedges; e++ )
	{
		(*edges)[e*(2+edge_dim)+0] = theEdges[e].vertices[0];
		(*edges)[e*(2+edge_dim)+1] = theEdges[e].vertices[1];
		(*edges)[e*(2+edge_dim)+2] = theEdges[e].faces[0];
		(*edges)[e*(2+edge_dim)+3] = theEdges[e].faces[1];
	}
}

void ParallelGather( int *vec, int len_per_proc )
{
#ifdef PARALLEL
	int *recv=NULL;
	if( par_info.my_id == BASE_TASK )
		recv=  (int *)malloc( sizeof(int) * len_per_proc * par_info.nprocs );
	MPI_Gather( vec+par_info.my_id*len_per_proc, len_per_proc, MPI_INT, recv, len_per_proc, MPI_INT, BASE_TASK, MPI_COMM_WORLD );

	if( par_info.my_id == BASE_TASK )
	{
		memcpy( vec, recv, sizeof(int) * len_per_proc * par_info.nprocs);
		free(recv);
	}
#endif
}

void ParallelGather( double *vec, int len_per_proc )
{
#ifdef PARALLEL
	double *recv=NULL;
	if( par_info.my_id == BASE_TASK )
		recv=  (double *)malloc( sizeof(double) * len_per_proc * par_info.nprocs );
	MPI_Gather( vec+par_info.my_id*len_per_proc, len_per_proc, MPI_DOUBLE, recv, len_per_proc, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );

	if( par_info.my_id == BASE_TASK )
	{
		memcpy( vec, recv, sizeof(double) * len_per_proc * par_info.nprocs);
		free(recv);
	}
#endif
}

void ParallelSum( int *vec, int len )
{
#ifdef PARALLEL
	int *vin=(int*)malloc(sizeof(int)*len);
	memset(vin,0,sizeof(int)*len);
	MPI_Allreduce( vec, vin, len, MPI_INT, MPI_SUM, MPI_COMM_WORLD ); 
	memcpy( vec, vin, sizeof(int) * len );
	free(vin);
#endif
}

void ParallelSum( double *vec, int len )
{
#ifdef PARALLEL
	double *vin=(double*)malloc(sizeof(double)*len);
	memset(vin,0,sizeof(double)*len);
	MPI_Allreduce( vec, vin, len, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD ); 
	memcpy( vec, vin, sizeof(double) * len );
	free(vin);
#endif
}

void ParallelSyncComplexes( pcomplex **allComplexes, int ncomplexes )
{
#ifdef PARALLEL
	/* I'm sure this is awful and needs to be rewritten but for now it just needs to work, -alex */

	int ispace = 100;
	int fspace = 100;

	int *ibuffer = (int *)malloc( sizeof(int) * ispace );
	double *fbuffer = (double *)malloc( sizeof(double) * fspace );

	for( int p = 0; p < par_info.nprocs; p++ )
	{
		int nrec = 0;
		if( p == par_info.my_id )
			nrec = par_info.nc;
		MPI_Bcast( &nrec, 1, MPI_INT, p, MPI_COMM_WORLD );
	
		if( nrec > 0 )
		{	
			int *list = (int *)malloc( sizeof(int) * nrec );
			
			if( p == par_info.my_id )
				memcpy( list, par_info.complexes, sizeof(int) * nrec );

			MPI_Bcast( list, nrec, MPI_INT, p, MPI_COMM_WORLD );

			// broadcast positions and gradients.
		
			int len = 0;

			for( int cx = 0; cx < nrec; cx++ )
				len += allComplexes[list[cx]]->packLenF();
		
			if( len > fspace )
			{
				fspace = len;
				fbuffer = (double *)realloc( fbuffer, sizeof(double) * fspace );

			}
			
			if( p == par_info.my_id )
			{
				int io = 0;
				for( int cx = 0; cx < par_info.nc; cx++ )
				{
					allComplexes[list[cx]]->packageF(fbuffer+io);
					io += allComplexes[list[cx]]->packLenF();
				}
			}

			MPI_Bcast( fbuffer, len, MPI_DOUBLE, p, MPI_COMM_WORLD );
				
			int io = 0;

			for( int cx = 0; cx < nrec; cx++ )
			{
				allComplexes[list[cx]]->unpackF(fbuffer+io);
				io += allComplexes[list[cx]]->packLenF();
			}
			
			len = 0;
			for( int cx = 0; cx < nrec; cx++ )
				len += allComplexes[list[cx]]->packLenI();
		
			if( len > ispace )
			{
				ispace = len;
				ibuffer = (int *)realloc( ibuffer, sizeof(double) * ispace );

			}
			
			if( p == par_info.my_id )
			{
				int io = 0;
				for( int cx = 0; cx < par_info.nc; cx++ )
				{
					allComplexes[list[cx]]->packageI(ibuffer+io);
					io += allComplexes[list[cx]]->packLenI();
				}
			}

			MPI_Bcast( ibuffer, len, MPI_INT, p, MPI_COMM_WORLD );
				
			io = 0;

			for( int cx = 0; cx < nrec; cx++ )
			{
				allComplexes[list[cx]]->unpackI(ibuffer+io);
				io += allComplexes[list[cx]]->packLenI();
			}

			free(list);	
		}
	} 	

	free(ibuffer);
	free(fbuffer);
#endif
}

void FullSyncGenQ( double *total_vec, int surface_id )
{
#ifdef PARALLEL
	int vspace = 100*3;

	double *vec = (double *)malloc( sizeof(double) * vspace );
	int    *ivec = (int *)malloc( sizeof(int) * vspace );

	for( int p = 0; p < par_info.nprocs; p++ )
	{
		int s = surface_id;
		int NQ = par_info.NQ[s];

		MPI_Bcast( &NQ, 1, MPI_INT, p, MPI_COMM_WORLD );

		if( NQ*3 > vspace )
		{
			vspace = NQ * 3;
			vec = (double *)realloc( vec, sizeof(double) * vspace );
			ivec = (int *)realloc( ivec, sizeof(double) * vspace );
		}
		if( p == par_info.my_id )
			memcpy( ivec, par_info.genQ[s], sizeof(int) * par_info.NQ[s] );
			
		MPI_Bcast( ivec, NQ, MPI_INT, p, MPI_COMM_WORLD );

		if( p == par_info.my_id )
		{
			for( int x = 0; x < NQ; x++ )
				vec[x] = total_vec[par_info.genQ[s][x]];
		}
		
		MPI_Bcast( vec, NQ, MPI_DOUBLE, p, MPI_COMM_WORLD );

		for( int x = 0; x < NQ; x++ )
		{
			total_vec[ivec[x]] = vec[x];
		}
	}

	free(vec);
	free(ivec);
#endif
}

void GatherVertices( double *total_vec, int surface_id)
{
#ifdef PARALLEL
	int vspace = 100*3;

	double *vec = (double *)malloc( sizeof(double) * vspace );
	int    *ivec = (int *)malloc( sizeof(int) * vspace );

	for( int p = 0; p < par_info.nprocs; p++ )
	{
		int nv = par_info.nv[surface_id];

		MPI_Bcast( &nv, 1, MPI_INT, p, MPI_COMM_WORLD );

		if( nv*3 > vspace )
		{
			vspace = nv * 3;
			vec = (double *)realloc( vec, sizeof(double) * vspace );
			ivec = (int *)realloc( ivec, sizeof(double) * vspace );
		}
		if( p == par_info.my_id )
			memcpy( ivec, par_info.verts[surface_id], sizeof(int) * par_info.nv[surface_id] );
			
		MPI_Bcast( ivec, nv, MPI_INT, p, MPI_COMM_WORLD );

		if( p == par_info.my_id )
		{
			for( int x = 0; x < nv; x++ )
			{
				vec[x*3+0] = total_vec[par_info.verts[surface_id][x]*3+0];
				vec[x*3+1] = total_vec[par_info.verts[surface_id][x]*3+1];
				vec[x*3+2] = total_vec[par_info.verts[surface_id][x]*3+2];
			}
		}
		
		MPI_Bcast( vec, nv*3, MPI_DOUBLE, p, MPI_COMM_WORLD );

		for( int x = 0; x < nv; x++ )
		{
			total_vec[ivec[x]*3+0] = vec[3*x+0];
			total_vec[ivec[x]*3+1] = vec[3*x+1];
			total_vec[ivec[x]*3+2] = vec[3*x+2];
		}
	}

	free(vec);
	free(ivec);
#endif
}

void FullSyncVertices( double *total_vec, int surface_id )
{
#ifdef PARALLEL
	int vspace = 100*3;

	double *vec = (double *)malloc( sizeof(double) * vspace );
	int    *ivec = (int *)malloc( sizeof(int) * vspace );

	for( int p = 0; p < par_info.nprocs; p++ )
	{
		int nv = par_info.nv[surface_id];

		MPI_Bcast( &nv, 1, MPI_INT, p, MPI_COMM_WORLD );

		if( nv*3 > vspace )
		{
			vspace = nv * 3;
			vec = (double *)realloc( vec, sizeof(double) * vspace );
			ivec = (int *)realloc( ivec, sizeof(double) * vspace );
		}
		if( p == par_info.my_id )
			memcpy( ivec, par_info.verts[surface_id], sizeof(int) * par_info.nv[surface_id] );
			
		MPI_Bcast( ivec, nv, MPI_INT, p, MPI_COMM_WORLD );

		if( p == par_info.my_id )
		{
			for( int x = 0; x < nv; x++ )
			{
				vec[x*3+0] = total_vec[par_info.verts[surface_id][x]*3+0];
				vec[x*3+1] = total_vec[par_info.verts[surface_id][x]*3+1];
				vec[x*3+2] = total_vec[par_info.verts[surface_id][x]*3+2];
			}
		}
		
		MPI_Bcast( vec, nv*3, MPI_DOUBLE, p, MPI_COMM_WORLD );

		for( int x = 0; x < nv; x++ )
		{
			total_vec[ivec[x]*3+0] = vec[3*x+0];
			total_vec[ivec[x]*3+1] = vec[3*x+1];
			total_vec[ivec[x]*3+2] = vec[3*x+2];
		}
	}

	free(vec);
	free(ivec);
#endif
}

void ParallelBroadcast( double *vec, int len )
{
#ifdef PARALLEL
	MPI_Bcast( vec, len, MPI_DOUBLE, BASE_TASK, MPI_COMM_WORLD );
#endif
}

void ParallelBroadcast( int *vec, int len )
{
#ifdef PARALLEL
	MPI_Bcast( vec, len, MPI_INT, BASE_TASK, MPI_COMM_WORLD );
#endif
}

void CartMatVecIncrScale( double *vec_out, double *vec_in, double *Mat, double scale, int nv, double *alphas, int surface_id)
{
	char trans = 'N';
	int n = nv;
	double one = 1.0;
	double zero = 0.0;
	int incrxy = 3;

#ifdef USE_DGEMV
	for( int c = 0; c < 3; c++ )
	{
		double sv = scale / alphas[c]/alphas[c];
		dgemv( &trans, &n, &n, &sv, Mat, &n, vec_in+c, &incrxy, &one, vec_out+c, &incrxy );
	}
#else
	for( int vx = 0; vx < par_info.nv[surface_id]; vx++ )
	{
		int v1 = par_info.verts[surface_id][vx];

		for( int v2 = 0; v2 < nv; v2++ )
		{
			vec_out[3*v1+0] += Mat[v1*nv+v2] * vec_in[3*v2+0] * scale / alphas[0]/alphas[0];
			vec_out[3*v1+1] += Mat[v1*nv+v2] * vec_in[3*v2+1] * scale / alphas[0]/alphas[0];
			vec_out[3*v1+2] += Mat[v1*nv+v2] * vec_in[3*v2+2] * scale / alphas[0]/alphas[0];
		}
	}
#ifdef PARALLEL
	FullSyncVertices( vec_out, surface_id );
#endif
#endif

}

void SparseCartMatVecIncrScale( double *vec_out, double *vec_in, double *Mat, double scale, int nv, int *sparse_use, int nv_use, double *alphas, int surface_id )
{
	char trans = 'N';
	int n = nv;
	double one = 1.0;
	double zero = 0.0;
	int incrxy = 3;

	double *mcopy = (double *)malloc( sizeof(double) * nv_use*3 );
	double *vout  = (double *)malloc( sizeof(double) * par_info.nv[surface_id]*3 );

	for( int v2x = 0; v2x < nv_use; v2x++ )
	{
		int v2 = sparse_use[v2x];
		mcopy[v2x*3+0] = vec_in[3*v2+0];
		mcopy[v2x*3+1] = vec_in[3*v2+1];
		mcopy[v2x*3+2] = vec_in[3*v2+2];
	}

	for( int vx = 0; vx < par_info.nv[surface_id]; vx++ )
	{
		vout[3*vx+0] = 0;
		vout[3*vx+1] = 0;
		vout[3*vx+2] = 0;

		for( int v2x = 0; v2x < nv_use; v2x++ )
		{
			vout[3*vx+0] += Mat[vx*nv_use+v2x] * mcopy[3*v2x+0] * scale / (alphas[0]*alphas[0]);
			vout[3*vx+1] += Mat[vx*nv_use+v2x] * mcopy[3*v2x+1] * scale / (alphas[1]*alphas[1]);
			vout[3*vx+2] += Mat[vx*nv_use+v2x] * mcopy[3*v2x+2] * scale / (alphas[2]*alphas[2]);
		}
	}


	for( int vx = 0; vx < par_info.nv[surface_id]; vx++ )
	{
		int v1 = par_info.verts[surface_id][vx];
		
		vec_out[3*v1+0] += vout[3*vx+0];
		vec_out[3*v1+1] += vout[3*vx+1];
		vec_out[3*v1+2] += vout[3*vx+2];
	}


#ifdef PARALLEL
	FullSyncVertices( vec_out, surface_id );
#endif

	free( mcopy);
	free(vout);
}

void GenQMatVecIncrScale( double *vec_out, double *vec_in, SparseMatrix *Mat, double scale )
{
	double *mcopy = (double *)malloc( sizeof(double) * Mat->nsource );
	double *vout  = (double *)malloc( sizeof(double) * Mat->nneed);
	memset( vout, 0, sizeof(double)*Mat->nneed);
	
	Mat->compress_source_vector( vec_in, mcopy, 1 );
	Mat->mult( vout, mcopy, 1 );
	for( int v = 0; v < Mat->nneed; v++ )
		vout[v] *= scale;
	Mat->expand_target_vector( vec_out, vout,1 );
	
	free(mcopy);
	free(vout);
}

void PartialSyncGenQ( double *pp )
{
#ifdef PARALLEL
	if( par_info.nprocs > 1 )
	{
		printf("sync of genq not implemented.\n");
		MPI_Barrier(MPI_COMM_WORLD);
		MPI_Finalize();
		exit(1);
	}
#endif
}

void AltSparseCartMatVecIncrScale( double *vec_out, double *vec_in, SparseMatrix *Mat, double scale, double *alphas, int surface_id )
{
#ifdef PARALLEL
	FullSyncVertices( vec_in, surface_id );
//	PartialSyncVertices( vec_in, surface_id );
#endif

	double *mcopy = (double *)malloc( sizeof(double) * Mat->nsource*3 );

	double *vout  = (double *)malloc( sizeof(double) * Mat->nneed*3 );
	memset( vout, 0, sizeof(double)*3*Mat->nneed);

	Mat->compress_source_vector( vec_in,   mcopy, 3 );
	Mat->compress_source_vector( vec_in+1, mcopy+1, 3 );
	Mat->compress_source_vector( vec_in+2, mcopy+2, 3 );


	Mat->mult( vout, mcopy, 3 );
	Mat->mult( vout+1, mcopy+1, 3 );
	Mat->mult( vout+2, mcopy+2, 3 );
	
	
	for( int c = 0; c < 3; c++ )
	for( int v = 0; v < Mat->nneed; v++ )
		vout[v*3+c] *= scale / (alphas[c]*alphas[c]);

	Mat->expand_target_vector( vec_out, vout, 3 );
	Mat->expand_target_vector( vec_out+1, vout+1, 3 );
	Mat->expand_target_vector( vec_out+2, vout+2, 3 );
	


#ifdef PARALLEL
	FullSyncVertices( vec_out, surface_id );
//	PartialSyncVertices( vec_out );
#endif
	free(mcopy);
	free(vout);
}

// when a node finds that its complex has moved to another node, it transfers it to that node through the master communication node.


void SyncParInfo( void )
{
	/*
 *	Master node decides who gets which complexes.
 * 	*/

	/*
 *	Loop through our complexes and determine if any have moved to another node.
 * 	*/

	
}

void PartialGenVertices( double *total_vec, int do_sum, int surface_id )
{
#ifdef PARALLEL
//	Communicate only the necessary information to other processors. unblocked.

	// the data we are sending.

	int nv = par_info.nv[surface_id];

	double *send_buffers[par_info.nprocs];
	double *receive_buffers[par_info.nprocs];

	MPI_Request requests[par_info.nprocs];
	MPI_Request recv_requests[par_info.nprocs];
	MPI_Status status[par_info.nprocs];

	int **send_to = par_info.send_to[surface_id];
	int **get_from = par_info.get_from[surface_id];
	int *nsend_to = par_info.nsend_to[surface_id];
	int *nget_from = par_info.nget_from[surface_id];

	if( do_sum )
	{
		send_to = par_info.psum_send_to[surface_id];
		nsend_to = par_info.psum_nsend_to[surface_id];
		get_from = par_info.psum_get_from[surface_id];
		nget_from = par_info.psum_nget_from[surface_id];
	}


	for( int p = 0; p < par_info.nprocs; p++ )
	{
		send_buffers[p] = NULL;
		requests[p] = MPI_REQUEST_NULL;
		if( p == par_info.my_id ) 
			continue;
		send_buffers[p] = (double *)malloc( sizeof(double) * 3 * nsend_to[p] );
		
		double *send = send_buffers[p];
		for( int x = 0; x < nsend_to[p]; x++ )
		{
			int fetch = send_to[p][x];

			send[3*x+0] = total_vec[3*fetch+0];
			send[3*x+1] = total_vec[3*fetch+1];
			send[3*x+2] = total_vec[3*fetch+2];
		}

//		printf("ID %d ", surface_id ); printf("PROCESS %d placed ISEND to %d of %d units.\n", par_info.my_id, p, nsend_to[p]*3 );
		MPI_Isend( send, 3 * nsend_to[p], MPI_DOUBLE, p, par_info.my_id * par_info.nprocs + p, MPI_COMM_WORLD, requests+p ); 
	}
	
	for( int p = 0; p < par_info.nprocs; p++ )
	{
		receive_buffers[p] = NULL;
		recv_requests[p] = MPI_REQUEST_NULL;
		if( p == par_info.my_id ) 
			continue;
		receive_buffers[p] = (double *)malloc( sizeof(double) * 3 * nget_from[p] );
		
		double *receive = receive_buffers[p];

//		printf("ID %d ", surface_id );printf("PROCESS %d placed IRECV from %d of %d units.\n", par_info.my_id, p, nget_from[p]*3 );
		MPI_Irecv( receive, 3 * nget_from[p], MPI_DOUBLE, p, p * par_info.nprocs + par_info.my_id, MPI_COMM_WORLD, recv_requests+p ); 

	}
		
	MPI_Barrier(MPI_COMM_WORLD);

	int not_done = 1;
	while( not_done )
	{
		int p=-1;
//		printf("ID %d ", surface_id );printf("Process %d is waiting.\n", par_info.my_id );
		MPI_Waitany( par_info.nprocs, recv_requests, &p, status );

		if( p >=0 )
		{
			double *receive = receive_buffers[p];

//			if( do_sum )
//			printf("Receiving %d from %d.\n", nget_from[p], p );

			for( int x = 0; x < nget_from[p]; x++ )
			{
				int fetch = get_from[p][x];

				if( do_sum )
				{	
					total_vec[3*fetch+0] += receive[3*x+0];
					total_vec[3*fetch+1] += receive[3*x+1];
					total_vec[3*fetch+2] += receive[3*x+2];
				}
				else
				{
					total_vec[3*fetch+0] = receive[3*x+0];
					total_vec[3*fetch+1] = receive[3*x+1];
					total_vec[3*fetch+2] = receive[3*x+2];
				}
			}

			free(receive);
			recv_requests[p] = MPI_REQUEST_NULL;
		}
		else
			not_done = 0;		
	}
/*	
	while( not_done )
	{
		int p=-1;
//		printf("Process %d is waiting.\n", par_info.my_id );
		MPI_Waitany( par_info.nprocs, recv_requests, &p, status );

		if( p >=0 )
		{
			double *receive = receive_buffers[p];

			for( int x = 0; x < par_info.nget_from[p]; x++ )
			{
				int fetch = par_info.get_from[p][x];
	
				total_vec[3*fetch+0]=receive[3*x+0];
				total_vec[3*fetch+1]=receive[3*x+1];
				total_vec[3*fetch+2]=receive[3*x+2];
			}

			free(receive);
			recv_requests[p] = MPI_REQUEST_NULL;
		}
		else
			not_done = 0;		
	}
*/	
	not_done = 1;
	while( not_done )
	{
		int p=-1;
		MPI_Waitany( par_info.nprocs, requests, &p, status );

		if( p >=0 )
		{
			free(send_buffers[p]);
		}
		else
			not_done = 0;		
	}
/*
	for( int p = 0; p < par_info.nprocs; p++ )
	{
		if( p == par_info.my_id ) 
			continue;
	}
*/
#endif
}

void PartialSumVertices( double *total_vec, int surface_id )
{
	PartialGenVertices( total_vec, 1 /* sum */, surface_id  );
}

void PartialSyncVertices( double *total_vec, int surface_id )
{
	PartialGenVertices( total_vec, 0 /* sum */, surface_id );
}

void surface::lipidSync(void)
{
#ifdef PARALLEL
	int *list = (int *)malloc( sizeof(int) * nt );
	double *fl = (double *)malloc( sizeof(double) * nt * 2 * bilayerComp.nlipidTypes);

	for( int p = 0; p < par_info.nprocs; p++ )
	{
		if( p == par_info.my_id )
		{
			for( int fx = 0; fx < par_info.nf[surface_id]; fx++ )
			{
				int f = par_info.faces[surface_id][fx];
				int t;
				if( f < nf_faces )
					t = theFormulas[f*nf_g_q_p].tri;
				else
					t = theIrregularFormulas[(f-nf_faces)*nf_irr_pts].tri;

				list[fx] = t;

				for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
				{
					fl[fx*bilayerComp.nlipidTypes*2+lx] = theTriangles[t].composition.innerLeaflet[lx];
					fl[fx*bilayerComp.nlipidTypes*2+bilayerComp.nlipidTypes+lx] = theTriangles[t].composition.outerLeaflet[lx];
				}
			}
		}

		int nvals = par_info.nf[surface_id];
		
		MPI_Bcast( &nvals, 1, MPI_INT, p, MPI_COMM_WORLD );		
		MPI_Bcast( list, nvals, MPI_INT, p, MPI_COMM_WORLD ); 
		MPI_Bcast( fl, nvals*bilayerComp.nlipidTypes*2, MPI_DOUBLE, p, MPI_COMM_WORLD ); 
	
		for( int tx = 0; tx < nvals; tx++ )
		{
			int t = list[tx];

			for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
			{
				theTriangles[t].composition.innerLeaflet[lx]	=	fl[tx*bilayerComp.nlipidTypes*2+lx];
				theTriangles[t].composition.outerLeaflet[lx]	=	fl[tx*bilayerComp.nlipidTypes*2+bilayerComp.nlipidTypes+lx];
			}
		}
	}
	free(list);
	free(fl);
#endif
}

void surface::lipidBroadcast(void)
{
#ifdef PARALLEL
	double *fl = (double *)malloc( sizeof(double) * nt * 2 * bilayerComp.nlipidTypes );

	if( par_info.my_id == BASE_TASK )
	{
		for( int t = 0; t < nt; t++ )
		{
			for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
			{	
				fl[t*2*bilayerComp.nlipidTypes+lx]                         = theTriangles[t].composition.innerLeaflet[lx];
				fl[t*2*bilayerComp.nlipidTypes+bilayerComp.nlipidTypes+lx] = theTriangles[t].composition.outerLeaflet[lx];
			}
		}	
	}

	ParallelBroadcast( fl, nt*2*bilayerComp.nlipidTypes );

	for( int t = 0; t < nt; t++ )
	{
		for( int lx = 0; lx < bilayerComp.nlipidTypes; lx++ )
		{	
			theTriangles[t].composition.innerLeaflet[lx]=fl[t*2*bilayerComp.nlipidTypes+lx];
			theTriangles[t].composition.outerLeaflet[lx]=fl[t*2*bilayerComp.nlipidTypes+bilayerComp.nlipidTypes+lx];
		}
	}
	
	free(fl);
#endif
}

