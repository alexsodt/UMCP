#include "interp.h"
#include "pcomplex.h"
#include <math.h>
#include <string.h>
#include "units.h"
#include "p_p.h"
#include "parallel.h"
#define FD_V_MESH
//#define FD_V_P
//#define FD_MESH_QDOT
#define FD_T_MESH
//#define FD_T_P

#define DO_PP


static double FD_TOL = 1e-3;
#if 0 // pre sim
void surface::debug_dynamics( double *r, force_set *theForceSet,  double *Minv, pcomplex **allComplexes, int ncomplex )
{
	double pmag = 10;
	double *mesh_p = (double *)malloc( sizeof(double) * 3 * nv );

	for( int v = 0; v < nv; v++ )
	{
		mesh_p[3*v+0] = pmag * 2 * ((double)rand()/(double)RAND_MAX-0.5);	
		mesh_p[3*v+1] = pmag * 2 * ((double)rand()/(double)RAND_MAX-0.5);	
		mesh_p[3*v+2] = pmag * 2 * ((double)rand()/(double)RAND_MAX-0.5);	
	}	

	double pmagp = 10;

	for( int c = 0; c < ncomplex; c++ )
	{
		for( int s = 0; s < allComplexes[c]->nattach; s++ )
		{
			allComplexes[c]->p[2*s+0] = pmagp * 2 * ((double)rand()/(double)RAND_MAX-0.5);
			allComplexes[c]->p[2*s+1] = pmagp * 2 * ((double)rand()/(double)RAND_MAX-0.5);
		}
	}

	int nattach = 0;
	for( int c = 0; c < ncomplex; c++ )
		nattach += allComplexes[c]->nattach;

	double *mesh_g_V = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	double *mesh_g_T = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	double *mesh_qdot = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	double *mesh_qdot0 = (double *)malloc( sizeof(double) * 3 * (nv+1) );
	double *mesh_qdot_temp = (double *)malloc( sizeof(double) * 3 * nv );
	double *pgrad_V = (double *)malloc( sizeof(double) * 2*nattach );
	double *pgrad_T = (double *)malloc( sizeof(double) * 2*nattach );

	memset( mesh_g_V, 0, sizeof(double) * 3 * (nv+1) );
	memset( mesh_g_T, 0, sizeof(double) * 3 * (nv+1) );
	
	memset( pgrad_V, 0, sizeof(double) * 2 * nattach );
	memset( pgrad_T, 0, sizeof(double) * 2 * nattach );

	grad( r, mesh_g_V );
	
	for( int i = 0; i < nv; i++ )
	{	
		mesh_qdot0[3*i+0] = 0;
		mesh_qdot0[3*i+1] = 0;
		mesh_qdot0[3*i+2] = 0;
	
		for( int k = 0; k < nv; k++ )
		{
			mesh_qdot0[3*i+0] += Minv[i*nv+k] * mesh_p[3*k+0];		
			mesh_qdot0[3*i+1] += Minv[i*nv+k] * mesh_p[3*k+1];		
			mesh_qdot0[3*i+2] += Minv[i*nv+k] * mesh_p[3*k+2];		
		}
	}

	memcpy( mesh_qdot, mesh_qdot0, sizeof(double) * 3 * nv );

	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->compute_qdot( this, r, mesh_qdot0, mesh_qdot );


	double *mesh_der_qdot = (double *)malloc( sizeof(double) * 3 * nv );
	memset( mesh_der_qdot, 0, sizeof(double) * 3 * nv );

	printf("BROKEN HERE.\n");
	exit(1);
//	for( int c = 0; c < ncomplex; c++ )	
//		allComplexes[c]->getMeshQxdot( this, r, Minv, mesh_p, mesh_qdot, mesh_qdot0, mesh_der_qdot ); 

	int sglobal = 0;
	
	for( int c = 0; c < ncomplex; c++ )
		memset( allComplexes[c]->save_grad, 0, sizeof(double) * 2 * allComplexes[c]->nattach );

#ifdef DO_PP
//	PP_G( this, r, allComplexes, ncomplex, mesh_g_V );
	ParallelSyncComplexes( allComplexes, ncomplex );
	Boxed_PP_G( this, r, allComplexes, ncomplex, mesh_g_V );		
	ParallelSyncComplexes( allComplexes, ncomplex );

	for( int c = 0; c < ncomplex; c++ )
	{
		for( int s = 0; s < allComplexes[c]->nattach;  s++ )
		{
			pgrad_V[sglobal*2+0] += allComplexes[c]->save_grad[2*s+0];
			pgrad_V[sglobal*2+1] += allComplexes[c]->save_grad[2*s+1];

			sglobal += 1;
		}
	}
#endif

	sglobal = 0;
	for( int c = 0; c < ncomplex; c++ )	
	{
		allComplexes[c]->debug_dynamics( this, r,
			mesh_qdot, mesh_qdot0,
			mesh_g_V,
			mesh_g_T,
			pgrad_V+sglobal*2,
			pgrad_T+sglobal*2 );
		sglobal += allComplexes[c]->nattach;
	}
	double eps_mesh = 0.000001;
	double eps_uv = 1e-7;

#ifdef FD_V_MESH
	//check with finite difference.
	printf("FINITE DIFFERENCE CHECK of dV/dq mesh.\n");	 	


	for( int v = 0; v < nv; v++ )
	{
		double fdg_c[3] = { 0,0,0};
		for( int c = 0; c < 3; c++ )
		{
			double epm[2] = {0,0};

			double meshd[2] = { -eps_mesh, eps_mesh };
			for( int dx = 0; dx < 2; dx++ )
			{
				r[3*v+c] += meshd[dx];
				epm[dx] = energy( r, NULL );
				for( int cx = 0; cx < ncomplex; cx++ )
				{	
					epm[dx] += allComplexes[cx]->V( this, r );
					epm[dx] += allComplexes[cx]->AttachV( this, r );
				}
#ifdef DO_PP
				double pp_v = Boxed_PP_V( this, r, allComplexes, ncomplex );
				ParallelSum( &pp_v, 1 );
				epm[dx] += pp_v;
#endif
				r[3*v+c] -= meshd[dx];
			}

			double fdg = (epm[1]-epm[0]) / ( meshd[1] - meshd[0]);

			fdg_c[c] = fdg;
		}

		double rel_err = sqrt( 
				( fdg_c[0] - mesh_g_V[3*v+0])*(fdg_c[0] - mesh_g_V[3*v+0]) +
				( fdg_c[1] - mesh_g_V[3*v+1])*(fdg_c[1] - mesh_g_V[3*v+1]) +
				( fdg_c[2] - mesh_g_V[3*v+2])*(fdg_c[2] - mesh_g_V[3*v+2]) ) / sqrt( fdg_c[0]*fdg_c[0]+fdg_c[1]*fdg_c[1]+fdg_c[2]*fdg_c[2] +
												     mesh_g_V[3*v+0]*mesh_g_V[3*v+0]+
												     mesh_g_V[3*v+1]*mesh_g_V[3*v+1]+
												     mesh_g_V[3*v+2]*mesh_g_V[3*v+2] + 1e-10 );
		if( rel_err < FD_TOL )
			printf(" OK ");
		else
			printf(" CHECK ");
		printf("Finite difference: %le %le %le G: %le %le %le diff %le %le %le\n",
			fdg_c[0], fdg_c[1], fdg_c[2],  mesh_g_V[3*v+0], mesh_g_V[3*v+1], mesh_g_V[3*v+2],
			fabs( fdg_c[0] - mesh_g_V[3*v+0]),	
			fabs( fdg_c[1] - mesh_g_V[3*v+1]),	
			fabs( fdg_c[2] - mesh_g_V[3*v+2]) );	
	}
#endif
	
#ifdef FD_V_P
	printf("FINITE DIFFERENCE CHECK of dV/dq particle.\n");	 	

	sglobal = 0;

	for( int c = 0; c < ncomplex; c++ )
	{
		printf("Complex %d\n",c ); //, %s\n", c );
		for( int s = 0; s < allComplexes[c]->nattach; s++, sglobal++ )
		{
			double fdg_uv[2] = { 0,0};
			for( int uv = 0; uv < 2; uv++ )
			{
				double epm[2] = {0,0};

				double meshd[2] = { -eps_uv, eps_uv };
				for( int dx = 0; dx < 2; dx++ )
				{
					allComplexes[c]->puv[2*s+uv] += meshd[dx];
					epm[dx] = energy( r, NULL );
					for( int c = 0; c < ncomplex; c++ )
					{
						epm[dx] += allComplexes[c]->V( this, r );
						epm[dx] += allComplexes[c]->AttachV( this, r );
					}
#ifdef DO_PP
					double pp_v = Boxed_PP_V( this, r, allComplexes, ncomplex );
					ParallelSum( &pp_v, 1 );
					epm[dx] += pp_v;
#endif
					allComplexes[c]->puv[2*s+uv] -= meshd[dx];
				}
			
				double fdg = (epm[1]-epm[0]) / ( meshd[1] - meshd[0]);

				fdg_uv[uv] = fdg;
			}
			
			double rel_err = sqrt( 
				( fdg_uv[0] - pgrad_V[2*sglobal+0])*(fdg_uv[0] - pgrad_V[2*sglobal+0]) +
				( fdg_uv[1] - pgrad_V[2*sglobal+1])*(fdg_uv[1] - pgrad_V[2*sglobal+1]) 
			) / sqrt( fdg_uv[0]*fdg_uv[0]+fdg_uv[1]*fdg_uv[1]+
												     pgrad_V[2*sglobal+0]*pgrad_V[2*sglobal+0]+
												     pgrad_V[2*sglobal+1]*pgrad_V[2*sglobal+1]
												     + 1e-10 );
		if( rel_err < FD_TOL )
			printf(" OK ");
		else
			printf(" CHECK ");
		
			printf("Finite difference: %le %le G: %le %le diff %le %le\n",
				fdg_uv[0], fdg_uv[1],  
				pgrad_V[2*sglobal+0], pgrad_V[2*sglobal+1], 
				fabs( fdg_uv[0] - pgrad_V[2*sglobal+0]),	
				fabs( fdg_uv[1] - pgrad_V[2*sglobal+1]) );	
		}
	}
#endif

#ifdef FD_MESH_QDOT
	//check with finite difference.
	printf("FINITE DIFFERENCE CHECK of d \ dot{q_x} / d mesh.\n");	 	

	for( int v = 0; v < nv; v++ )
	{
		double fdg_c[3] = { 0,0,0};
		for( int c = 0; c < 3; c++ )
		{
			double epm[2] = {0,0};

			double meshd[2] = { -eps_mesh, eps_mesh };
			for( int dx = 0; dx < 2; dx++ )
			{
				r[3*v+c] += meshd[dx];
				memset( mesh_qdot_temp, 0, sizeof(double ) * 3 * nv );
				memcpy( mesh_qdot, mesh_qdot0, sizeof(double) * 3 * nv );

				for( int c = 0; c < ncomplex; c++ )
					allComplexes[c]->compute_qdot( this, r, mesh_qdot0, mesh_qdot_temp );

				for( int c1 = 0; c1 < nv; c1++ )
				for( int c2 = 0; c2 < nv; c2++ )
				{
					mesh_qdot[3*c1+0] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+0];
					mesh_qdot[3*c1+1] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+1];
					mesh_qdot[3*c1+2] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+2];
				}

				epm[dx] = mesh_qdot[3*v+0];

				r[3*v+c] -= meshd[dx];
			}

			double fdg = (epm[1]-epm[0]) / ( meshd[1] - meshd[0]);

			fdg_c[c] = fdg;
		}
	
		printf("Finite difference: %le %le %le G: %le %le %le diff %le %le %le\n",
			fdg_c[0], fdg_c[1], fdg_c[2],  mesh_der_qdot[3*v+0], mesh_der_qdot[3*v+1], mesh_der_qdot[3*v+2],
			fabs( fdg_c[0] - mesh_der_qdot[3*v+0]),	
			fabs( fdg_c[1] - mesh_der_qdot[3*v+1]),	
			fabs( fdg_c[2] - mesh_der_qdot[3*v+2]) );	
	}
#endif
	
#ifdef FD_T_MESH
	//check with finite difference.
	printf("FINITE DIFFERENCE CHECK of dT/dq mesh.\n");	 	

	for( int v = 0; v < nv; v++ )
	{
		double fdg_c[3] = { 0,0,0};
		for( int c = 0; c < 3; c++ )
		{
			double epm[2] = {0,0};

			double meshd[2] = { -eps_mesh, eps_mesh };
			for( int dx = 0; dx < 2; dx++ )
			{
				r[3*v+c] += meshd[dx];
				memset( mesh_qdot_temp, 0, sizeof(double ) * 3 * nv );
				memcpy( mesh_qdot, mesh_qdot0, sizeof(double) * 3 * nv );

				for( int cx = 0; cx < ncomplex; cx++ )
					allComplexes[cx]->compute_qdot( this, r, mesh_qdot0, mesh_qdot_temp );

				for( int c1 = 0; c1 < nv; c1++ )
				for( int c2 = 0; c2 < nv; c2++ )
				{
					mesh_qdot[3*c1+0] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+0];
					mesh_qdot[3*c1+1] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+1];
					mesh_qdot[3*c1+2] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+2];
				}

				
	
				epm[dx] = evaluate_T( mesh_qdot, mesh_p,NULL,NULL,r+3*nv );
				for( int cx = 0; cx < ncomplex; cx++ )
					epm[dx] += allComplexes[cx]->T(this,r);
				r[3*v+c] -= meshd[dx];
			}

			double fdg = (epm[1]-epm[0]) / ( meshd[1] - meshd[0]);

			fdg_c[c] = fdg;
		}
		
		double rel_err = sqrt( 
				( fdg_c[0] - mesh_g_T[3*v+0])*(fdg_c[0] - mesh_g_T[3*v+0]) +
				( fdg_c[1] - mesh_g_T[3*v+1])*(fdg_c[1] - mesh_g_T[3*v+1]) +
				( fdg_c[2] - mesh_g_T[3*v+2])*(fdg_c[2] - mesh_g_T[3*v+2]) ) / sqrt( fdg_c[0]*fdg_c[0]+fdg_c[1]*fdg_c[1]+fdg_c[2]*fdg_c[2] +
												     mesh_g_T[3*v+0]*mesh_g_T[3*v+0]+
												     mesh_g_T[3*v+1]*mesh_g_T[3*v+1]+
												     mesh_g_T[3*v+2]*mesh_g_T[3*v+2] + 1e-10 );
		if( rel_err < FD_TOL )
			printf(" OK ");
		else
			printf(" CHECK ");
	
		printf("Finite difference: %le %le %le G: %le %le %le diff %le %le %le\n",
			fdg_c[0], fdg_c[1], fdg_c[2],  mesh_g_T[3*v+0], mesh_g_T[3*v+1], mesh_g_T[3*v+2],
			fabs( fdg_c[0] - mesh_g_T[3*v+0]),	
			fabs( fdg_c[1] - mesh_g_T[3*v+1]),	
			fabs( fdg_c[2] - mesh_g_T[3*v+2]) );	
	}
#endif

#ifdef FD_T_P	
	printf("FINITE DIFFERENCE CHECK of dT/dq particle.\n");	 	

	sglobal = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		printf("Complex %d\n",c ); //, %s\n", c );
		for( int s = 0; s < allComplexes[c]->nattach; s++, sglobal++ )
		{
			double fdg_uv[2] = { 0,0};
			for( int uv = 0; uv < 2; uv++ )
			{
				double epm[2] = {0,0};

				double meshd[2] = { -eps_uv, eps_uv  };
				for( int dx = 0; dx < 2; dx++ )
				{
					allComplexes[c]->puv[2*s+uv] += meshd[dx];
					allComplexes[c]->grad_puv[2*s+uv] += meshd[dx];
					
					memset( mesh_qdot_temp, 0, sizeof(double ) * 3 * nv );
					memcpy( mesh_qdot, mesh_qdot0, sizeof(double) * 3 * nv );
					for( int c = 0; c < ncomplex; c++ )
						allComplexes[c]->compute_qdot( this, r,  mesh_qdot0, mesh_qdot_temp );
				
					for( int c1 = 0; c1 < nv; c1++ )
					for( int c2 = 0; c2 < nv; c2++ )
					{
						mesh_qdot[3*c1+0] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+0];
						mesh_qdot[3*c1+1] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+1];
						mesh_qdot[3*c1+2] += Minv[c1*nv+c2] * mesh_qdot_temp[3*c2+2];
					}

					epm[dx] = evaluate_T(  mesh_qdot, mesh_p );
					for( int c = 0; c < ncomplex; c++ )
						epm[dx] += allComplexes[c]->T();
					allComplexes[c]->puv[2*s+uv] -= meshd[dx];
					allComplexes[c]->grad_puv[2*s+uv] -= meshd[dx];
				}
			
				double fdg = (epm[1]-epm[0]) / ( meshd[1] - meshd[0]);

				fdg_uv[uv] = fdg;
			}
		
			double rel_err = sqrt( 
				( fdg_uv[0] - pgrad_T[2*sglobal+0])*(fdg_uv[0] - pgrad_T[2*sglobal+0]) +
				( fdg_uv[1] - pgrad_T[2*sglobal+1])*(fdg_uv[1] - pgrad_T[2*sglobal+1]) 
			) / sqrt( fdg_uv[0]*fdg_uv[0]+fdg_uv[1]*fdg_uv[1]+
												     pgrad_T[2*sglobal+0]*pgrad_T[2*sglobal+0]+
												     pgrad_T[2*sglobal+1]*pgrad_T[2*sglobal+1]
												     + 1e-10 );
		if( rel_err < FD_TOL )
			printf(" OK ");
		else
			printf(" CHECK ");
		
			printf("Finite difference: %le %le G: %le %le diff %le %le\n",
				fdg_uv[0], fdg_uv[1],  
				pgrad_T[2*sglobal+0], pgrad_T[2*sglobal+1], 
				fabs( fdg_uv[0] - pgrad_T[2*sglobal+0]),	
				fabs( fdg_uv[1] - pgrad_T[2*sglobal+1]) );	
		}
	}
#endif
	
}
#endif

void surface::timestep_analysis( double *r, force_set *theForceSet,  double *Minv, pcomplex **allComplexes, int ncomplex, double approx_timestep )
{
#if 0   //pre-simulation, change me
	// for now do the diagonals. could diagonalize the whole system...

	int nscan = 5;
	int n_boltzmann = 3;

	double min_TS = 1;
	char dof[256];
	sprintf(dof, "UNKNOWN");

#if 1 
	printf("Analyzing mesh.\n");
	for( int q_mesh = 0; q_mesh < 3 * nv; q_mesh++ )
	{
		int iq = q_mesh/3;

		double vals[nscan];

		double eff_mass = 1.0 / Minv[iq*nv+iq];

		// n kT == m/2 qdot * qdot
		double qdot = sqrt( n_boltzmann * 0.592 / (eff_mass/2) );
		// that's a "reasonable" qdot.
		double eps_q = qdot * approx_timestep * AKMA_TIME;
		
//		printf("eps_q: %le\n", eps_q);

		double pe = 0;
		for( int dq = 0; dq < nscan; dq++ )
		{
			r[q_mesh] += dq * eps_q;
				
			vals[dq] = energy( r, NULL );

			for( int c = 0; c < ncomplex; c++ )
			{	
				vals[dq] += allComplexes[c]->V( this, r );
				vals[dq] += allComplexes[c]->AttachV( this, r );
			}

//			printf(" %le", vals[dq] );

			pe = vals[dq];


			r[q_mesh] -= dq * eps_q;				
		}

		double Gs[nscan-1];
		for(int dq = 0; dq < nscan-1; dq++ )
		{
			Gs[dq] = (vals[dq+1]-vals[dq])/eps_q;
//			printf(" %le", Gs[dq] );
		}
//		printf("\n");

		double approx_k = 0;
		for( int dq = 0; dq < nscan-2; dq++ )
			approx_k += (Gs[dq+1]-Gs[dq])/eps_q;
		approx_k /= (nscan-2);


		// V(q) ~ k q^2

		int steps_per_period = 10;

		double freq = sqrt( approx_k / eff_mass );
	
		double period = 2 * M_PI / freq;

		double ts = period/10 / AKMA_TIME;

//		printf("approx_k %le D.O.F. TS %le\n", approx_k, ts / AKMA_TIME ); 

		if( ts < min_TS )
		{
			min_TS = ts;
			 
			const char *xyz = "xyz";
			sprintf(dof, "MESH %d %c", iq, xyz[q_mesh%3] ); 
		}
	}

	printf("MESH info %s dt %le\n", dof, min_TS );
#endif
	printf("Analyzing complexes.\n");

	for( int pc = 0; pc < ncomplex; pc++ )
	{
		for( int s = 0; s < allComplexes[pc]->nattach; s++ )
		{
			double P_uv[4];
			
			int f = allComplexes[pc]->fs[s];
			double u = allComplexes[pc]->puv[2*s+0];
			double v = allComplexes[pc]->puv[2*s+1];

			fetchPuv( f, u, v, P_uv, r );
			
			P_uv[0] *= allComplexes[pc]->mass[s];
			P_uv[1] *= allComplexes[pc]->mass[s];
			P_uv[2] *= allComplexes[pc]->mass[s];
			P_uv[3] *= allComplexes[pc]->mass[s];
	
			double Pdet = P_uv[0] * P_uv[3] - P_uv[1] * P_uv[2];
			double Pinv0[4] = { P_uv[3]/Pdet, -P_uv[2]/Pdet, -P_uv[1]/Pdet, P_uv[0]/Pdet };

			for( int iu = 0; iu < 2; iu++ )
			{
				double vals[nscan];
				double eff_mass = 1.0 / Pinv0[iu*2+iu];
		
				// n kT == m/2 qdot * qdot
				double qdot = sqrt( n_boltzmann * 0.592 / (eff_mass/2) );
				// that's a "reasonable" qdot.
				double eps_q = qdot * approx_timestep * AKMA_TIME;
				
		//		printf("eps_q: %le\n", eps_q);
		
				double pe = 0;
				for( int dq = 0; dq < nscan; dq++ )
				{
					allComplexes[pc]->puv[2*s+iu] += dq * eps_q;
					
					vals[dq] = energy( r, NULL );
		
					for( int c = 0; c < ncomplex; c++ )
					{	
						vals[dq] += allComplexes[c]->V( this, r );
						vals[dq] += allComplexes[c]->AttachV( this, r );
					}
				
#ifdef DO_PP
					double pp_v = Boxed_PP_V( this, r, allComplexes, ncomplex );
					ParallelSum( &pp_v, 1 );
					vals[dq] += pp_v;
#endif
					pe = vals[dq];
		
					allComplexes[pc]->puv[2*s+iu] -= dq * eps_q;
					
				}
		
				double Gs[nscan-1];
				for(int dq = 0; dq < nscan-1; dq++ )
				{
					Gs[dq] = (vals[dq+1]-vals[dq])/eps_q;
		//			printf(" %le", Gs[dq] );
				}
		//		printf("\n");
		
				double approx_k = 0;
				for( int dq = 0; dq < nscan-2; dq++ )
					approx_k += (Gs[dq+1]-Gs[dq])/eps_q;
				approx_k /= (nscan-2);
		
		
				// V(q) ~ k q^2
		
				int steps_per_period = 10;
		
				double freq = sqrt( approx_k / eff_mass );
			
				double period = 2 * M_PI / freq;
		
				double ts = period/10 / AKMA_TIME;
		
		//		printf("approx_k %le D.O.F. TS %le\n", approx_k, ts / AKMA_TIME ); 
		
				if( ts < min_TS )
				{
					min_TS = ts;
					 
					const char *xyz = "uv";
					sprintf(dof, "COMPLEX %d sub-particle %d %c", pc, s, xyz[iu] ); 
				}
			}
		}
	}

	printf("Timestep %le limited by %s.\n", min_TS, dof );
#endif
}































