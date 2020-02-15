#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mutil.h"
#include "interp.h"
#include "parallel.h"
const double F_CUTOFF = 0.5;

int rand_dbg( void )
{
	return RAND_MAX/2;
}

void surface::stashf( void )
{
	for( int t = 0; t < nt; t++ )
		theTriangles[t].f_lipids_stashed = theTriangles[t].f_lipids; 
}

void surface::set_g0_from_f(int f)
{
	if( f < nf_faces )
	{
		int tri = theFormulas[f*nf_g_q_p].tri;

		double c0_i  = 0;
		double c0_o  = 0;
		double atot_i = 0;
		double atot_o = 0;

		for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
		{
			c0_o += bilayerComp.c0[x] * theTriangles[tri].composition.outerLeaflet[x];
			atot_o += theTriangles[tri].composition.outerLeaflet[x];

			// opposite sign for inner leaflet, taken into account below.
			c0_i += bilayerComp.c0[x] * theTriangles[tri].composition.innerLeaflet[x];
			atot_i += theTriangles[tri].composition.innerLeaflet[x];
		}
		
		c0_i /= atot_i;
		c0_o /= atot_o;

		double atot = (atot_o+atot_i)/2;

		theTriangles[tri].composition.A_inst = atot;

		for( int q = 0; q < nf_g_q_p; q++ )
		{	// needs to come out with the change in g0 being proportional to "innerLeaflet[xBab]" (with the same proportionality!)
			theFormulas[f*nf_g_q_p+q].g0 = theFormulas[f*nf_g_q_p].g0_base * (atot / theTriangles[tri].composition.A0); 
			theFormulas[f*nf_g_q_p+q].c0 = (c0_o - c0_i)/2;
		}
	}
	else
	{
		int df = f-nf_faces;
		int tri = theIrregularFormulas[df*nf_irr_pts].tri;
		
		double c0_i  = 0;
		double c0_o  = 0;
		double atot_i = 0;
		double atot_o = 0;

		for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
		{
			c0_o += bilayerComp.c0[x] *  theTriangles[tri].composition.outerLeaflet[x];
			atot_o += theTriangles[tri].composition.outerLeaflet[x];

			// opposite sign for inner leaflet, taken into account below
			c0_i -= bilayerComp.c0[x] * theTriangles[tri].composition.innerLeaflet[x];
			atot_i += theTriangles[tri].composition.innerLeaflet[x];
		}
		
		double atot = (atot_o+atot_i)/2;
		theTriangles[tri].composition.A_inst = atot;

		c0_o /= atot_o;
		c0_i /= atot_i;
		
//		if( tri == 0 ) printf("c0: %lf A: %lf A0: %lf\n", c0, atot, theTriangles[tri].composition.A0 );

		for( int q = 0; q < nf_irr_pts; q++ )
		{
			theIrregularFormulas[df*nf_irr_pts+q].g0 = theIrregularFormulas[df*nf_irr_pts+q].g0_base * (atot / theTriangles[tri].composition.A0); 
			theIrregularFormulas[df*nf_irr_pts+q].c0 = c0; 
		}
	}
/*
	if( f < nf_faces )
	{
		int tri = theFormulas[f*nf_g_q_p].tri;

		for( int q = 0; q < nf_g_q_p; q++ )
			theFormulas[f*nf_g_q_p+q].g0 = theFormulas[f*nf_g_q_p].g0_base * theTriangles[tri].f_lipids; 
	}
	else
	{
		int df = f-nf_faces;
		int tri = theIrregularFormulas[df*nf_irr_pts].tri;

		for( int q = 0; q < nf_irr_pts; q++ )
			theIrregularFormulas[df*nf_irr_pts+q].g0 = theIrregularFormulas[df*nf_irr_pts+q].g0_base * theTriangles[tri].f_lipids; 
	}
*/
}

void surface::set_g0_from_f(void)
{
	for( int f = 0; f < nt; f++ )
		set_g0_from_f(f);
}

void surface::unstashf( void )
{
	for( int t = 0; t < nt; t++ )
		theTriangles[t].f_lipids = theTriangles[t].f_lipids_stashed; 
	
	set_g0_from_f();	
}



void surface::local_lipidMCMove( double *r, pcomplex **allComplexes, int ncomplex, double dt, double beta, double DC, double PBC[3][3], int swap_only, int *tracer)
{
	// move lipids between adjacent faces using montecarlo moves.
	// for now ignore the complexes.
	
//	double E0 = energy(r,NULL);
//	stashf();
	// let's set this based on dt.
	double nL_move_average = 1./10.0; // times the area.

	FullSyncVertices(r, surface_id);

	int nacc = 0;
	int nrej = 0;

	double lipids_moved = 0;

	if( par_info.my_id == BASE_TASK )
	{
		//for( int fx = 0; fx < par_info.nf; fx++ )
		//{
		//	int f = par_info.faces[fx];
		for( int f = 0; f < nt; f++ )
		{
			int t1;
	
			if( f < nf_faces ) 
				t1 = theFormulas[f*nf_g_q_p].tri;
			else
				t1 = theIrregularFormulas[(f-nf_faces)*nf_irr_pts].tri;
			
			for( int bt = 0; bt < 3; bt++ )
			{
				int t2 = theTriangles[t1].border_tri[bt];
	
				if( t2 < t1 ) continue;
				double E0 = 0;
	
				int f1 = theTriangles[t1].f;
				int f2 = theTriangles[t2].f;
	
				for( int leaf = 0; leaf < 2; leaf++ )
				{
	
					double *comp1 = theTriangles[t1].composition.innerLeaflet;
					double *comp2 = theTriangles[t2].composition.innerLeaflet;
		
					if( leaf == 1 )
					{
						comp1 = theTriangles[t1].composition.outerLeaflet;
						comp2 = theTriangles[t2].composition.outerLeaflet;
					}
	
					int from_1 = rand() % bilayerComp.nlipidTypes;
	
					for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
					{
						double nL1 = comp1[from_1];
						int n_here = lround(nL1/bilayerComp.APL[from_1]);

						int l_select_1 = -1;

						if( nL1 > 0 ) l_select_1 = rand() % n_here;

						// how many should move? 	
						// compute distance between triangle centers.
						double r1[3];
						double r2[3],nrmjunk[3];
						evaluateRNRM( f1, 1.0/3.0, 1.0/3.0, r1, nrmjunk, r );
						evaluateRNRM( f2, 1.0/3.0, 1.0/3.0, r2, nrmjunk, r );
						double dr[3] = { r2[0]-r1[0],r2[1]-r1[1],r2[2]-r1[2] };
						double put[3]={0,0,0};
						MinImage3D( dr, PBC_vec, put, r+3*nv );
						double del_r = normalize(dr);

						double dist = sqrt( DC * dt );
						double n_move = dist / del_r;

						while( n_move > 0 )
						{
							double pr = rand() / (double)RAND_MAX;

							if( pr > n_move ) break;

							int from_2 = from_1;
				
							if( swap_only )
								from_2 = rand() % bilayerComp.nlipidTypes; 

							E0 =  faceEnergy( f1, r, NULL, 0 );
							E0 += faceEnergy( f2, r, NULL, 0 );
		
							// move lipid material between triangles.
			
							double nL2 = comp2[from_2];
				
	//						double dL = 2*((double)rand() / (double)RAND_MAX - 0.5 ) * (theTriangles[t1].composition.A0) * nL_move_average;
							double dL = bilayerComp.APL[from_1] * (rand() %2 == 0 ? 1 : -1);
							double dL2 = -dL;
	
							if( swap_only )
								dL2 = bilayerComp.APL[from_2] * (dL/bilayerComp.APL[from_1]);
			
	
							int l_select_2 = -1;
							if( nL2 > 0 ) l_select_2 = rand() % lround(nL2/bilayerComp.APL[from_2]);
	
							int move_tracer_1 = (theTriangles[t1].composition.tracer && from_1 == 0 && l_select_1 == 0 ); 
							int move_tracer_2 = (theTriangles[t2].composition.tracer && from_2 == 0 && l_select_2 == 0 ); 
		
							comp1[from_1] -= dL;
							comp2[from_1] += dL;
							if( swap_only )
							{
								comp1[from_2] += dL2;
								comp2[from_2] -= dL2;
							}
							double sum1 =0,sum2=0;
	
							for( int y = 0; y < bilayerComp.nlipidTypes; y++ ) sum1 += comp1[y];
							for( int y = 0; y < bilayerComp.nlipidTypes; y++ ) sum2 += comp2[y];
	
							if( comp1[from_1] < 0 || comp2[from_1] < 0 || 
							    comp2[from_2] < 0 || comp1[from_2] < 0 || sum1<1e-3 || sum2 < 1e-3)
							{
								comp1[from_1] += dL;
								comp2[from_1] -= dL;
								if( swap_only )
								{
									comp1[from_2] -= dL2;
									comp2[from_2] += dL2;
								}
							}	
							else
							{
								set_g0_from_f(f1);				
								set_g0_from_f(f2);				
				
								double E1 =  faceEnergy( f1, r, NULL, 0 );
								E1 += faceEnergy( f2, r, NULL, 0 );
	#if 0 
								if( !strcasecmp( bilayerComp.names[x], "DOPE") && leaf == 1)	
								{
									double k;
									double c0_1 = c( f1, 1.0/3.0, 1.0/3.0, r, &k );
									double c0_2 = c( f2, 1.0/3.0, 1.0/3.0, r, &k );
									printf("Moving %lf lipid area %s with c0 %lf from curvature %le to %le dE: %le\n", dL, bilayerComp.names[x], bilayerComp.c0[x],
									c0_1, c0_2, E1-E0 );
	
	#if 0
									if( E1 < E0 && c0_2 > c0_1 && dL > 0)
									{
										double E1 =  faceEnergy( f1, r, NULL, 0 );
										E1 += faceEnergy( f2, r, NULL, 0 );
										comp1[x] += dL;
										comp2[x] -= dL;
										set_g0_from_f(f1);				
										set_g0_from_f(f2);				
										double E0 =  faceEnergy( f1, r, NULL, 0 );
										E0 += faceEnergy( f2, r, NULL, 0 );
										printf("Confirm: dE: %le\n", E1-E0 );
										exit(1);
									}
	#endif
								}
	
	#endif
								double p = exp(-beta*(E1-E0));
				
								double rn = ((double)rand())/(double)RAND_MAX;
				
								double rat1 = theTriangles[t1].composition.A_inst / theTriangles[t1].composition.A0;
								double rat2 = theTriangles[t2].composition.A_inst / theTriangles[t2].composition.A0;
				
								if( rn < p && rat1 > 0.25 && rat2 > 0.25 )
								{
	//								printf("Moving %le (out of %le %le) %s from %d to %d.\n", dL, comp1[x], comp2[x], bilayerComp.names[x], t1, t2 );
									nacc++;
									if( move_tracer_1 )
									{
										theTriangles[t1].composition.tracer = 0;
										theTriangles[t2].composition.tracer = 1;
										
										if( tracer ) *tracer = t2;
									}
									if( move_tracer_2 )
									{
										theTriangles[t2].composition.tracer = 0;
										theTriangles[t1].composition.tracer = 1;
	
										if( tracer ) *tracer = t1;
									}
									
									lipids_moved += 1;
								}
								else
								{
									nrej++;
							
	
									comp1[from_1] += dL;
									comp2[from_1] -= dL;
									if( swap_only )
									{
										comp1[from_2] -= dL2;
										comp2[from_2] += dL2;
									}
								
									set_g0_from_f(f1);				
									set_g0_from_f(f2);				
								}
							}
				
							n_move -= 1;
						}
					}
				}
			}
		}
	}
//	printf("Lipids moved: %lf\n", lipids_moved );
/*
	for( int t1 = 0; t1 < nt; t1++ )
	{
		for( int bt = 0; bt < 3; bt++ )
		{
			int t2 = theTriangles[t1].border_tri[bt];

			if( t2 < t1 ) continue;

			double E0 = 0;

			int f1 = theTriangles[t1].f;
			int f2 = theTriangles[t2].f;

			E0 =  faceEnergy( f1, r, NULL, 0 );
			E0 += faceEnergy( f2, r, NULL, 0 );
			
			// move lipid material between triangles.

			double nL1 = theTriangles[t1].nlipids * theTriangles[t1].f_lipids;
			double nL2 = theTriangles[t2].nlipids * theTriangles[t2].f_lipids;	

			double dL = 2*((double)rand() / (double)RAND_MAX - 0.5 ) * theTriangles[t1].nlipids * nL_move_average;
	
			theTriangles[t1].f_lipids = (nL1 + dL) / theTriangles[t1].nlipids;
			theTriangles[t2].f_lipids = (nL2 - dL) / theTriangles[t2].nlipids;

			if( theTriangles[t1].f_lipids < F_CUTOFF || 
			    theTriangles[t2].f_lipids < F_CUTOFF )
			{
				theTriangles[t1].f_lipids = nL1 / theTriangles[t1].nlipids;
				theTriangles[t2].f_lipids = nL2 / theTriangles[t2].nlipids;
			}	
			else
			{
				set_g0_from_f(f1);				
				set_g0_from_f(f2);				
	
				double E1 =  faceEnergy( f1, r, NULL, 0 );
				E1 += faceEnergy( f2, r, NULL, 0 );
	
				double p = exp(-beta*(E1-E0));

				double rn = ((double)rand())/(double)RAND_MAX;
	
				if( rn < p )
				{
					nacc++;
				}
				else
				{
					nrej++;
					theTriangles[t1].f_lipids = nL1 / theTriangles[t1].nlipids;
					theTriangles[t2].f_lipids = nL2 / theTriangles[t2].nlipids;
				
					set_g0_from_f(f1);				
					set_g0_from_f(f2);				
				}
			}
		}
	}
	
	printf("nacc: %d nrej: %d\n", nacc, nrej );
	*/
	
#ifdef PARALLEL
//	lipidSync();
	lipidBroadcast();
	set_g0_from_f();
#endif
/*
	printf("Triangle 0 inner leaflet has composition: ");
	for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
		printf(" %lf %s", theTriangles[0].composition.innerLeaflet[x], bilayerComp.names[x] );  
	printf("\n");
	printf("Triangle 0 outer leaflet composition: ");
	for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
		printf(" %lf %s", theTriangles[0].composition.outerLeaflet[x], bilayerComp.names[x] );  
	printf("\n");
*/
}

void surface::measureLipidCurvature( double *r, int pre_equil )
{
	double inst_sum_c[2*bilayerComp.nlipidTypes];
	double inst_num_c[2*bilayerComp.nlipidTypes];

	memset( inst_sum_c, 0, sizeof(double) * 2 * bilayerComp.nlipidTypes );
	memset( inst_num_c, 0, sizeof(double) * 2 * bilayerComp.nlipidTypes );

	for( int t = 0; t < nt; t++ )
	{
		int f = theTriangles[t].f;
		
		double k;
		double cmid = c(f,1.0/3.0,1.0/3.0, r,&k);	

		for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
		{
			inst_sum_c[2*x] += theTriangles[t].composition.innerLeaflet[x] * (cmid * -1);
			inst_num_c[2*x] += theTriangles[t].composition.innerLeaflet[x];

			inst_sum_c[2*x+1] += theTriangles[t].composition.outerLeaflet[x] * cmid;
			inst_num_c[2*x+1] += theTriangles[t].composition.outerLeaflet[x];
		}	
	}

	printf("(Inst) averaged lipid curvature:" );

	for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
	{
		printf(" %s %le (inner) %le (outer)",
			bilayerComp.names[x], inst_sum_c[2*x+0] / inst_num_c[2*x+0],
					      inst_sum_c[2*x+1] / inst_num_c[2*x+1] );
		if( !pre_equil )
		{
			bilayerComp.sum_C[2*x] += inst_sum_c[2*x+0];
			bilayerComp.num_C[2*x] += inst_num_c[2*x+0];

			bilayerComp.sum_C[2*x+1] += inst_sum_c[2*x+1];
			bilayerComp.num_C[2*x+1] += inst_num_c[2*x+1];
		}
	}
	printf("\n");
	
	if( ! pre_equil )
	{
		printf("(Running) averaged lipid curvature:");
		for( int x = 0; x < bilayerComp.nlipidTypes; x++ )
		{
			printf(" %s %le (inner) %le (outer)",
				bilayerComp.names[x], 
					bilayerComp.sum_C[2*x+0] / bilayerComp.num_C[2*x+0],
					bilayerComp.sum_C[2*x+1] / bilayerComp.num_C[2*x+1]
						       );
		}
		printf("\n");
	}
}


