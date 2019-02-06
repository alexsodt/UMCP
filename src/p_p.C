/*
inter-complex particle-particle interactions
*/

#include "interp.h"
#include "pcomplex.h"
#include "p_p.h"
#include <math.h>
#include "parallel.h"
#include <string.h>
#include "global_boxing.h"
#include "units.h"
#include "mutil.h"

//#define DISABLE_PP
//#define DISABLE_ELASTIC
int enable_elastic_interior = 1;
static double cut = 1.122462048; // 2^(1/6)	
static double eps = 1;
static double SMOOTH = 0.9;
static double SMOOTH_THRESH = 1e-2;

static int *elastic_exclusion_list = NULL;
static int nexcluded = 0;

/*
int pos_pot = 12;
int neg_pot = 6;
*/
int pos_pot = 4;
int neg_pot = 2;
static double use_max_sigma = 0;

double rmin_mult = pow(neg_pot,1.0/(neg_pot-pos_pot)) * pow(pos_pot,1.0/(pos_pot-neg_pot) );
	
void local_setup( surface *theSurface, pcomplex **allComplexes, int ncomplex )
{
	double max_sigma = theSurface->PBC_vec[0][0] / 100;

	for( int c1 = 0; c1 < ncomplex; c1++ )
	{
		for( int p1 = 0; p1 < allComplexes[c1]->nsites; p1++ )
		{
			double sigma1 = allComplexes[c1]->sigma[p1];
			if( sigma1 > max_sigma )
				max_sigma = sigma1;
			sigma1 = allComplexes[c1]->att_sigma[p1];
			if( sigma1 > max_sigma )
				max_sigma = sigma1;
		}
	}

	use_max_sigma = max_sigma;

	setup_global_boxing(2 * max_sigma, theSurface->PBC_vec  );
}

/* potential */
double PP_V( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex )
{
	double *alphas = rsurf + theSurface->nv*3;

#ifdef DISABLE_PP
	return 0;
#endif

	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->setrall(theSurface, rsurf);

	double v = 0;

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c1 = par_info.complexes[cx];
		for( int c2 = c1+1; c2 < ncomplex; c2++ )
		{
			for( int p1 = 0; p1 < allComplexes[c1]->nsites; p1++ )
			{
				double max_sigma1 = allComplexes[c1]->sigma[p1];

				if( allComplexes[c1]->att_sigma[p1] > max_sigma1 )
					max_sigma1 = allComplexes[c1]->att_sigma[p1];

				double sigma1 = allComplexes[c1]->sigma[p1];
				double att_sigma1 = allComplexes[c1]->att_sigma[p1];

				if( !(max_sigma1>0) ) continue;

				for( int p2 = 0; p2 < allComplexes[c2]->nsites; p2++ )
				{	
					double max_sigma2 = allComplexes[c2]->sigma[p2];

					if( allComplexes[c2]->att_sigma[p2] > max_sigma2 )
						max_sigma2 = allComplexes[c2]->att_sigma[p2];

					double sigma2 = allComplexes[c2]->sigma[p2];
					double att_sigma2 = allComplexes[c2]->att_sigma[p2];

					if( !(max_sigma2>0) ) continue;
	
					double dr[3] = { 
						allComplexes[c1]->rall[3*p1+0] - allComplexes[c2]->rall[3*p2+0],
						allComplexes[c1]->rall[3*p1+1] - allComplexes[c2]->rall[3*p2+1],
						allComplexes[c1]->rall[3*p1+2] - allComplexes[c2]->rall[3*p2+2] };
	
					theSurface->wrapPBC( dr, alphas );
	
					double sigma = sigma1+sigma2;
					double att_sigma = att_sigma1+att_sigma2;
					double att_eps = (allComplexes[c1]->att_eps[p1] + allComplexes[c2]->att_eps[p2])/2;

					if( allComplexes[c2]->isElastic() && allComplexes[c1]->isElastic() && sigma1>0 && sigma2>0)
					{
						double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
						if( r < sigma && enable_elastic_interior )
						{
							double r2 = (sigma/(rmin_mult*r))*(sigma/(rmin_mult*r));
							double r4 = r2*r2;
							double r6 = r2*r4;

							v += eps + 4.0 * eps * (r4-r2);
							printf("ELASTIC POT %le\n", eps + 4.0 * eps * (r4-r2) );
						}
					}
					else
					{
						if( sigma1 > 0 && sigma2 > 0 )
						{
							double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
							double r2 = (sigma/r)*(sigma/r);
							double r4 = r2*r2;
							double r6 = r2*r4;
	
							if( r < rmin_mult * sigma )
								v += eps + 4.0 * eps * (r4-r2);
						}
						else if( att_sigma1 > 0 && att_sigma2 > 0 )
						{
							double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
							double r0 = SMOOTH * att_sigma;
							double r1 = att_sigma;
							
							if( r < SMOOTH * att_sigma * (1+SMOOTH_THRESH) )
								v += att_eps; 	
							else if( r < r1 * (1-SMOOTH_THRESH) )
								v += att_eps * (1 - 1 / (1 + exp( (r1-r0)/(r-r0) - (r1-r0)/(r1-r) ) ));
						}
					}
//						v += eps + 4.0 * eps * (r6*r6-r6);
				}
			}
		}
	}

	return v;
}

/* gradient */
double PP_G( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, double *mesh_grad )
{
	double *alphas = rsurf + theSurface->nv*3;
#ifdef DISABLE_PP
	return 0;
#endif

	double v = 0;

	int *ourp = (int *)malloc( sizeof(int) * ncomplex );
	memset( ourp, 0, sizeof(int) * ncomplex );

	for( int cx = 0; cx < par_info.nc; cx++ )
		ourp[par_info.complexes[cx]] = 1;

//	for( int c1 = 0; c1 < ncomplex; c1++ )
	
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c1 = par_info.complexes[cx];

		for( int c2 = 0; c2 < ncomplex; c2++ )
		{
			if( c1 == c2 ) 
				continue;
			if( ourp[c1] && ourp[c2] && c2 < c1 ) 
				continue;
			for( int p1 = 0; p1 < allComplexes[c1]->nsites; p1++ )
			{
				double max_sigma1 = allComplexes[c1]->sigma[p1];

				if( allComplexes[c1]->att_sigma[p1] > max_sigma1 )
					max_sigma1 = allComplexes[c1]->att_sigma[p1];

				double sigma1 = allComplexes[c1]->sigma[p1];
				double att_sigma1 = allComplexes[c1]->att_sigma[p1];

				if( !(max_sigma1>0) ) continue;

				for( int p2 = 0; p2 < allComplexes[c2]->nsites; p2++ )
				{	
					double max_sigma2 = allComplexes[c2]->sigma[p2];

					if( allComplexes[c2]->att_sigma[p2] > max_sigma2 )
						max_sigma2 = allComplexes[c2]->att_sigma[p2];

					double sigma2 = allComplexes[c2]->sigma[p2];
					double att_sigma2 = allComplexes[c2]->att_sigma[p2];

					if( !(max_sigma2>0) ) continue;
	
					double dr[3] = { 
						allComplexes[c1]->rall[3*p1+0] - allComplexes[c2]->rall[3*p2+0],
						allComplexes[c1]->rall[3*p1+1] - allComplexes[c2]->rall[3*p2+1],
						allComplexes[c1]->rall[3*p1+2] - allComplexes[c2]->rall[3*p2+2] };
	
					theSurface->wrapPBC( dr, alphas );
	
					double sigma = sigma1+sigma2;
					double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

					if( allComplexes[c2]->isElastic() && allComplexes[c1]->isElastic() && sigma1 > 0 && sigma2 > 0 )
					{
						if( r < sigma && enable_elastic_interior )
						{
							double r2 = (sigma/(rmin_mult*r))*(sigma/(rmin_mult*r));
							double r4 = r2*r2;
							double r6 = r2*r4;

							double local_fac = 1.0;
							if( ourp[c1] != ourp[c2] ) local_fac = 0.5; // calculation is duplicated on another process.
							v += (eps + 4.0 * eps * (r4-r2)) * local_fac;
							printf("ELASTIC FRC, V %le\n", eps + 4.0 * eps * (r4-r2) * local_fac );
			
							double dvdr = 4 * eps * (-4 * r4 + 2 * r2)/r;
							
							double null_f[3] = { 0,0,0}; // zero explicit force on the normal.
			
							double f1[3] = {  dvdr*dr[0]/r,  dvdr*dr[1]/r,  dvdr*dr[2]/r };
							double f2[3] = { -f1[0], -f1[1], -f1[2] };
		
	//						printf("%d %d and %d %d dist %le\n", p1, c1, p2, c2, r );
		
							if( ourp[c1] )
							{
								allComplexes[c1]->save_grad[3*p1+0] += f1[0];
								allComplexes[c1]->save_grad[3*p1+1] += f1[1];
								allComplexes[c1]->save_grad[3*p1+2] += f1[2];
							}
						
							if( ourp[c2] )
							{
								allComplexes[c2]->save_grad[3*p2+0] += f2[0];
								allComplexes[c2]->save_grad[3*p2+1] += f2[1];
								allComplexes[c2]->save_grad[3*p2+2] += f2[2];
							}
						}
					}
					else
					{
						double dvdr = 0;
						double local_fac = 1.0;
						if( ourp[c1] != ourp[c2] ) local_fac = 0.5; // calculation is duplicated on another process.

						if( sigma1 > 0 && sigma2 > 0 )
						{
							if( r < rmin_mult * sigma )
							{
								double r2 = (sigma/r)*(sigma/r);
								double r4 = r2*r2;
								double r6 = r2*r4;
			
								v += (eps + 4.0 * eps * (r4-r2)) * local_fac;
			
								dvdr = 4 * eps * (-4 * r4 + 2 * r2)/r;
							}		
						}
						else if( att_sigma1 > 0 && att_sigma2 > 0 )
						{
							double att_sigma = att_sigma1+att_sigma2;
							double att_eps = (allComplexes[c1]->att_eps[p1] + allComplexes[c2]->att_eps[p2])/2;
							double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
							double r0 = SMOOTH * att_sigma;
							double r1 = att_sigma;
							
							if( r < SMOOTH * att_sigma * (1+SMOOTH_THRESH) )
								v += att_eps * local_fac; 	
							else if( r < r1 * (1-SMOOTH_THRESH) )
							{
								double eval = exp( (r1-r0)/(r-r0) - (r1-r0)/(r1-r) );
								v += att_eps * (1 - 1 / (1 + exp( (r1-r0)/(r-r0) - (r1-r0)/(r1-r) ) )) * local_fac;
	
								dvdr = att_eps * eval * ( -(r1-r0)/(r-r0)/(r-r0)-(r1-r0)/(r1-r)/(r1-r))/(1+eval)/(1+eval); 
							}
						}

						double null_f[3] = { 0,0,0}; // zero explicit force on the normal.
		
						double f1[3] = {  dvdr*dr[0]/r,  dvdr*dr[1]/r,  dvdr*dr[2]/r };
						double f2[3] = { -f1[0], -f1[1], -f1[2] };
	
						if( fabs(dvdr) > 0 )
						{	
							if( p1 < allComplexes[c1]->nattach)
							{
								if( ourp[c1] )
								{
								theSurface->pointGradient( allComplexes[c1]->grad_fs[p1], 
										allComplexes[c1]->grad_puv[2*p1+0], 
										allComplexes[c1]->grad_puv[2*p1+1],
											rsurf, mesh_grad, allComplexes[c1]->save_grad+2*p1, f1, null_f ); 
								}
							}
							else
							{
								if( ourp[c1] )
								{
									allComplexes[c1]->save_grad[3*p1+0] += f1[0];
									allComplexes[c1]->save_grad[3*p1+1] += f1[1];
									allComplexes[c1]->save_grad[3*p1+2] += f1[2];
								}
							}
						
							if( p2 < allComplexes[c2]->nattach )
							{
								if( ourp[c2] )
								{
									theSurface->pointGradient( allComplexes[c2]->grad_fs[p2], 
										allComplexes[c2]->grad_puv[2*p2+0], 
										allComplexes[c2]->grad_puv[2*p2+1],
											rsurf, mesh_grad, allComplexes[c2]->save_grad+2*p2, f2, null_f ); 
								}
							}
							else
							{
								if( ourp[c2] )
								{
									allComplexes[c2]->save_grad[3*p2+0] += f2[0];
									allComplexes[c2]->save_grad[3*p2+1] += f2[1];
									allComplexes[c2]->save_grad[3*p2+2] += f2[2];
								}
							}
						}
					}
				}
			}
		}
	}
	free(ourp);

	return v;
}

/* gradient */
double Boxed_PP_V( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex )
{
	double *alphas = rsurf + theSurface->nv*3;
#ifdef DISABLE_PP
	return 0;
#endif

	for( int c = 0; c < ncomplex; c++ )
		allComplexes[c]->setrall(theSurface, rsurf);

	if( ! global_boxing_init || !boxing )
		local_setup( theSurface, allComplexes, ncomplex );

	boxing->setPBC( theSurface->PBC_vec, alphas );

	int ntotp = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		ntotp += allComplexes[c]->nsites;
	}
	int *to_clear = (int *)malloc( sizeof(int) * ntotp );
	double v = 0;

	int *ourp = (int *)malloc( sizeof(int) * ncomplex );
	memset( ourp, 0, sizeof(int) * ncomplex );
	int *nearlist = (int *)malloc( sizeof(int) * ntotp );	

	for( int cx = 0; cx < par_info.nc; cx++ )
		ourp[par_info.complexes[cx]] = 1;



	int *complex_for_id = (int *)malloc( sizeof(int) * ntotp );
	int *subp_for_id = (int *)malloc( sizeof(int) * ntotp );

	int id = 0;
	for( int c1 = 0; c1 < ncomplex; c1++ )
	{
		for( int p = 0; p < allComplexes[c1]->nsites; p++ )
		{
			complex_for_id[id] = c1;
			subp_for_id[id] = p;
			to_clear[id] = boxing->addp(allComplexes[c1]->rall+3*p, id );
			id++;
		}
	}

//	for( int c1 = 0; c1 < ncomplex; c1++ )
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c1 = par_info.complexes[cx];
		

		for( int p1 = 0; p1 < allComplexes[c1]->nsites; p1++ )
		{
			double max_sigma1 = allComplexes[c1]->sigma[p1];

			if( allComplexes[c1]->att_sigma[p1] > max_sigma1 )
				max_sigma1 = allComplexes[c1]->att_sigma[p1];

			double sigma1 = allComplexes[c1]->sigma[p1];
			double att_sigma1 = allComplexes[c1]->att_sigma[p1];

			if( !(max_sigma1>0) ) continue;
	
			int np  = boxing->getNearPts( allComplexes[c1]->rall+3*p1, nearlist, 1.25 * rmin_mult * (use_max_sigma + max_sigma1) );

			for( int x = 0; x < np; x++ )
			{
				int c2 = complex_for_id[nearlist[x]];
				int p2 = subp_for_id[nearlist[x]];
			
				if( c1 == c2 ) 
					continue;
				if( c2 < c1 ) 
					continue;
	
				double max_sigma2 = allComplexes[c2]->sigma[p2];

				if( allComplexes[c2]->att_sigma[p2] > max_sigma2 )
					max_sigma2 = allComplexes[c2]->att_sigma[p2];

				double sigma2 = allComplexes[c2]->sigma[p2];
				double att_sigma2 = allComplexes[c2]->att_sigma[p2];

				if( !(max_sigma2>0) ) continue;
	
				double dr[3] = { 
					allComplexes[c1]->rall[3*p1+0] - allComplexes[c2]->rall[3*p2+0],
					allComplexes[c1]->rall[3*p1+1] - allComplexes[c2]->rall[3*p2+1],
					allComplexes[c1]->rall[3*p1+2] - allComplexes[c2]->rall[3*p2+2] };

				theSurface->wrapPBC( dr, alphas );

				double sigma = sigma1+sigma2;
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

				double att_sigma = att_sigma1+att_sigma2;
				double att_eps = (allComplexes[c1]->att_eps[p1] + allComplexes[c2]->att_eps[p2])/2;

				if( allComplexes[c2]->isElastic() && allComplexes[c1]->isElastic() && sigma1>0 && sigma2>0)
				{
					double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
					if( r < sigma && enable_elastic_interior )
					{
						double r2 = (sigma/(rmin_mult*r))*(sigma/(rmin_mult*r));
						double r4 = r2*r2;
						double r6 = r2*r4;

						v += eps + 4.0 * eps * (r4-r2);
					}
				}
				else
				{
					if( sigma1 > 0 && sigma2 > 0 )
					{
						double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
						double r2 = (sigma/r)*(sigma/r);
						double r4 = r2*r2;
						double r6 = r2*r4;
	
						if( r < rmin_mult * sigma )
							v += eps + 4.0 * eps * (r4-r2);
					}
					else if( att_sigma1 > 0 && att_sigma2 > 0 )
					{
						double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
						double r0 = SMOOTH * att_sigma;
						double r1 = att_sigma;
						
						if( r < SMOOTH * r1 * (1+SMOOTH_THRESH) )
							v += att_eps; 	
						else if( r < r1 * (1-SMOOTH_THRESH) )
							v += att_eps * (1- 1 / (1 + exp( (r1-r0)/(r-r0) - (r1-r0)/(r1-r) ) ));
					}
				}
			}
		}
	}
	
	free(nearlist);
	free(ourp);
	free(complex_for_id);
	free(subp_for_id);
	
	boxing->clearBoxing(to_clear, ntotp);
	free(to_clear);
	return v;
}

/* gradient */
double Boxed_PP_G( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, double *mesh_grad )
{
	double *alphas = rsurf + theSurface->nv*3;
#ifdef DISABLE_PP
	return 0;
#endif
	if( ! global_boxing_init || !boxing )
		local_setup( theSurface, allComplexes, ncomplex );
	boxing->setPBC( theSurface->PBC_vec, alphas );


	int ntotp = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		ntotp += allComplexes[c]->nsites;
	}
	int *to_clear = (int *)malloc( sizeof(int) * ntotp );
	double v = 0;

	int *ourp = (int *)malloc( sizeof(int) * ncomplex );
	memset( ourp, 0, sizeof(int) * ncomplex );
	int *nearlist = (int *)malloc( sizeof(int) * ntotp );	

	for( int cx = 0; cx < par_info.nc; cx++ )
		ourp[par_info.complexes[cx]] = 1;



	int *complex_for_id = (int *)malloc( sizeof(int) * ntotp );
	int *subp_for_id = (int *)malloc( sizeof(int) * ntotp );

	int id = 0;
	for( int c1 = 0; c1 < ncomplex; c1++ )
	{
		for( int p = 0; p < allComplexes[c1]->nsites; p++ )
		{
			complex_for_id[id] = c1;
			subp_for_id[id] = p;
			to_clear[id] = boxing->addp(allComplexes[c1]->rall+3*p, id );
			id++;
		}
	}

//	for( int c1 = 0; c1 < ncomplex; c1++ )
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c1 = par_info.complexes[cx];
		

		for( int p1 = 0; p1 < allComplexes[c1]->nsites; p1++ )
		{
			double max_sigma1 = allComplexes[c1]->sigma[p1];

			if( allComplexes[c1]->att_sigma[p1] > max_sigma1 )
				max_sigma1 = allComplexes[c1]->att_sigma[p1];

			double sigma1 = allComplexes[c1]->sigma[p1];
			double att_sigma1 = allComplexes[c1]->att_sigma[p1];

			if( !(max_sigma1>0) ) continue;
	
			int np  = boxing->getNearPts( allComplexes[c1]->rall+3*p1, nearlist, 1.25 * rmin_mult * (use_max_sigma + max_sigma1) );

			for( int x = 0; x < np; x++ )
			{
				int c2 = complex_for_id[nearlist[x]];
				int p2 = subp_for_id[nearlist[x]];
			
				if( c1 == c2 ) 
					continue;
				if( ourp[c1] && ourp[c2] && c2 < c1 ) 
					continue;
/*
			for( int c2 = 0; c2 < ncomplex; c2++ )
			
				if( c1 == c2 ) 
					continue;
				if( ourp[c1] && ourp[c2] && c2 < c1 ) 
					continue;
				for( int p1 = 0; p1 < allComplexes[c1]->nsites; p1++ )
				
*/
	
				double max_sigma2 = allComplexes[c2]->sigma[p2];

				if( allComplexes[c2]->att_sigma[p2] > max_sigma2 )
					max_sigma2 = allComplexes[c2]->att_sigma[p2];

				double sigma2 = allComplexes[c2]->sigma[p2];
				double att_sigma2 = allComplexes[c2]->att_sigma[p2];

				if( !(max_sigma2>0) ) continue;
	
				double dr[3] = { 
					allComplexes[c1]->rall[3*p1+0] - allComplexes[c2]->rall[3*p2+0],
					allComplexes[c1]->rall[3*p1+1] - allComplexes[c2]->rall[3*p2+1],
					allComplexes[c1]->rall[3*p1+2] - allComplexes[c2]->rall[3*p2+2] };

				theSurface->wrapPBC( dr, alphas );

				double att_sigma = att_sigma1+att_sigma2;
				double att_eps = (allComplexes[c1]->att_eps[p1] + allComplexes[c2]->att_eps[p2])/2;
				double sigma = sigma1+sigma2;
				double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);

				if( allComplexes[c2]->isElastic() && allComplexes[c1]->isElastic() && sigma1 > 0 && sigma2 > 0 )
				{
					if( r < sigma && enable_elastic_interior )
					{

						double r2 = (sigma/(rmin_mult*r))*(sigma/(rmin_mult*r));
						double r4 = r2*r2;
						double r6 = r2*r4;

						double local_fac = 1.0;
						if( ourp[c1] != ourp[c2] ) local_fac = 0.5; // calculation is duplicated on another process.
						v += (eps + 4.0 * eps * (r4-r2)) * local_fac;
			
						double dvdr = 4 * eps * (-4 * r4 + 2 * r2)/r;
						
						double null_f[3] = { 0,0,0}; // zero explicit force on the normal.
			
						double f1[3] = {  dvdr*dr[0]/r,  dvdr*dr[1]/r,  dvdr*dr[2]/r };
						double f2[3] = { -f1[0], -f1[1], -f1[2] };
		
	//					printf("%d %d and %d %d dist %le\n", p1, c1, p2, c2, r );
		
						if( ourp[c1] )
						{
							allComplexes[c1]->save_grad[3*p1+0] += f1[0];
							allComplexes[c1]->save_grad[3*p1+1] += f1[1];
							allComplexes[c1]->save_grad[3*p1+2] += f1[2];
						}
					
						if( ourp[c2] )
						{
							allComplexes[c2]->save_grad[3*p2+0] += f2[0];
							allComplexes[c2]->save_grad[3*p2+1] += f2[1];
							allComplexes[c2]->save_grad[3*p2+2] += f2[2];
						}
					}
				}
				else
				{
					double dvdr = 0;
					double local_fac = 1.0;
					if( ourp[c1] != ourp[c2] ) local_fac = 0.5; // calculation is duplicated on another process.

					if( sigma1 > 0 && sigma2 > 0 )
					{
						if( r < rmin_mult * sigma )
						{
							double r2 = (sigma/r)*(sigma/r);
							double r4 = r2*r2;
							double r6 = r2*r4;
			
							v += (eps + 4.0 * eps * (r4-r2)) * local_fac;
			
							dvdr = 4 * eps * (-4 * r4 + 2 * r2)/r;
						}
					}
					else if( att_sigma1 > 0 && att_sigma2 > 0 )
					{
						double att_sigma = att_sigma1+att_sigma2;
						double att_eps = (allComplexes[c1]->att_eps[p1] + allComplexes[c2]->att_eps[p2])/2;
						double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
						double r0 = SMOOTH * att_sigma;
						double r1 = att_sigma;
						
						if( r < SMOOTH * att_sigma * (1+SMOOTH_THRESH) )
							v += att_eps * local_fac; 	
						else if( r < r1 * (1-SMOOTH_THRESH) )
						{
							double eval = exp( (r1-r0)/(r-r0) - (r1-r0)/(r1-r) );
							v += att_eps * (1.0 - 1.0/ (1 + eval )) * local_fac;

							dvdr = att_eps * eval * ( -(r1-r0)/(r-r0)/(r-r0)-(r1-r0)/(r1-r)/(r1-r))/(1+eval)/(1+eval); 
							if( !(dvdr < 1 ) && !(dvdr > 0 ) )
							{
								printf("DBG\n");
							}
						}
					}

					double null_f[3] = { 0,0,0}; // zero explicit force on the normal.
		
					double f1[3] = {  dvdr*dr[0]/r,  dvdr*dr[1]/r,  dvdr*dr[2]/r };
					double f2[3] = { -f1[0], -f1[1], -f1[2] };
	
					if( fabs(dvdr) > 0 )
					{	
						if( p1 < allComplexes[c1]->nattach)
						{
							if( ourp[c1] )
							{
							theSurface->pointGradient( allComplexes[c1]->grad_fs[p1], 
									allComplexes[c1]->grad_puv[2*p1+0], 
									allComplexes[c1]->grad_puv[2*p1+1],
										rsurf, mesh_grad, allComplexes[c1]->save_grad+2*p1, f1, null_f ); 
							}
						}
						else
						{
							if( ourp[c1] )
							{
								allComplexes[c1]->save_grad[3*p1+0] += f1[0];
								allComplexes[c1]->save_grad[3*p1+1] += f1[1];
								allComplexes[c1]->save_grad[3*p1+2] += f1[2];
							}
						}
					
						if( p2 < allComplexes[c2]->nattach )
						{
							if( ourp[c2] )
							{
								theSurface->pointGradient( allComplexes[c2]->grad_fs[p2], 
									allComplexes[c2]->grad_puv[2*p2+0], 
									allComplexes[c2]->grad_puv[2*p2+1],
										rsurf, mesh_grad, allComplexes[c2]->save_grad+2*p2, f2, null_f ); 
							}
						}
						else
						{
							if( ourp[c2] )
							{
								allComplexes[c2]->save_grad[3*p2+0] += f2[0];
								allComplexes[c2]->save_grad[3*p2+1] += f2[1];
								allComplexes[c2]->save_grad[3*p2+2] += f2[2];
							}
						}
					}
				}
			}
		}
	}
	
	free(nearlist);
	free(ourp);
	free(complex_for_id);
	free(subp_for_id);
	
	boxing->clearBoxing(to_clear, ntotp);
	free(to_clear);
	return v;
}

double handleElasticCollisions( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, double time_step )
{
#ifdef DISABLE_ELASTIC
	return 0.0;
#endif

	double *alphas = rsurf + theSurface->nv*3;
	if( ! global_boxing_init || !boxing )
		local_setup( theSurface, allComplexes, ncomplex );
	boxing->setPBC( theSurface->PBC_vec, alphas );

	int ntotp = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		ntotp += allComplexes[c]->nsites;
	}
	int *to_clear = (int *)malloc( sizeof(int) * ntotp );
	double v = 0;

	int *ourp = (int *)malloc( sizeof(int) * ncomplex );
	memset( ourp, 0, sizeof(int) * ncomplex );
	int *nearlist = (int *)malloc( sizeof(int) * ntotp );	

	for( int cx = 0; cx < par_info.nc; cx++ )
		ourp[par_info.complexes[cx]] = 1;



	int *complex_for_id = (int *)malloc( sizeof(int) * ntotp );
	int *subp_for_id = (int *)malloc( sizeof(int) * ntotp );

	int id = 0;

	for( int c1 = 0; c1 < ncomplex; c1++ )
	{
		for( int p = 0; p < allComplexes[c1]->nsites; p++ )
		{
			complex_for_id[id] = c1;
			subp_for_id[id] = p;
			to_clear[id] = boxing->addp(allComplexes[c1]->rall+3*p, id );
			id++;
		}
	}

	double n_not_ok = 0;

//	for( int c1 = 0; c1 < ncomplex; c1++ )
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c1 = par_info.complexes[cx];

		if( ! allComplexes[c1]->isElastic() ) 
			continue;

		// starts at nattach: we only do solution elastic collisions.
		for( int p1 = allComplexes[c1]->nattach; p1 < allComplexes[c1]->nsites; p1++ )
		{
			double sigma1 = allComplexes[c1]->sigma[p1];
	
			int np  = boxing->getNearPts( allComplexes[c1]->rall+3*p1, nearlist, 1.25 * rmin_mult * (use_max_sigma + sigma1) );

			for( int x = 0; x < np; x++ )
			{
				int c2 = complex_for_id[nearlist[x]];
				int p2 = subp_for_id[nearlist[x]];
		
				if( ! allComplexes[c2]->isElastic() ) 
					continue;
			
				if( c1 == c2 ) 
					continue;

				if( ourp[c1] && ourp[c2] && c2 < c1 ) 
					continue;

				if( p2 >=  allComplexes[c1]->nattach )
				{
					// solution-solution collision: only thing currently implemented (1/2019-AJS)	
					// two objects, r1 + v1 with mass m1, r2 + v2 with mass m2
					
					double *r1 = allComplexes[c1]->rall+3*p1;
					double *r2 = allComplexes[c2]->rall+3*p2;
					double m1 = allComplexes[c1]->mass[p1];
					double m2 = allComplexes[c2]->mass[p2];

					double dr[3] = { r1[0]-r2[0], r1[1]-r2[1],r1[2]-r2[2]};
					theSurface->wrapPBC( dr, alphas );
					double dr2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
	
					double p1r = allComplexes[c1]->sigma[p1];
					double p2r = allComplexes[c2]->sigma[p2];
					double r02 = (p1r+p2r)*(p1r+p2r);

					if( dr2 < r02 )
					{

						double v1[3] = { 
							allComplexes[c1]->qdot[3*p1+0], 
							allComplexes[c1]->qdot[3*p1+1], 
							allComplexes[c1]->qdot[3*p1+2] };
						double v2[3] = { 
							allComplexes[c2]->qdot[3*p2+0], 
							allComplexes[c2]->qdot[3*p2+1], 
							allComplexes[c2]->qdot[3*p2+2] }; 
						double dv[3] = { v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2] };
						double dv2 = 
							(dv[0])*(dv[0]) +
							(dv[1])*(dv[1]) +
							(dv[2])*(dv[2]);
						if( dv2 > 1e-30 )
						{

							double drdv = dr[0]*dv[0]+dr[1]*dv[1]+dr[2]*dv[2];
		
							if( 4*drdv*drdv-4*dv2*(dr2-r02) >= 0 )		
							{
								// time of collision.
								double t1 = (-0.5*sqrt(4*drdv*drdv-4*dv2*(dr2-r02))-drdv)/dv2;
								double t2 = (+0.5*sqrt(4*drdv*drdv-4*dv2*(dr2-r02))-drdv)/dv2;
								
								int t1_ok = t1 < 0 && t1  > -time_step;
								int t2_ok = t2 < 0 && t2  > -time_step;
		
								double use_t = t1;
								int ok = 0;
								if( t1_ok && t2_ok )
								{
									if( t1 > t2 )		
										use_t = t1;
									else
										use_t = t2;
									ok =1;
								}
								else if( t1_ok )
								{
									ok = 1;
									use_t = t1;
								}
								else if( t2_ok )
								{
									ok = 1;
									use_t = t2;
								}

								if( ok )
								{		
									double pos1_at_c[3] = { r1[0] + use_t * v1[0],
												r1[1] + use_t * v1[1],
												r1[2] + use_t * v1[2] };
	
									double pos2_at_c[3] = { r2[0] + use_t * v2[0],
												r2[1] + use_t * v2[1],
												r2[2] + use_t * v2[2] };
	
									// direction of dv.
									double dr_at_c[3] = { 
											pos1_at_c[0] - pos2_at_c[0],
											pos1_at_c[1] - pos2_at_c[1],
											pos1_at_c[2] - pos2_at_c[2] };
									theSurface->wrapPBC( dr_at_c, alphas );
									normalize(dr_at_c);
									// solve for dKE
	
									double lambda2 = (2 * m1 * (dr_at_c[0] * dv[0] +dr_at_c[1]*dv[1] + dr_at_c[2] * dv[2])/(m1+m2));
									double lambda1 = -lambda2 * (m2/m1);


									double v1_after[3] = { v1[0] + lambda1 * dr_at_c[0],
											      v1[1] + lambda1 * dr_at_c[1],			
											      v1[2] + lambda1 * dr_at_c[2] };			
									double v2_after[3] = { v2[0] + lambda2 * dr_at_c[0],
											      v2[1] + lambda2 * dr_at_c[1],			
											      v2[2] + lambda2 * dr_at_c[2] };			
	
									double KE_before = 0.5 * m1 * (v1[0]*v1[0]+v1[1]*v1[1]+v1[2]*v1[2]) +0.5 * m2 * (v2[0]*v2[0]+v2[1]*v2[1]+v2[2]*v2[2]);
									double KE_after = 0.5 * m1 * (v1_after[0]*v1_after[0]+v1_after[1]*v1_after[1]+v1_after[2]*v1_after[2]) +0.5 * m2 * (v2_after[0]*v2_after[0]+v2_after[1]*v2_after[1]+v2_after[2]*v2_after[2]);

//									printf("COLLISION %d %d vafter: %le %le %le and %le %le %le KEBA %le/%le\n", c1, c2, 
//											v1_after[0], v1_after[1], v1_after[2],	
//											v2_after[0], v2_after[1], v2_after[2], KE_before, KE_after );	

									if( ourp[c1] )
									{
										allComplexes[c1]->save_grad[3*p1+0] += -lambda1 * m1 * dr_at_c[0] / time_step;
										allComplexes[c1]->save_grad[3*p1+1] += -lambda1 * m1 * dr_at_c[1] / time_step;
										allComplexes[c1]->save_grad[3*p1+2] += -lambda1 * m1 * dr_at_c[2] / time_step;
									}
									if( ourp[c2] )
									{
										allComplexes[c2]->save_grad[3*p2+0] += -lambda2 * m2 * dr_at_c[0] / time_step;
										allComplexes[c2]->save_grad[3*p2+1] += -lambda2 * m2 * dr_at_c[1] / time_step;
										allComplexes[c2]->save_grad[3*p2+2] += -lambda2 * m2 * dr_at_c[2] / time_step;
									}
									printf("COLLISION %d %d\n", c1, c2 );
								}
								else
								{
									printf("%d %d dr: %le %le %le v1: %le %le %le v2: %le %le %le REJECTED collision with t1: %le t2: %le step %le\n", c1, c2,  
										dr[0],dr[1],dr[2], v1[0],v1[1],v1[2],v2[0],v2[1],v2[2],
										t1, t2, time_step );
									double sigma = p1r+p2r;
									double r = sqrt(dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2]);
				
									if( r < rmin_mult * sigma )
									{
										double r2 = (sigma/r)*(sigma/r);
										double r4 = r2*r2;
										double r6 = r2*r4;
						
				//						v += eps + 4.0 * eps * (r6*r6-r6);
										double local_fac = 1.0;
										if( ourp[c1] != ourp[c2] ) local_fac = 0.5; // calculation is duplicated on another process.
										v += (eps + 4.0 * eps * (r4-r2)) * local_fac;
						
										double dvdr = 4 * eps * (-4 * r4 + 2 * r2)/r;
										
										double null_f[3] = { 0,0,0}; // zero explicit force on the normal.
						
										double f1[3] = {  dvdr*dr[0]/r,  dvdr*dr[1]/r,  dvdr*dr[2]/r };
										double f2[3] = { -f1[0], -f1[1], -f1[2] };
					
				//						printf("%d %d and %d %d dist %le\n", p1, c1, p2, c2, r );
					
										if( ourp[c1] )
										{
											n_not_ok += 0.5; 
											allComplexes[c1]->save_grad[3*p1+0] += f1[0];
											allComplexes[c1]->save_grad[3*p1+1] += f1[1];
											allComplexes[c1]->save_grad[3*p1+2] += f1[2];
										}
										if( ourp[c2] )
										{
											n_not_ok += 0.5; 
											allComplexes[c2]->save_grad[3*p2+0] += f2[0];
											allComplexes[c2]->save_grad[3*p2+1] += f2[1];
											allComplexes[c2]->save_grad[3*p2+2] += f2[2];
										}
									}
								}
							}
						}
					}
				}				
			}
		}
	}
	
#ifdef PARALLEL
	ParallelSum(&n_not_ok, 1);
#endif
	
	printf("NNOT: %lf\n", n_not_ok );
	free(nearlist);
	free(ourp);
	free(complex_for_id);
	free(subp_for_id);
	
	boxing->clearBoxing(to_clear, ntotp);
	free(to_clear);

	return v;
}


#define OVERRIDE_BOXING
double nElasticCollisions( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex )
{
	double *alphas = rsurf + theSurface->nv*3;
	if( ! global_boxing_init || !boxing )
		local_setup( theSurface, allComplexes, ncomplex );
	boxing->setPBC( theSurface->PBC_vec, alphas );

	int ntotp = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		ntotp += allComplexes[c]->nsites;
	}
	int *to_clear = (int *)malloc( sizeof(int) * ntotp );
	double v = 0;

	int *ourp = (int *)malloc( sizeof(int) * ncomplex );
	memset( ourp, 0, sizeof(int) * ncomplex );
	int *nearlist = (int *)malloc( sizeof(int) * ntotp );	

	for( int cx = 0; cx < par_info.nc; cx++ )
		ourp[par_info.complexes[cx]] = 1;



	int *complex_for_id = (int *)malloc( sizeof(int) * ntotp );
	int *subp_for_id = (int *)malloc( sizeof(int) * ntotp );

	int id = 0;

	for( int c1 = 0; c1 < ncomplex; c1++ )
	{
		for( int p = 0; p < allComplexes[c1]->nsites; p++ )
		{
			complex_for_id[id] = c1;
			subp_for_id[id] = p;
			to_clear[id] = boxing->addp(allComplexes[c1]->rall+3*p, id );
			id++;
		}
	}

	double ncol = 0;

	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c1 = par_info.complexes[cx];

		if( ! allComplexes[c1]->isElastic() ) 
			continue;

		// starts at nattach: we only do solution elastic collisions.
		for( int p1 = allComplexes[c1]->nattach; p1 < allComplexes[c1]->nsites; p1++ )
		{
			double sigma1 = allComplexes[c1]->sigma[p1];

#ifdef OVERRIDE_BOXING
			for( int c2 = 0; c2 < ncomplex; c2++ )
			for( int p2 = allComplexes[c2]->nattach; p2 < allComplexes[c2]->nsites; p2++ )
			{
#else	
			int np  = boxing->getNearPts( allComplexes[c1]->rall+3*p1, nearlist, 1.25 * rmin_mult * (use_max_sigma + sigma1) );

			for( int x = 0; x < np; x++ )
			{
				int c2 = complex_for_id[nearlist[x]];
				int p2 = subp_for_id[nearlist[x]];
#endif		
				if( ! allComplexes[c2]->isElastic() ) 
					continue;
			
				if( c1 == c2 ) 
					continue;

				if( ourp[c1] && ourp[c2] && c2 < c1 ) 
					continue;

				if( p2 >=  allComplexes[c1]->nattach )
				{
					// solution-solution collision: only thing currently implemented (1/2019-AJS)	
					// two objects, r1 + v1 with mass m1, r2 + v2 with mass m2
					
					double *r1 = allComplexes[c1]->rall+3*p1;
					double *r2 = allComplexes[c2]->rall+3*p2;
					double m1 = allComplexes[c1]->mass[p1];
					double m2 = allComplexes[c2]->mass[p2];

					double dr[3] = { r1[0]-r2[0], r1[1]-r2[1],r1[2]-r2[2]};
					theSurface->wrapPBC( dr, alphas );
					double dr2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
	
					double p1r = allComplexes[c1]->sigma[p1];
					double p2r = allComplexes[c2]->sigma[p2];
					double r02 = (p1r+p2r)*(p1r+p2r);

					if( dr2 < r02 )
					{
						if( ourp[c2] )
							ncol += 1.0;
						else
							ncol += 0.5;
					}
				}				
			}
		}
	}
	
	free(nearlist);
	free(ourp);
	free(complex_for_id);
	free(subp_for_id);
	
	boxing->clearBoxing(to_clear, ntotp);
	free(to_clear);

#ifdef DISABLE_ELASTIC
	return 0;
#else
	return ncol;
#endif
}

double timePrecedingElasticCollision( surface *theSurface, double *rsurf, pcomplex **allComplexes, int ncomplex, double time_step, int *out_col_c1, int *out_col_p1, int *out_col_c2, int *out_col_p2 )
{
#ifdef DISABLE_ELASTIC
	return 0.0;
#endif

	double *alphas = rsurf + theSurface->nv*3;
	if( ! global_boxing_init || !boxing )
		local_setup( theSurface, allComplexes, ncomplex );
	boxing->setPBC( theSurface->PBC_vec, alphas );

	int ntotp = 0;
	for( int c = 0; c < ncomplex; c++ )
	{
		ntotp += allComplexes[c]->nsites;
	}
	int *to_clear = (int *)malloc( sizeof(int) * ntotp );
	double v = 0;

	int *ourp = (int *)malloc( sizeof(int) * ncomplex );
	memset( ourp, 0, sizeof(int) * ncomplex );
	int *nearlist = (int *)malloc( sizeof(int) * ntotp );	

	for( int cx = 0; cx < par_info.nc; cx++ )
		ourp[par_info.complexes[cx]] = 1;



	int *complex_for_id = (int *)malloc( sizeof(int) * ntotp );
	int *subp_for_id = (int *)malloc( sizeof(int) * ntotp );

	int id = 0;

	for( int c1 = 0; c1 < ncomplex; c1++ )
	{
		for( int p = 0; p < allComplexes[c1]->nsites; p++ )
		{
			complex_for_id[id] = c1;
			subp_for_id[id] = p;
			to_clear[id] = boxing->addp(allComplexes[c1]->rall+3*p, id );
			id++;
		}
	}

	double n_not_ok = 0;

	double time_preceding_collision = 0; 

	int col_c1=-1,col_p1=-1,col_c2=-1,col_p2=-1;


//	for( int c1 = 0; c1 < ncomplex; c1++ )
	for( int cx = 0; cx < par_info.nc; cx++ )
	{
		int c1 = par_info.complexes[cx];

		if( ! allComplexes[c1]->isElastic() ) 
			continue;

		// starts at nattach: we only do solution elastic collisions.
		for( int p1 = allComplexes[c1]->nattach; p1 < allComplexes[c1]->nsites; p1++ )
		{
			double sigma1 = allComplexes[c1]->sigma[p1];
	
			int np  = boxing->getNearPts( allComplexes[c1]->rall+3*p1, nearlist, 1.25 * rmin_mult * (use_max_sigma + sigma1) );

			for( int x = 0; x < np; x++ )
			{
				int c2 = complex_for_id[nearlist[x]];
				int p2 = subp_for_id[nearlist[x]];
		
				if( ! allComplexes[c2]->isElastic() ) 
					continue;
			
				if( c1 == c2 ) 
					continue;

				if( ourp[c1] && ourp[c2] && c2 < c1 ) 
					continue;

				if( p2 >=  allComplexes[c1]->nattach )
				{
					// solution-solution collision: only thing currently implemented (1/2019-AJS)	
					// two objects, r1 + v1 with mass m1, r2 + v2 with mass m2
					
					double *r1 = allComplexes[c1]->rall+3*p1;
					double *r2 = allComplexes[c2]->rall+3*p2;
					double m1 = allComplexes[c1]->mass[p1];
					double m2 = allComplexes[c2]->mass[p2];

					double dr[3] = { r1[0]-r2[0], r1[1]-r2[1],r1[2]-r2[2]};
					theSurface->wrapPBC( dr, alphas );
					double dr2 = dr[0]*dr[0]+dr[1]*dr[1]+dr[2]*dr[2];
	
					double p1r = allComplexes[c1]->sigma[p1];
					double p2r = allComplexes[c2]->sigma[p2];
					double r02 = (p1r+p2r)*(p1r+p2r);

					if( dr2 < r02 )
					{

						double v1[3] = { 
							allComplexes[c1]->qdot[3*p1+0], 
							allComplexes[c1]->qdot[3*p1+1], 
							allComplexes[c1]->qdot[3*p1+2] };
						double v2[3] = { 
							allComplexes[c2]->qdot[3*p2+0], 
							allComplexes[c2]->qdot[3*p2+1], 
							allComplexes[c2]->qdot[3*p2+2] }; 
						double dv[3] = { v1[0]-v2[0], v1[1]-v2[1], v1[2]-v2[2] };
						double dv2 = 
							(dv[0])*(dv[0]) +
							(dv[1])*(dv[1]) +
							(dv[2])*(dv[2]);
						if( dv2 > 1e-30 )
						{

							double drdv = dr[0]*dv[0]+dr[1]*dv[1]+dr[2]*dv[2];
		
							if( 4*drdv*drdv-4*dv2*(dr2-r02) >= 0 )		
							{
								// time of collision.
								double t1 = (-0.5*sqrt(4*drdv*drdv-4*dv2*(dr2-r02))-drdv)/dv2;
								double t2 = (+0.5*sqrt(4*drdv*drdv-4*dv2*(dr2-r02))-drdv)/dv2;
								
								int t1_ok = t1 < 0 && t1  > -time_step;
								int t2_ok = t2 < 0 && t2  > -time_step;
		
								double use_t = t1;
								int ok = 0;
								if( t1_ok && t2_ok )
								{
									if( t1 > t2 )		
										use_t = t1;
									else
										use_t = t2;
									ok =1;
								}
								else if( t1_ok )
								{
									ok = 1;
									use_t = t1;
								}
								else if( t2_ok )
								{
									ok = 1;
									use_t = t2;
								}


								if( ok )
								{
									if( -use_t > time_preceding_collision )
									{
										time_preceding_collision = -use_t;

										col_c1 = c1;
										col_p1 = p1;
										col_c2 = c2;
										col_p2 = p2;
									} 
								}
							}
						}
					}
				}				
			}
		}
	}
	
#ifdef PARALLEL
	double min_times[par_info.nprocs];
	double culprits[par_info.nprocs*4];
	int my_id = par_info.my_id;
	min_times[my_id] = time_preceding_collision;
	culprits[4*my_id+0] = col_c1;
	culprits[4*my_id+1] = col_p1;
	culprits[4*my_id+2] = col_c2;
	culprits[4*my_id+3] = col_p2;

	ParallelGather( min_times, 1 );
	ParallelGather( culprits, 4 );

	if( par_info.my_id == BASE_TASK )
	{
		for( int px = 0; px < par_info.nprocs; px++ )
		{
			if( min_times[px] > time_preceding_collision )
			{
				time_preceding_collision = min_times[px];
				col_c1 = culprits[4*px+0];
				col_p1 = culprits[4*px+1];
				col_c2 = culprits[4*px+2];
				col_p2 = culprits[4*px+3];
			}
		}
	}

	ParallelBroadcast( & time_preceding_collision,1 );

	ParallelBroadcast( &col_c1,1 );
	ParallelBroadcast( &col_c2,1 );
	ParallelBroadcast( &col_p1,1 );
	ParallelBroadcast( &col_p2,1 );	 
#endif

	*out_col_c1 = col_c1;
	*out_col_p1 = col_p1;
	*out_col_c2 = col_c2;
	*out_col_p2 = col_p2;

	// find out each processes' min time and who the culprits are.

	free(nearlist);
	free(ourp);
	free(complex_for_id);
	free(subp_for_id);
	
	boxing->clearBoxing(to_clear, ntotp);
	free(to_clear);

	return time_preceding_collision;
}
