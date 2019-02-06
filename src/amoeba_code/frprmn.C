#include <stdio.h>
#include <math.h>
#define NRANSI
#include "nrutil.h"

int nr_iter_max = 500000;


#define EPS 1.0e-14
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []))
{
	void linmin(double p[], double xi[], int n, double *fret,
		double (*func)(double []));
       void dlinmin(double p[], double xi[], int n, double *fret, 
                 double (*func)(double []), void (*dfunc)(double [], double []));
	int j,its;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=nr_iter_max;its++) {
		*iter=its;
//		dlinmin(p,xi,n,fret,func,dfunc);
		linmin(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp= *fret;
		(*dfunc)(p,xi);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
//	printf("Too many iterations in frprmn\n");
//	nrerror("Too many iterations in frprmn");
}
#undef EPS
#undef FREEALL
#undef NRANSI
