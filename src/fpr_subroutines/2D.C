#include <stdio.h>
#include <cmath>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_vector.h>
#include "2Drelated.h"
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <math.h>
#include <cfloat>
#include "2D.h"

// DDtableindex: index of where we store information about a 2D reaction.

static double tolKvalue = 1e-3; // correct this. depends on units, machine precision etc.
static double tolDvalue = 1e-3;
static double rtol      = 1e-10;
//
//  copied and reorganized slightly from Margaret Johnson's code.
//

twoDReactionPTable *newTable( void )
{
	twoDReactionPTable *theTable = (twoDReactionPTable *)malloc( sizeof(twoDReactionPTable) );

	theTable->nEntries = 0;
	theTable->nEntriesSpace = 1;
	theTable->TBLID = (double *)malloc(  sizeof(double) * 2 * theTable->nEntriesSpace );
	theTable->vecLen = (int *)malloc( sizeof(int) * theTable->nEntriesSpace );
	theTable->contsur = (gsl_matrix **)malloc(  sizeof(gsl_matrix *) * theTable->nEntriesSpace );
	theTable->contnorm = (gsl_matrix **)malloc( sizeof(gsl_matrix *) * theTable->nEntriesSpace );
	theTable->contpir = (gsl_matrix **)malloc( sizeof(gsl_matrix *) * theTable->nEntriesSpace );

	return theTable;
}

static twoDReactionPTable *rxnTable = NULL;

double get_2D_2D_rxn_prob( double R1, double kr, double bindrad, double Dtot, double deltat, double Rmax, double prev_sep, double ps_prev, double *p0_ratio )
{
	if( !rxnTable )
		rxnTable = newTable();

	double ktemp = kr / 2 / bindrad;
	int uniquetableindex = -1;
	int tableexistflag = 0;

	for( int idd = 0; idd < rxnTable->nEntries; idd++ )
	{
		if( fabs(rxnTable->TBLID[2*idd+0]-ktemp) < tolKvalue && fabs(rxnTable->TBLID[2*idd+1]-Dtot) < tolDvalue )
		{
			tableexistflag = 1;
			uniquetableindex = idd;
			break;
		}
	}
	
	int vL=0;

	if( uniquetableindex == -1 )
	{
		// compute table.

		if( rxnTable->nEntries == rxnTable->nEntriesSpace )
		{
			rxnTable->TBLID = (double *)realloc( rxnTable->TBLID, sizeof(double) * 2 * rxnTable->nEntriesSpace );
			rxnTable->vecLen = (int *)realloc( rxnTable->vecLen, sizeof(int) * rxnTable->nEntriesSpace );
			rxnTable->contsur = (gsl_matrix **)realloc( rxnTable->contsur, sizeof(gsl_matrix *) * rxnTable->nEntriesSpace );
			rxnTable->contnorm = (gsl_matrix **)realloc( rxnTable->contnorm, sizeof(gsl_matrix *) * rxnTable->nEntriesSpace );
			rxnTable->contpir = (gsl_matrix **)realloc( rxnTable->contpir, sizeof(gsl_matrix *) * rxnTable->nEntriesSpace );
		} 

		uniquetableindex = rxnTable->nEntries;

		rxnTable->TBLID[2*uniquetableindex+0] = ktemp;
		rxnTable->TBLID[2*uniquetableindex+1] = Dtot;

		rxnTable->vecLen[uniquetableindex] = sizelookup(bindrad, Dtot, deltat, Rmax);

		vL = rxnTable->vecLen[uniquetableindex];
		rxnTable->contsur[uniquetableindex] = gsl_matrix_alloc(2,vL); 
		rxnTable->contnorm[uniquetableindex] = gsl_matrix_alloc(2,vL); 
		rxnTable->contpir[uniquetableindex] = gsl_matrix_alloc(vL,vL); 
		      
		DDmatrixcreate(rxnTable->contsur[uniquetableindex], rxnTable->contnorm[uniquetableindex], rxnTable->contpir[uniquetableindex], bindrad, Dtot, ktemp, deltat, Rmax);

		rxnTable->nEntries++;
	}
		    
	double probvec1 = DDpsur(rxnTable->contsur[uniquetableindex], Dtot, deltat, R1, bindrad);


	if( prev_sep > 0 )
	{
		*p0_ratio = DDpirr_pfree_ratio_ps(
				rxnTable->contpir[uniquetableindex], rxnTable->contsur[uniquetableindex], rxnTable->contnorm[uniquetableindex], 
				R1, Dtot, deltat, prev_sep, ps_prev, rtol, bindrad);
		if( *p0_ratio < 0 || *p0_ratio > 1e100 )		
		{
			printf("funny p0_ratio.\n");
		}
	}
	else
		*p0_ratio = 1.0;
	return probvec1;
}


