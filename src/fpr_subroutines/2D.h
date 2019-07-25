#ifndef __2Dh__
#define __2Dh__

typedef struct
{
	int nEntries;
	int nEntriesSpace;
	int DDtableindex;
	double *TBLID;
	gsl_matrix **contsur;
	gsl_matrix **contnorm;
	gsl_matrix **contpir;
	int *vecLen;
}  twoDReactionPTable;

twoDReactionPTable *newTable( void );
double get_2D_2D_rxn_prob( double R1, double kr, double bindrad, double Dtot, double deltat, double Rmax );

#endif
