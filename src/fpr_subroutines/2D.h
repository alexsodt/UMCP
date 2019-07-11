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

#endif
