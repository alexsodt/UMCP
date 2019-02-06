
void setupSpline( double XMIN, double XMAX, int NBINS, int index, int periodic=0 );
double evaluateSpline( double x, int index);
void SolveSpline( int index );
void AddPointToSpline( double x, double f, int index );
double splined2FdX2( double x, int index );
double splinedFdX( double x, int index );
