
double splevl(double xb, long n, double x[], double f[], double **c,
		double *dfb, double *ddfb, int *ier);

int spline(double x[], double f[], long n, double **c);

int smooth(long ntab, double x[], double f[], double **c, int np,
	double xp[], double fp[]);

int smoothmod(long ntab, double x[], double f[], double **c, int np,
	double xp[], double fp[]);


double smoothscalar(long ntab, double x[], double f[], double **c, double xp);

