
#include "pnl/pnl_complex.h"
#include "nrutil.h"


int fft(long n, double **cg, int iflg);

int dft(long n, double cg[][2], double cf[][2], int iflg);

void TableIFT(int choice, int model, double rf, double dt, long N, 
                   double dx, double aa, double parameters[], 
				   double inv[], double logk[]);


void TableIFRT(int choice, int model, double rf,double  dt,long N, 
                   double b, double aa, double parameters[],
				   double inv[], double logk[]);

double find_umax_cf(int model,double  rf,double dt, double aa, double parameters[]);

void TableDensity(int choice, int model, double rf,double  dt, long N, 
                   double b, double parameters[],
				   double dens[], double logk[]);

