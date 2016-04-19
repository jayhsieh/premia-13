#ifndef ESM_FUNC
#define ESM_FUNC

void ESM_update_const_char(double Kappa,double sigma, double delta, double d);

void Moments_ESM(double Vs, double Vt, double Kappa, double sigma, double delta, double d, double *mean, double *variance);

void values_all_ESM(int M,double Vs, double Vt,double Kappa,double sigma,double delta, double d, double epsilon, double  h, int * N, double * values);


void values_all_AESM(int M,double Vmax, int NS ,double Kappa,double sigma,double delta, double d,double epsilon,
					 double * mean, double * variance, double * h, int * N, double ** values);

double inverse_ESM(double u, double h, int N, double * val);

#endif
