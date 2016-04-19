
#include "pnl/pnl_complex.h"
#include "nrutil.h"


dcomplex cfrncall(int model,double  rf,double dt,dcomplex g,double aa,double parameters[]);

dcomplex cfrn(int model,double rf, double dt,dcomplex g, double parameters[]);

dcomplex cfCDF(int model,double dt, dcomplex g, double aa, double parameters[]);

dcomplex cfrnshifted(int model,double  aa,double rf,double dt,dcomplex g, double parameters[]);

dcomplex cfLevy(int model, double dt, dcomplex g, double parameters[]);

dcomplex cfGauss(double sg, dcomplex g);

double  MomentsLevy(int model, double  rf,int moment, double dt, double parameters[]);

dcomplex cfrnstandardized(int model,double rf, double dt, dcomplex g, double parameters[]);

//NIG
dcomplex cfNig(double alpha, double beta,double delta, dcomplex g);

//'meixner
dcomplex cfMeixner(double alpha, double beta,double delta,dcomplex g);


dcomplex cfVarianceGamma(double sg, double nu,double theta,dcomplex g);


    
//'cgmy
dcomplex cfCgmy(double ccc, double ggg,double mmm, double yyy,dcomplex g);




//'de

dcomplex cfDe(double sg, double lambda, double p, double eta1, double eta2, dcomplex g);


///'jd
dcomplex cfMerton(double sg, double alpha, double lambda, double delta, dcomplex g);

///Compute tail bounds using moment x is in log terms
double  BoundUpperTailLevy(int model, double x, double  rf, double dt, int maxmoment, double parameters[]);

double  BoundLowerTailLevy(int model, double x, double  rf, double dt, int maxmoment, double parameters[]);

