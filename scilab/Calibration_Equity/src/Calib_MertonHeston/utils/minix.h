#include "costFunction.h"
#include "gradFunction.h"

double costFunction1(int dimx,double *x);
void gradFunction1(int dimx,double *x,double *grad);
void x2minix(double *x,double *minix);
void minix2x(double *x,double *minix);
int initminix(int dimx,int *nbd,double *xmin,double *xmax,int *mininbd,double *minixmin,double *minixmax);
void FreeMinix(void);
