#ifndef  _LINSYS_H
#define _LINSYS_H

int gmres(double **H,double **band,double x[],double b[],int m,int precond,int n,int max_iter,double tol,double *pivots);
int bicgstab(double **band,double x[],double b[],int n,int max_iter,double tol,int precond,double *pivots);
void  Diagonal_Precond(double **band, int n,double *pivots);
void  ILU_Precond(double **band, int n,double *pivots);
void mv_product(int n,double **band,double b[],double x[]);
#endif
