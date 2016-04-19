#include <stdio.h>
#include <math.h>
#include <time.h>
#include "objFunction.h"
#include "gradFiniteDiff.h"
#include "spline.h"
#include "inout.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.

The program RAFSIG computes the bicubic spline degrees of freedom 
of the volatility Sigma on a grid n2xm2, from its degrees of freedom 
on a coarser grid n1xm1.

The volatility on the grid n1xm1 is given by the file name_in_sigma1_ddl. 
This file contains: 
n1, m1, y_0, ..., y_n1, t_0, ..., t_m1, ddl --> (n1+3)(m1+3) values.

The grid n2xm2 is obtained by splitting the grid n1xm1 in nbsplit_y 
subdivisions in y and nbsplit_T subdivisions in T. The subdivisions 
are regular on the grid [K,T]:
 Grid1 [Y,T] --> Grid1 [K,T] splitted regularly w.r.t nbsplit_y and nbsplit_T
 --> Grid2 [K,T] --> Grid2 [Y,T]

The volatility on the grid n2xm2 is saved in the file name_out_sigma2_ddl:
n2, m2, y_0, ..., y_n2, t_0, ..., t_m2, ddl --> (n2+3)(m2+3) values.

The input file RAFSIG.IN contains the following data:
################################# FILES IN/OUT #################################
name_in_sigma1_ddl : name of the file containing the ddl
nbsplit_y : number (integer) of subdivisions in y
nbsplit_T : number (integer) of subdivisions in T
name_out_sigma2_ddl : name of the out file containing the new ddl
*/

main(){ 
  
  /* declarations */
  double *y_fineGrid1, *T_fineGrid1, *y_fineGrid2, *T_fineGrid2;
  int n1, m1, n2, m2;
  double *y_coarseGrid1, *T_coarseGrid1, *y_coarseGrid2, *T_coarseGrid2;
  double *sigma1, *sigma2;

  int nbsplit_y,nbsplit_T;
  char *nbsplit_y_char,*nbsplit_T_char;

  char *name_in_sigma1_ddl, *name_out_sigma2_ddl;

  char *name_in_rafsig;

  int i,j,k,l;

  name_in_rafsig = "rafsig.in";

  loadParamRafsig(name_in_rafsig,&name_in_sigma1_ddl,&nbsplit_y_char,&nbsplit_T_char,&name_out_sigma2_ddl);

  /* initialization of n1, m1, y_coarseGrid1, T_coarseGrid1 and sigma1 from 
     the file name_in_sigma1_ddl */
  /* nbParamSigma = (n+3)*(m+3) */
 
  loadSigmaddl(name_in_sigma1_ddl,&n1,&m1,&y_coarseGrid1,&T_coarseGrid1,&sigma1);
 
  nbsplit_y = stringToInt(nbsplit_y_char);
  nbsplit_T = stringToInt(nbsplit_T_char);
  n2 = nbsplit_y*n1;
  m2 = nbsplit_T*m1;
  y_coarseGrid2 = (double *) malloc((n2+1)*sizeof(double));
  for (i=0;i<n1;i++)
    for (j=i*nbsplit_y;j<(i+1)*nbsplit_y;j++)
      y_coarseGrid2[j] = log(exp(y_coarseGrid1[i]) + (j-nbsplit_y*i)*(exp(y_coarseGrid1[i+1])-exp(y_coarseGrid1[i]))/nbsplit_y);
  y_coarseGrid2[n2] =y_coarseGrid1[n1]; 

  T_coarseGrid2 = (double *) malloc((m2+1)*sizeof(double));
  for (i=0;i<m1;i++)
    for (j=i*nbsplit_T;j<(i+1)*nbsplit_T;j++)
      T_coarseGrid2[j] = T_coarseGrid1[i] + (j-nbsplit_T*i)*(T_coarseGrid1[i+1]-T_coarseGrid1[i])/nbsplit_T;
  T_coarseGrid2[m2] =T_coarseGrid1[m1]; 
  
  ddlGrid1ToGrid2(&sigma2,sigma1,n1,m1,n2,m2,y_coarseGrid1,T_coarseGrid1,y_coarseGrid2,T_coarseGrid2);  
  
  saveSigmaddl(name_out_sigma2_ddl,sigma2,n2,m2,y_coarseGrid2,T_coarseGrid2);
  
  /* free memory */
  free(sigma1);
  free(sigma2);
  free(y_coarseGrid1);
  free(T_coarseGrid1);
  free(y_coarseGrid2);
  free(T_coarseGrid2);

}
