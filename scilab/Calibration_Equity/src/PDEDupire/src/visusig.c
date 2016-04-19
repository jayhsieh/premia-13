#include <stdio.h>
#include <math.h>
#include <time.h>

#include "inout.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.

The program VISUSIG computes the volatility Sigma on a rectangular 
grid from its bicubic spline degrees of freedom. It is useful to 
visualize Sigma.

The volatility on the grid nxm is given by the file name_in_sigma_ddl. 
This file contains: 
n, m, y_0, ..., y_n, t_0, ..., t_m, ddl --> (n+3)(m+3) values.

The rectangular grid is given by the file name_in_sigma_visu containing: 
n_visu, m_visu, S_0, ..., S_{n_visu}, t_0, ..., t_{m_visu}.

The values of Sigma on the rectangular grid are saved in the file 
name_out_sigma_visu: 
n_visu, m_visu, S_0, ..., S_{n_visu}, t_0, ..., t_{m_visu}, 
Sigma_{0,0}..Sigma_{0,m_visu} ... Sigma_{n_visu,0}..Sigma_{n_visu,m_visu}
where Sigma_{i,j} = Sigma_{S_i,t_j}

The input file VISUSIG.IN contains the following data:
################################# FILES IN/OUT #################################
name_in_sigma_ddl : name of the file containing the ddl
name_in_sigma_visu : name of the file containing the visualization parameters
name_out_sigma_visu : name of the out file containing sigma 
*/

main(){ 
  
  /* declarations */
  int n, m;
  double *Y_coarseGrid, *T_coarseGrid;
  double *sigma;

  char *name_in_sigma_ddl, *name_in_sigma_visu, *name_out_sigma_visu;
  
  char *name_in_visusig;

  int i,j,k,l;

  int Nsigma_visu,Msigma_visu;
  double *Ysigma_visu, *Ksigma_visu, *Tsigma_visu;
  double **sigmaVisuGrid;

  name_in_visusig = "visusig.in";

  loadParamVisusig(name_in_visusig,&name_in_sigma_ddl,&name_in_sigma_visu,&name_out_sigma_visu);

  /* initialization of n, m, Y_coarseGrid, T_coarseGrid and sigma from 
     the file name_in_sigma_ddl */
  loadSigmaddl(name_in_sigma_ddl,&n,&m,&Y_coarseGrid,&T_coarseGrid,&sigma);

  loadSigmaGrid(name_in_sigma_visu,&Nsigma_visu,&Msigma_visu,&Ksigma_visu,&Tsigma_visu);

  Ysigma_visu = (double *) malloc((Nsigma_visu+1)*sizeof(double));
  for (i=0;i<Nsigma_visu+1;i++)
    Ysigma_visu[i] = log(Ksigma_visu[i]);

  sigmaVisuGrid = (double **) malloc((Nsigma_visu+1)*sizeof(double *));
  for (i=0;i<Nsigma_visu+1;i++)
    sigmaVisuGrid[i] = (double *) malloc((Msigma_visu+1)*sizeof(double));

  sigmaddlToSigmaFineGrid(sigmaVisuGrid,n,m,Y_coarseGrid,T_coarseGrid,sigma,Ysigma_visu,Tsigma_visu,Nsigma_visu,Msigma_visu);
  
  savePriceOrSigma(name_out_sigma_visu,sigmaVisuGrid,Nsigma_visu,Msigma_visu,Ksigma_visu,Tsigma_visu);
  
  printf("\n--> Sigma saved in the file \"%s\"\n\n",name_out_sigma_visu);

  /* free memory */
  free(sigma);
  free(Y_coarseGrid);
  free(T_coarseGrid);
  free(Ksigma_visu);
  free(Tsigma_visu);

}
