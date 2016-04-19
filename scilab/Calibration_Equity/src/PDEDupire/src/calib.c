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

The program CALIB computes the volatility Sigma from given values of 
S_0, r, q and option prices, by minimizing the following cost function: 
  F(Sigma) = G(Sigma)                         + lambda F_1(Sigma)
           = 1/2 || V(Sigma) - V_{tilde} ||^2 + lambda || nabla Sigma ||^2
with : 
- V(Sigma) denotes the computed option prices with PDE Dupire
- V_{tilde} the option prices given by the market
- ||.|| the Euclidean norm
- lambda the coefficient of regularization

The data prices are contained in the file name_in_data which has three 
columns: K_i T_j V_{i,j} with V_{i,j} = V(K_i,T_j)

The minimization algorithm is defined by choice_optim:
- 1 --> Quasi-Newton without constraints (implemented by V. Barette)
- 2 --> Quasi-Newton with bounds (implemented by Nocedal et al.)
In each case, some parameters concerning the optimization can be 
defined by the file name_in_optim.
If choice_optim = 1 (QN), name_in_optim contains:
 gradtol : tolerance on the relative gradient 
 steptol : tolerance on the relative change of x
 verbosity : level of printed information (0 --> 3)
 saveSuccessiveXinFile : save succesive x0 in the file data.out (0 ou 1)
 maxCounter : maximum number of iterations
 lambda : coeff of F1 in the obj function : F = G+lambda*F1
If choice_optim = 2 (QN with bounds), name_in_optim contains:
 pgtol : tolerance on the infinity norm of the projected gradient
 factr : tolerance on the change of the objective function
 iprint : level of printed information
 maxCounter : maximum number of iterations
 sigma_min : minimal bound on sigma
 sigma_max : maximal bound on sigma
 lambda : coeff of F1 in the obj function : F = G+lambda*F1

The initial volatility is given by its bicubic spline degrees of freedom 
on a coarse grid nxm, given by the file name_in_sigmainit_ddl. This file 
contains: 
n, m, y_0, ..., y_n, t_0, ..., t_m, ddl --> (n+3)(m+3) values.

At the end of the calibration algorithm, the degrees of freedom are 
saved in the file name_out_sigmaest_ddl.

In order to visualize the volatilities, the initial and the estimated 
Sigma can also be computed on a rectangular grid given by the file 
name_in_sigma_visu. This file contains:
n_visu, m_visu, S_0, ..., S_{n_visu}, t_0, ..., t_{m_visu}

The values of the initial and the estimated Sigma on the rectangular grid 
are saved in the files name_out_sigmainit_visu and name_out_sigmaest_visu 
respectively. These files contain:
n_visu, m_visu, S_0, ..., S_{n_visu}, t_0, ..., t_{m_visu}, 
Sigma_{0,0}..Sigma_{0,m_visu} ... Sigma_{n_visu,0}..Sigma_{n_visu,m_visu}
where Sigma_{i,j} = Sigma_{S_i,t_j}

The input file CALIB.IN contains the following data:
###################### VARIABLES S_0, R, Q AND OPTIONTYPE ######################
S_0 : price of the underlying asset at t=t_0
r : risk-free rate     
q : dividends (continuously compounded) 
optionType : type of the option (1 for call, 0 for put)  
##################################### DUPIRE ###################################
t_0 : time origin
T_max : maximal maturity (in year)
y_min : --> S_min=exp(y_min)     
y_max : --> S_max=exp(y_max)
N : number of space steps of the fine grid
M : number of time steps of the fine grid
gridType : type of the y_fineGrid (0 for regular, 1 for tanh)
theta : type of scheme used
################################## OPTIMIZATION ################################
choice_optim : choice of the optimization method (1-->QN, 2-->QN with bounds)
################################# FILES IN/OUT #################################
name_in_optim : name of the file containing the optimization parameters
name_in_data : name of the file containing the data prices
name_in_sigmainit_ddl : name of the file containing the degrees of freedom of
  the initial sigma
name_out_sigmaest_ddl : name of the file containing the degrees of freedom of
  the estimated sigma
name_in_sigma_visu : name of the file containing the visualization parameters
  Type "_" if no visualization is required
name_out_sigmainit_visu : name of the out file containing the initial sigma 
  to be visualized (Type "_" is no visualization is required)
name_out_sigmaest_visu : name of the out file containing the estimated sigma
  to be visualized (Type "_" is no visualization is required)
*/

double S_0,r,q;
int optionType;
int N,M;
double theta;

double lambda;

struct marketData *marketOptionPrices;

double *y_fineGrid, *T_fineGrid;
int n,m;
double *y_coarseGrid, *T_coarseGrid;

double costDupire(double *x){

  return F(x,lambda,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid,n,m,N,M,S_0,r,q,theta,optionType,marketOptionPrices);

}

void gradDupire(double *x, double *grad){

return computeGrad_F_finite_diff(grad,x,lambda,n,m,N,M,y_coarseGrid,T_coarseGrid,y_fineGrid,T_fineGrid,r,q,S_0,theta,optionType,marketOptionPrices);
 
}

int calib(){
  
  /* declarations */
  double t_0,T_max,y_min,y_max;
  int gridType, choice_optim;
  char *name_in_optim, *name_in_data;
  char *name_in_sigmainit_ddl, *name_out_sigmaest_ddl;
  char *name_in_sigma_visu, *name_out_sigmainit_visu, *name_out_sigmaest_visu;

  char *name_in_calib;

  FILE *fic_in_sigmainit_ddl;
  FILE *fic_in_sigma_visu;
  FILE *fic_out_sigmainit_visu;

  double gradtol,steptol;
  int maxCounter,verbosity,saveSuccessiveXinFile;

  int nbParamSigma;
  double *sigma_init, *sigma_est;

  int i,j,k,l;

  int Nsigma_visu,Msigma_visu;
  double *Ksigma_visu, *Ysigma_visu, *Tsigma_visu;
  double **sigmainitVisuGrid, **sigmaestVisuGrid;

  int m_nbcorrec, nbwa, iprint;
  double costf, factr, pgtol;
  double sigma_min, sigma_max;
  double *lower_bound, *upper_bound, *gradient;
  double *wa;
  int *nbd, *iwa;
  char task[60], csave[60];
  int lsave[4];
  int isave[44];
  double dsave[29];

  printf("\nCalibration using a PDE Dupire simulator\n\n");

  name_in_calib = "calib.in";
  
  /* initialization of S_0, r, q,..., name_out_sigmaest_visu from the file 
     name_in_calib */
  loadParamCalib(name_in_calib,&S_0,&r,&q,&optionType,&t_0,&T_max,&y_min,&y_max,&N,&M,&gridType,&theta,&choice_optim,&name_in_optim,&name_in_data,&name_in_sigmainit_ddl,&name_out_sigmaest_ddl,&name_in_sigma_visu,&name_out_sigmainit_visu,&name_out_sigmaest_visu);

  /* initialization of data prices from the file name_in_data */
  printf("Data are defined from the file \"%s\"\n\n",name_in_data);
  loadDataPrices(name_in_data,&marketOptionPrices);

  /* initialization of n, m, y_coarseGrid, T_coarseGrid and sigma_init from 
     the file name_in_sigmainit_ddl */
  printf("Sigma_init is defined from its degrees of freedom (see the file \"%s\")\n\n",name_in_sigmainit_ddl);
  loadSigmaddl(name_in_sigmainit_ddl,&n,&m,&y_coarseGrid,&T_coarseGrid,&sigma_init);
  nbParamSigma = (n+3)*(m+3); // i.e (n+1)*(m+1)+2*(m+1)+2*(n+1)+4

  /* y_fineGrid and T_fineGrid are used by costDupire and gradDupire */
  y_fineGrid = (double *) malloc((N+1)*sizeof(double));
  T_fineGrid = (double *) malloc((M+1)*sizeof(double));
  buildFineGrid(y_fineGrid,T_fineGrid,N,M,t_0,T_max,y_min,y_max,S_0,gridType);

  /* memory allocation */
  sigma_est = (double *) malloc(nbParamSigma*sizeof(double));

  if (choice_optim==1) // Quasi-Newton without constraints
    {

      /* initialization of gradtol, steptol, verbosity, saveSuccessiveXinFile, 
	 maxCounter and lambda from the file name_in_optim */
      loadParamOptim(name_in_optim,&gradtol,&steptol,&verbosity,&saveSuccessiveXinFile,&maxCounter,&lambda);

      /* optimization */
      printf("-------------------------------------------");
      printf("-------------------------------------------\n");
      printf("Optimization: computation of Sigma_est using a QuasiNewton algorithm\n");
      printf("-------------------------------------------");
      printf("-------------------------------------------\n");

      QuasiNewton(nbParamSigma,sigma_init,costDupire,gradDupire,sigma_est,gradtol,steptol,maxCounter,verbosity,saveSuccessiveXinFile);

      printf("-------------------------------------------");
      printf("-------------------------------------------\n");

    }
  else if (choice_optim==2) // Quasi-Newton with bounds
    {

      m_nbcorrec = 5; // number of corrections used in the limited memory matrix

      /* initialization of pgtol, factr, iprint, maxCounter, sigma_min, sigma_max
	 and lambda from the file name_in_optim */
      loadParamOptim2(name_in_optim,&pgtol,&factr,&iprint,&maxCounter,&sigma_min,&sigma_max,&lambda);

      /* optimization */
      printf("-------------------------------------------");
      printf("-------------------------------------------\n");
      printf("Optimization: computation of Sigma_est using a QuasiNewton algorithm with bounds\n");
      printf("-------------------------------------------");
      printf("-------------------------------------------\n");

      /* memory allocation */
      lower_bound = (double *) malloc(nbParamSigma*sizeof(double));  
      upper_bound = (double *) malloc(nbParamSigma*sizeof(double));  
      nbd = (int *) malloc(nbParamSigma*sizeof(int));  
      gradient = (double *) malloc(nbParamSigma*sizeof(double));  
      nbwa = 2*m_nbcorrec*nbParamSigma+4*nbParamSigma+12*m_nbcorrec*m_nbcorrec+12*m_nbcorrec;
      wa = (double *) malloc(nbwa*sizeof(double));  
      iwa = (int *) malloc(3*nbParamSigma*sizeof(int));  

      for (i=0;i<nbParamSigma;i++)
	{
	  sigma_est[i] = sigma_init[i];
	  lower_bound[i] = sigma_min;
	  upper_bound[i] = sigma_max;
	  nbd[i] = 0; // nbd[i]=0 if sigma_est[i] is unbounded
	}

      // bounds for the values of sigma on the nodes (not for the derivatives)
      for (i=0;i<(n+1)*(m+1);i++)
	nbd[i] = 2; // nbd[i]=2 if sigma_est[i] has both lower and upper bounds

      inittask_(task); // task = 'START'
      
      setulb_(&nbParamSigma,&m_nbcorrec,sigma_est,lower_bound,upper_bound,nbd,&costf,gradient,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave);
      
      while (memcmp(task,"FG",2) == 0 || memcmp(task,"NEW_X",5) == 0)
	{
	  if (memcmp(task,"FG",2) == 0) // si les 2 premiers caracteres de task sont FG
	    {
	      costf = costDupire(sigma_est);
	      gradDupire(sigma_est,gradient);
	    }
	  if (memcmp(task,"NEW_X",5) == 0 && isave[29]>=maxCounter)
	    //isave[29] en C ou isave(30) en Fortran correspond au nb d'iter. courant
	    {
	      inittask2_(task); // task = 'STOP: Maximum number of iterations reached'
	    }
	  setulb_(&nbParamSigma,&m_nbcorrec,sigma_est,lower_bound,upper_bound,nbd,&costf,gradient,&factr,&pgtol,wa,iwa,task,&iprint,csave,lsave,isave,dsave);
	}

      printf("-------------------------------------------");
      printf("-------------------------------------------\n");
      
    }
  
  /* save the degrees of freedom of sigmaest in the file name_out_sigmaest_ddl */
  if (strcmp(name_out_sigmaest_ddl,"_") != 0){
    saveSigmaddl(name_out_sigmaest_ddl,sigma_est,n,m,y_coarseGrid,T_coarseGrid);
    printf("--> Sigma_est saved in the file \"%s\"\n\n",name_out_sigmaest_ddl);
  }

  if (strcmp(name_in_sigma_visu,"_") != 0){

    if (strcmp(name_out_sigmainit_visu,"_") != 0 || strcmp(name_out_sigmaest_visu,"_") != 0){
      loadSigmaGrid(name_in_sigma_visu,&Nsigma_visu,&Msigma_visu,&Ksigma_visu,&Tsigma_visu);
      Ysigma_visu = (double *) malloc((Nsigma_visu+1)*sizeof(double));
      for (i=0;i<Nsigma_visu+1;i++)
	Ysigma_visu[i] = log(Ksigma_visu[i]);
    }

    if (strcmp(name_out_sigmainit_visu,"_") != 0){
      printf("Computation of Sigma_init on the visualization grid defined in the file \"%s\"\n",name_in_sigma_visu);
      sigmainitVisuGrid = (double **) malloc((Nsigma_visu+1)*sizeof(double *));
      for (i=0;i<Nsigma_visu+1;i++)
	sigmainitVisuGrid[i] = (double *) malloc((Msigma_visu+1)*sizeof(double));
      sigmaddlToSigmaFineGrid(sigmainitVisuGrid,n,m,y_coarseGrid,T_coarseGrid,sigma_init,Ysigma_visu,Tsigma_visu,Nsigma_visu,Msigma_visu);
      savePriceOrSigma(name_out_sigmainit_visu,sigmainitVisuGrid,Nsigma_visu,Msigma_visu,Ksigma_visu,Tsigma_visu);
      printf("--> Sigma_init saved in the file \"%s\"\n\n",name_out_sigmainit_visu);
    }

    if (strcmp(name_out_sigmaest_visu,"_") != 0){
      printf("Computation of Sigma_est on the visualization grid defined in the file \"%s\"\n",name_in_sigma_visu);
      sigmaestVisuGrid = (double **) malloc((Nsigma_visu+1)*sizeof(double *));
      for (i=0;i<Nsigma_visu+1;i++)
	sigmaestVisuGrid[i] = (double *) malloc((Msigma_visu+1)*sizeof(double));
      sigmaddlToSigmaFineGrid(sigmaestVisuGrid,n,m,y_coarseGrid,T_coarseGrid,sigma_est,Ysigma_visu,Tsigma_visu,Nsigma_visu,Msigma_visu);
      savePriceOrSigma(name_out_sigmaest_visu,sigmaestVisuGrid,Nsigma_visu,Msigma_visu,Ksigma_visu,Tsigma_visu);
      printf("--> Sigma_est saved in the file \"%s\"\n\n",name_out_sigmaest_visu);
    }

  }

  /* free memory */
  free(sigma_init);
  free(sigma_est);
  free(y_coarseGrid);
  free(T_coarseGrid);
  free(y_fineGrid);
  free(T_fineGrid);
  return 0;
  
}
