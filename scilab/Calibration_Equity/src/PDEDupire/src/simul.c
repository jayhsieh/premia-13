#include <stdio.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include "DupirePDE.h"
#include "spline.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.

The program SIMUL computes option prices from given S_0, r, q and 
volatility (Sigma) parameters :
- using PDE Dupire ;
- using Black-Scholes formulae with Sigma = Constant ;
- using Black-Scholes formulae with Sigma = Sigma(t).

If Dupire is used, Sigma can be given by the function sigma_func 
(defined in spline.c) or by its bicubic spline degrees of freedom 
on a coarse grid nxm, given by the file name_ddl. This file contains: 
n, m, y_0, ..., y_n, t_0, ..., t_m, ddl --> (n+3)(m+3) values.
If Black-Scholes (with Sigma = Constant) is used, Sigma is given by sigmaCte.
If Black-Scholes (with Sigma = Sigma(t)) is used, the function intSigma 
(defined in DupirePDE.c) has to be defined:
  intSigma = 1/(T-t_0) int_{t_0}^{T} Sigma^2

The option prices can be computed:
- on a rectangular grid given by the file name_in_visu. This file contains:
n_visu, m_visu, K_0, ..., K_{n_visu}, t_0, ..., t_{m_visu}.
In this case, the prices are saved in the file name_out_visu:
n_visu, m_visu, K_0, ..., K_{n_visu}, t_0, ..., t_{m_visu}, 
V_{0,0}..V_{0,m_visu} ... V_{n_visu,0}..V_{n_visu,m_visu} 
where V_{i,j} = V(K_i,T_j)
- for strikes-maturities (K_i,T_j) given by the file name_in_data.
In this case, the prices are saved in the file name_out_data containing 
three columns: K_i T_j V_{i,j}

The input file SIMUL.IN contains the following data:
###################### VARIABLES S_0, R, Q AND OPTIONTYPE ######################
S_0 : price of the underlying asset at t=t_0
r : risk-free rate     
q : dividends (continuously compounded) 
optionType : type of the option (1 for call, 0 for put)  
########################### DUPIRE OR BLACK-SCHOLES ? ##########################
optionSimul :
  1 for PDE Dupire                --> t_0, T_max, y_min, y_max, N, M, gridType
                                  --> theta, name_ddl or sigma_func (spline.c)
  2 for Black-Scholes (sigma cte) --> sigmaCte
  3 for Black-Scholes (sigma(t))  --> intSigma (DupirePDE.c)
########################### IF OPTIONSIMUL=1 (DUPIRE) ##########################
t_0 : time origin
T_max : maximal maturity (in year)
y_min : --> S_min=exp(y_min)     
y_max : --> S_max=exp(y_max)
N : number of space steps of the fine grid
M : number of time steps of the fine grid
gridType : type of the y_fineGrid (0 for regular, 1 for tanh)
theta : type of scheme used
name_ddl : name of the file containing the degrees of freedom of sigma
  Type "_" if sigma has to be defined from the function sigma_func (spline.c)
############# IF OPTIONSIMUL=2 (B-S WITH SIGMA CONSTANT) #######################
sigmaCte : value of sigma
################################# FILES IN/OUT #################################
name_in_visu : name of the file containing the visualization parameters
  Type "_" if no visualization is required
name_out_visu : name of the out file containing the prices to be visualized
  used if name_in_visu is not given by "_"
name_in_data : name of the file containing Ki and Tj
  Type "_" if no data have to be computed
name_out_data : name of the out file containing Ki, Tj and the prices Vij
  used if name_in_data is not given by "_"
*/

main(){

  /* declarations */
  double S_0,r,q;
  int optionType,optionSimul;
  double t_0,T_max,y_min,y_max;
  int N,M,gridType;
  double theta;
  char *name_in_sigma_ddl;
  double sigmaCte;
  char *name_in_visu, *name_out_visu;
  char *name_in_data, *name_out_data;

  char *name_in_simul;

  double *y_fineGrid, *t_fineGrid; /* grids of the uniformly discretized 
				      values of y and t used for the 
				      finite-difference fine grid */
  double *k_fineGrid;
  double **sigma_fineGrid; /* sigma on the finite-difference fine grid */
  double **price_fineGrid; /* prices of the option today on the 
			      finite-difference fine grid */

  int n,m;
  double *y_coarseGrid, *T_coarseGrid;
  double *sigma_param;

  int i,j,indi,indj;
  double K,T,V,price,price1,price2;

  FILE *fic_in_visu;
  int Nprice_visu, Mprice_visu;
  double *Kprice_visu, *Tprice_visu;
  double **price_visuGrid; /* prices of the option today on the 
			      visualization grid */
  
  FILE *fic_in_data, *fic_out_data;

  name_in_simul = "simul.in";
  
  /* initialization of S_0, r, q,..., name_out_data from the file name_in_simul */
  loadParamSimul(name_in_simul,&S_0,&r,&q,&optionType,&optionSimul,&t_0,&T_max,&y_min,&y_max,&N,&M,&gridType,&theta,&name_in_sigma_ddl,&sigmaCte,&name_in_visu,&name_out_visu,&name_in_data,&name_out_data);

  /* errors */
  if (strcmp(name_in_visu,"_") == 0 && strcmp(name_in_data,"_") == 0){
    printf("\nNo computation of prices because name_in_visu = \"_\" and name_in_data = \"_\" (see the input file \"%s\")\n\n",name_in_simul);
    return;
  }
  if (strcmp(name_in_visu,"_") != 0 && strcmp(name_out_visu,"_") == 0){
    printf("\nThe file name_out_visu has to be initialized in the input file \"%s\"\n\n",name_in_simul);
    return;
  }
  if (strcmp(name_in_data,"_") != 0 && strcmp(name_out_data,"_") == 0){
    printf("\nThe file name_out_data has to be initialized in the input file \"%s\"\n\n",name_in_simul);
    return;
  }
  if (optionSimul != 1 && optionSimul != 2 && optionSimul != 3){
    printf("\nBad value of optionSimul in the input file \"%s\" (type 1, 2 or 3)\n\n",name_in_simul);
    return;
  }

  if (optionSimul == 1)
    printf("\nComputation of option prices using PDE Dupire\n\n");
  else if (optionSimul == 2){
    printf("\nComputation of option prices using Black-Scholes formulae\n\n");
    printf("Sigma = Constant = %lf\n\n",sigmaCte);
  }
  else if (optionSimul == 3){
    printf("\nComputation of option prices using Black-Scholes formulae\n\n");
    printf("Sigma = Sigma(t) is defined in DupirePDE.c via the function intSigma\n\n");
  }

  /* if optionSimul=1 (PDE Dupire), computation of prices on the 
     finite-difference fine grid */
  if (optionSimul == 1){ // PDE Dupire
    /* memory allocation */
    y_fineGrid = (double *) malloc((N+1)*sizeof(double));
    k_fineGrid = (double *) malloc((N+1)*sizeof(double));
    t_fineGrid = (double *) malloc((M+1)*sizeof(double));
    sigma_fineGrid = (double **) malloc((N+1)*sizeof(double *));
    for (i=0;i<=N;i++)
      sigma_fineGrid[i] = (double *) malloc((M+1)*sizeof(double));
    price_fineGrid = (double **) malloc((N+1)*sizeof(double *));
    for (i=0;i<=N;i++)
      price_fineGrid[i] = (double *) malloc((M+1)*sizeof(double));
    /* initialization of y_fineGrid and t_fineGrid */
    buildFineGrid(y_fineGrid,t_fineGrid,N,M,t_0,T_max,y_min,y_max,S_0,gridType);
    for (i=0;i<=N;i++)
      k_fineGrid[i] = exp(y_fineGrid[i]);
    /* initialization of sigma_fineGrid ... */
    if (strcmp(name_in_sigma_ddl,"_") == 0){
      /* ... from a discretization of sigma_func (defined in spline.c) */
      printf("Sigma is defined from a discretization of the function sigma_func defined in spline.c\n\n");
      discretizeSigma(sigma_func,sigma_fineGrid,N,M,y_fineGrid,t_fineGrid);
    } else {
      /* ... from the file name_in_sigma_ddl */
      printf("Sigma is defined from its degrees of freedom (see the file \"%s\")\n\n",name_in_sigma_ddl);
      loadSigmaddl(name_in_sigma_ddl,&n,&m,&y_coarseGrid,&T_coarseGrid,&sigma_param);
      sigmaddlToSigmaFineGrid(sigma_fineGrid,n,m,y_coarseGrid,T_coarseGrid,sigma_param,y_fineGrid,t_fineGrid,N,M);
    }
    /* numerical resolution of Dupire PDE : computation of price_fineGrid[i][j] 
       on the finite-difference fine grid (i=0...N, j=0...M) */
    printf("Computation of prices on the finite-difference fine grid\n\n");
    solve(optionType,price_fineGrid,S_0,N,M,r,q,theta,f,sigma_fineGrid,y_fineGrid,t_fineGrid);
  }
  
  if (strcmp(name_in_visu,"_") != 0){
    printf("Computation of prices on the visualization grid defined in the file \"%s\"\n",name_in_visu);
    /* initialization of Nprice_visu, Mprice_visu, Kprice_visu[0...Nprice_visu] 
       and Tprice_visu[0...Mprice_visu] from the file name_in_visu */
    fic_in_visu = fopen(name_in_visu,"r");
    fscanf(fic_in_visu,"%d",&Nprice_visu);
    fscanf(fic_in_visu,"%d",&Mprice_visu);
    Kprice_visu = (double *)malloc((Nprice_visu+1)*sizeof(double));
    Tprice_visu = (double *)malloc((Mprice_visu+1)*sizeof(double));
    for (i=0;i<Nprice_visu+1;i++)
      fscanf(fic_in_visu,"%lf",&(Kprice_visu[i]));
    for (i=0;i<Mprice_visu+1;i++)
      fscanf(fic_in_visu,"%lf",&(Tprice_visu[i]));
    /* memory allocation */
    price_visuGrid = (double **) malloc((Nprice_visu+1)*sizeof(double *));
    for (i=0;i<=Nprice_visu;i++)
      price_visuGrid[i] = (double *) malloc((Mprice_visu+1)*sizeof(double));
    /* computation of prices on the visualization grid */
    if (optionSimul == 1){ // PDE Dupire
      /* bilinear interpolation : computation of price_visuGrid[i][j] 
	 (i=0...Nprice_visu, j=0...Mprice_visu) from price_fineGrid[i][j] 
	 (i=0...N, j=0...M) */
      indj = 0;
      for (j=0;j<=Mprice_visu;j++){
	indj = find_index(Tprice_visu[j],t_fineGrid,M,indj);
	indi = 0;
	for (i=0;i<=Nprice_visu;i++){
	  indi = find_index(log(Kprice_visu[i]),y_fineGrid,N,indi);
	  price1 = price_fineGrid[indi+1][indj] - (y_fineGrid[indi+1]-log(Kprice_visu[i]))*(price_fineGrid[indi+1][indj]-price_fineGrid[indi][indj])/(y_fineGrid[indi+1]-y_fineGrid[indi]);
	  price2 = price_fineGrid[indi+1][indj+1] - (y_fineGrid[indi+1]-log(Kprice_visu[i]))*(price_fineGrid[indi+1][indj+1]-price_fineGrid[indi][indj+1])/(y_fineGrid[indi+1]-y_fineGrid[indi]);
	  price_visuGrid[i][j] = price2 - (t_fineGrid[indj+1]-Tprice_visu[j])*(price2-price1)/(t_fineGrid[indj+1]-t_fineGrid[indj]);
	}
      }
    } else if (optionSimul == 2){ // Black-Scholes (sigma cte)
      for (i=0;i<=Nprice_visu;i++)
	for (j=0;j<=Mprice_visu;j++)
	  price_visuGrid[i][j] = BSprice(optionType,S_0,t_0,sigmaCte,r,q,Tprice_visu[j],Kprice_visu[i]);
    } else if (optionSimul == 3){ // Black-Scholes (sigma(t))
      for (i=0;i<=Nprice_visu;i++)
	for (j=0;j<=Mprice_visu;j++)
	  price_visuGrid[i][j] = BSpriceTimeDepVol(S_0,t_0,r,q,Tprice_visu[j],Kprice_visu[i],intSigma);
    }
    /* save Nprice_visu, Mprice_visu, Kprice_visu[i], Tprice_visu[j] and 
       price_visuGrid[i][j] in the file name_out_visu */
    savePriceOrSigma(name_out_visu,price_visuGrid,Nprice_visu,Mprice_visu,Kprice_visu,Tprice_visu);
    printf("--> prices saved in the file \"%s\"\n\n",name_out_visu);
  }
  
  if (strcmp(name_in_data,"_") != 0){
    printf("Computation of prices for given strikes-maturities (Ki,Tj) defined in the file \"%s\"\n",name_in_data);
    fic_in_data = fopen(name_in_data,"r");
    fic_out_data = fopen(name_out_data,"w");
    while (fscanf(fic_in_data,"%lf %lf ",&K,&T) == 2){
      if (optionSimul == 1){ // PDE Dupire
	indj = find_index(T,t_fineGrid,M,0);
	indi = find_index(log(K),y_fineGrid,N,0);
	/* bilinear interpolation */
	price1 = price_fineGrid[indi+1][indj] - (y_fineGrid[indi+1]-log(K))*(price_fineGrid[indi+1][indj]-price_fineGrid[indi][indj])/(y_fineGrid[indi+1]-y_fineGrid[indi]);
	price2 = price_fineGrid[indi+1][indj+1] - (y_fineGrid[indi+1]-log(K))*(price_fineGrid[indi+1][indj+1]-price_fineGrid[indi][indj+1])/(y_fineGrid[indi+1]-y_fineGrid[indi]);
	price = price2 - (t_fineGrid[indj+1]-T)*(price2-price1)/(t_fineGrid[indj+1]-t_fineGrid[indj]);
      } else if (optionSimul == 2){ // Black-Scholes (sigma cte)
	price = BSprice(optionType,S_0,t_0,sigmaCte,r,q,T,K);
      } else if (optionSimul == 3){ // Black-Scholes (sigma(t))
	price = BSpriceTimeDepVol(S_0,t_0,r,q,T,K,intSigma);
      }
      fprintf(fic_out_data,"%lf %lf %lf\n",K,T,price);
    }
    printf("--> prices saved in the file \"%s\"\n\n",name_out_data);
  }
}
