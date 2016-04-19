#include <stdio.h>
#include <math.h>
#include "DupirePDE.h"
#include "inout.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.

The program IMPSIG computes the implied volatilities from given S_0, r, q, 
optionType, t_0 and option prices.

The data prices are contained in the file name_in_data which has three 
columns: K_i T_j V_{i,j} with V_{i,j} = V(K_i,T_j)

The implied volatilities are saved in the file name_out_sigma which has three 
columns: K_i T_j Sigma_{i,j} with Sigma_{i,j} = Implied_Sigma(K_i,T_j)

The input file IMPSIG.IN contains the following data:
#################### VARIABLES S_0, R, Q, OPTIONTYPE AND T_0 ###################
S_0 : price of the underlying asset at t=t_0
r : risk-free rate     
q : dividends (continuously compounded) 
optionType : type of the option (1 for call, 0 for put)  
t_0 : time origin
################################# FILES IN/OUT #################################
name_in_data : name of the file containing Ki, Tj and V_ij
name_out_sigma : name of the file containing Ki, Tj and impsigma_ij
*/

main(){

  /* declarations */
  double S_0,r,q,t_0;
  int optionType;
  char *name_in_impsig;
  char *name_in_data;
  char *name_out_sigma;
  double tol;

  int i,j;
  double K,T,V,sigma;

  FILE *fic_in_impsig;
  FILE *fic_out_sigma;
  FILE *fic_in_data;

  tol = 1.0e-6;
 
  name_in_impsig = "impsig.in";
 
  /* initialization of S_0, r, q,..., name_out_sigma from the file name_in_impsig */
  loadParamImpsig(name_in_impsig,&S_0,&r,&q,&optionType,&t_0,&name_in_data,&name_out_sigma);

  fic_in_data = fopen(name_in_data,"r");
  fic_out_sigma = fopen(name_out_sigma,"w");

 while (fscanf(fic_in_data,"%lf %lf %lf",&K,&T,&V) == 3){
   sigma = computeImpsig(S_0,r,q,optionType,t_0,K,T,V,tol);
   //printf("K = %lf, T = %lf, V = %lf, sigma = %lf\n",K,T,V,sigma);
   fprintf(fic_out_sigma,"%lf %lf %lf\n",K,T,sigma);
 }

 printf("\n--> implied volatilities saved in the file \"%s\"\n\n",name_out_sigma);
  
}
