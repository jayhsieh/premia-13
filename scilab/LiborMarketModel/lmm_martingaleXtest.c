/*************************************************************************************************************************************
 *
 *
 *
 *   Simulation in the libor market model of a  swaption 
 *          method:  Glasserman and Zhao Discrete Martingale - Monte Carlo 
 *
 *   Marc Barton-Smith
 *   June  2004
 *
 *
 *
 *
 ************************************************************************************************************************************/
#include<stdio.h>


#include"lmm_header.h"
#include"lmm_volatility.h"
#include"lmm_martingaleX.h"
#include"lmm_libor.h"
#include"lmm_random_generator.h"
#include"lmm_products.h"
#include"lmm_numerical.h"
#include"lmm_zero_bond.h"
#include"cumulfunc.h"

int main()
{
  Volatility *ptVol;
  Libor *ptLib;
  MartingaleX *ptX, *ptX2;
  Swaption *ptSwpt;
  ZeroBond * ptZb;

  int N,i,j;
  double* E;

  double  tenor=0.5;
  int numFac=2;
  int Nmc   =10000;
  int numMat;
  int numberTimeStep=10;
  double  dt=tenor/numberTimeStep;
  double  swaptionMat_start=3.;
  double  swapMat=7.;
  double  priceVal=0.20;
  double K=0.050;
  double T,p;
  double swapRate;
  double price=0.0;
  double *Bias;
  double *Caps;
  int nb_test;
  char *sorties="outputsTerm.dat";
  int  etat;
  FILE* fich;
  double q,x,bound;
  int status;
  double* maturities;
  double swaption_price;
  double swaption_maturity=3;
  double swap_maturity=7.;
  int nb_MC   =10000;
  int nb_time_step=10;
  double  strike=0.05;
  int nb_factors=1;
  
  
  fich=fopen(sorties, "w");

  numMat=(int)(swapMat/tenor);
  N=numMat;
  Bias=(double*)malloc(N*sizeof(double));
  
  Caps=(double*)malloc(N*sizeof(double));
  maturities=(double*)malloc(N*sizeof(double));

  //lmm_caplet_terminalX_pricer(Caps ,  maturities , numMat , Nmc , numFac , numberTimeStep , K ,  tenor);
  //lmm_swaption_payer_terminalX_pricer(&swaption_price ,  swaption_maturity , swap_maturity   ,  nb_MC ,  nb_factors ,  nb_time_step ,  strike , tenor);
  
  T=3;
  mallocLibor(&ptLib , numMat , tenor );
  mallocVolatility(&ptVol , numFac );
  mallocMartingaleX(&ptX, ptLib);
  mallocSwaption(&ptSwpt, swaptionMat_start,swapMat, 0.20,K,tenor);
 
  printf("STARTS - MARTINGALE TERMINAL METHOD.\n");
  printf("nb of Monte carlo drawings : %d \n", Nmc);
  etat=fprintf(fich,"nb of Monte carlo drawings : %d \n", Nmc);
  printf("PLEASE WAIT...\n");

	initLibor(ptLib);
	computeCapletTerminalX(Nmc,dt,ptLib,K,ptVol, Caps);
	printf("Caplet prices:\n");
	etat=fprintf(fich,"Caplet prices:\n");
	for(i=1;i<N;i++){etat=fprintf(fich,"Caplet(T=%lf,%lf,K=%lf)=%lf \n",ptLib->maturity[i],ptLib->maturity[i]+tenor,K, Caps[i]);}
	for(i=1;i<N;i++){printf("Caplet(Ti=%lf,%lf,K=%lf)=%lf \n",ptLib->maturity[i],ptLib->maturity[i]+tenor,K, Caps[i]);}
	printf("\n");
	etat=fprintf(fich,"\n");
	printf("Swaptions prices:\n");
	etat=fprintf(fich,"Swaptions prices:\n");
	nb_test=(int)((swapMat-swaptionMat_start)/tenor)-1;
	for(i=0;i<nb_test;i++)
	  {
		ptSwpt->swaptionMaturity = swaptionMat_start + tenor*i;
		initLibor(ptLib);
		p=computeSwaptionTerminalX(Nmc,dt,ptLib,ptSwpt,ptVol);
		etat=fprintf(fich,"Spt(T=%lf,%lf,K=%lf)=%lf \n", ptSwpt->swaptionMaturity, ptSwpt->swapMaturity,K,p);
		printf("Spt(T=%lf,%lf,K=%lf)=%lf \n", ptSwpt->swaptionMaturity, ptSwpt->swapMaturity,K,p);
	  }
	etat=fprintf(fich,"\n");
	etat=fclose(fich);
	
	//initLibor(ptLib);
	//computeBiasTerminalX(Nmc,dt,ptLib,K,ptVol, Bias);

  
	//initLibor(ptLib);
	// computeZCBondTerminalX(Nmc,dt,ptLib, T, ptVol);
   
	printf("Output prices were put in ouputsTerm.dat \n");
	printf("END\n");
	//computeZCBondTerminalX(Nmc,dt,ptLib,T,ptVol);
	//printLibor(ptLib);
  
  return(1);
}

