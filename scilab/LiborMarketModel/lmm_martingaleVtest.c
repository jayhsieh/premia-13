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
#include"lmm_martingaleV.h"
#include"lmm_libor.h"
#include"lmm_random_generator.h"
#include"lmm_products.h"
#include"lmm_numerical.h"
#include"lmm_zero_bond.h"


int main()
{
  Volatility *ptVol;
  Libor *ptLib;
  MartingaleV *ptV, *ptV2;
  Swaption *ptSwpt;
  ZeroBond * ptZb;


  int N,i,j;
  double* E;

  double  tenor=0.5;
  double t=5;
  int numFac=1;
  int Nmc   =10000;
  int numMat;
  int numberTimeStep=10;
  double dt=tenor/numberTimeStep;
  double swaptionMat=5.;
  double swapMat=7.;
  double swaptionMat_start=3.;
  double priceVal=0.20;
  double K=0.050;
  int nb_test;
  double T,p;
  double swapRate;
  double price=0.0;
  double *Bias;
  double *Caps;

  
  char *sorties="outputsSpot.dat";
  int  etat;
  FILE* fich;

  //
  double* maturities;
  double swaption_price;
  double swaption_maturity=3;
  double swap_maturity=7.;
  int nb_MC   =10000;
  int nb_time_step=10;
  double   strike=0.05;
  int nb_factors=1;

  maturities=(double*)malloc(N*sizeof(double));

  
  fich=fopen(sorties, "w");
  
  numMat=(int)(swapMat/tenor);


  N=numMat;
  Bias=(double*)malloc(N*sizeof(double));
  Caps=(double*)malloc((N+1)*sizeof(double));
 
  //  lmm_caplet_spotV_pricer(Caps ,  maturities , numMat , Nmc , numFac , numberTimeStep , K ,  tenor);
  lmm_swaption_payer_spotV_pricer(&swaption_price ,  swaption_maturity , swap_maturity   ,  nb_MC ,  nb_factors ,  nb_time_step ,  strike , tenor);

  T=3;
  
  mallocLibor( &ptLib , numMat , tenor );
  mallocVolatility( &ptVol , numFac );
  mallocMartingaleV( &ptV , ptLib );
  mallocSwaption( &ptSwpt , swaptionMat , swapMat , 0.20, K , tenor );

  
  initMartingaleV(ptV);  
  initLibor(ptLib);
  
  printf("DEBUT - MARTINGALE SPOT METHOD.\n"); 
  printf("nb of Monte carlo drawings : %d \n", Nmc);
  etat=fprintf(fich,"nb of Monte carlo drawings : %d \n", Nmc);
  printf("PLEASE WAIT...\n");

  computeCapletSpotV(Nmc,dt,ptLib,K,ptVol, Caps);
  printf("Caplet prices:\n");
  etat=fprintf(fich,"Caplets prices:\n");
  for(i=1;i<N;i++){etat=fprintf(fich,"Caplet(T=%lf,%lf,K=%lf)=%lf \n",ptLib->maturity[i],ptLib->maturity[i]+tenor,K, Caps[i]);}
  for(i=1;i<N;i++){printf("Caplet(Ti=%lf,%lf,K=%lf)=%lf \n",ptLib->maturity[i],ptLib->maturity[i]+tenor,K, Caps[i]);}
  printf("\n");
  etat=fprintf(fich,"\n");	
  printf("Swaptions prices:\n");
  etat=fprintf(fich,"Swaptions prices:\n");

  nb_test=(int)((swapMat-swaptionMat_start)/tenor)-1;
  for(i=0;i<=nb_test;i++)
  {
    ptSwpt->swaptionMaturity=swaptionMat_start + tenor*i;
    initLibor(ptLib);
    p=computeSwaptionSpotV(Nmc,dt,ptLib,ptSwpt,ptVol);
    etat=fprintf(fich,"Spt(T=%lf,%lf,K=%lf)=%lf \n", ptSwpt->swaptionMaturity, ptSwpt->swapMaturity,K,p);
    printf("Spt(T=%lf,%lf,K=%lf)=%lf \n", ptSwpt->swaptionMaturity, ptSwpt->swapMaturity,K,p);
  }
  etat=fprintf(fich,"\n");
  etat=fclose(fich);	
  
  // initLibor(ptLib);
  // computeBiasTerminalX(Nmc,dt,ptLib,K,ptVol, Bias);
  
  //initLibor(ptLib);
  //computeZCBondSpotV(Nmc,dt,ptLib, T, ptVol);
   
  printf("Ces sorties sont enregistrees dans ouputsSpot.dat \n");
  printf("END\n");
    
  return(1);
}
