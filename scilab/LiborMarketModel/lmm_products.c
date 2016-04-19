
#include"lmm_header.h"

int mallocSwaption(Swaption **ptSwpt ,double swaptionMat , double swapMat ,double priceVal ,double K , double tenor )
{
  Swaption *pt;
  pt=(Swaption *)malloc(sizeof(Swaption));
  
  pt->swaptionMaturity=swaptionMat;
  pt->swapMaturity=swapMat;
  pt->strike=K;
  pt->price=priceVal;
  pt->tenor=tenor;
  pt->numberOfDates=(int)(pt->swapMaturity/tenor);

  *ptSwpt=pt;
  return(1);
}

int freeSwaption(Swaption **ptSwpt)
{
  Swaption  *pt;

  pt=(*ptSwpt);
  *ptSwpt=NULL;
  
  return(1);
}



int printSwaption(Swaption *pt)
{

  printf("swaption maturity %f\n", pt->swaptionMaturity);
  printf("swap maturity  %f\n", pt->swapMaturity);
  printf("price  %f\n", pt->price);
  printf("swaption strike %f\n", pt->strike);
  printf("\n");
}
  
  
////////////////////////
int mallocCaplet(Caplet **ptCplt ,double Mat , double priceVal , double K)
{
  Caplet *pt;
  pt=(Caplet *)malloc(sizeof(Caplet));
  
  pt->maturity=Mat;
  pt->strike=K;
  pt->price=priceVal;
  
  *ptCplt=pt;
  return(1);
}

int freeCaplet(Caplet **ptCplt)
{
  Caplet  *pt;

  pt=(*ptCplt);
  *ptCplt=NULL;
  
  return(1);
}



int printCaplet(Caplet *pt)
{

  printf("option maturity %f\n", pt->maturity);
  printf("price  %f\n", pt->price);
  printf("strike %f\n", pt->strike);
  printf("\n");
}




