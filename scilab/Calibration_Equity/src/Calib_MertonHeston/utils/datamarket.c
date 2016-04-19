#include "datamarket.h"


DataMarket * CreateDataMarket(int N)
{
  int i;
  
  DataMarket * DM;
  DM          = (DataMarket*)malloc(sizeof(DataMarket));
  DM->nbdata  = N;
  DM->TypeOpt = (int*)malloc(sizeof(int)*N);
  DM->St0     = (double*)malloc(sizeof(double)*N);
  DM->K       = (double*)malloc(sizeof(double)*N);
  DM->T       = (double*)malloc(sizeof(double)*N);
  DM->r       = (double*)malloc(sizeof(double)*N);
  DM->d       = (double*)malloc(sizeof(double)*N);
  DM->wi      = (double*)malloc(sizeof(double)*N);
  DM->PrixObs = (double*)malloc(sizeof(double)*N);
  DM->SigmaImpObs = (double*)malloc(sizeof(double)*N);
  DM->SigmaImp = (double*)malloc(sizeof(double)*N);
  DM->PrixMod = (double*)malloc(sizeof(double)*N);

  // on initialise les poids a 1.
  for(i=0;i<N;i++) DM->wi[i] = 1.;
  
	
  return DM;
  
}


void AfficheDataMarket(DataMarket *DM)
{
  int i;
  
  printf("DataMarket : %d  donnees\n",DM->nbdata);
  printf("    TypeOpt,      St0,             T,             K,             PrixObs,          SigmaImpObs(%%)   \n");
  for(i=0;i<DM->nbdata;i++)
	{
	  printf("        %d     %f       %f       %f         %f         %f\n",DM->TypeOpt[i],DM->St0[i],DM->T[i],DM->K[i],DM->PrixObs[i],100.*DM->SigmaImpObs[i]);
	  
	  
	}
}

void FreeDataMarket(DataMarket *DM)
{
  //
  free(DM->TypeOpt);
  free(DM->St0);
  free(DM->K);
  free(DM->T);
  free(DM->r);
  free(DM->d);
  free(DM->wi);
  free(DM->PrixObs);
  free(DM->SigmaImpObs);
  free(DM->SigmaImp);
  free(DM->PrixMod);
  free(DM);
  //
  DM = NULL;
  //
}



