
#include"lmm_header.h"


int mallocZeroBond( ZeroBond **ptZb ,int numOfMat, double tenorVal )
{
  int i;
  ZeroBond *pt;
  pt=(ZeroBond *)malloc(sizeof(ZeroBond));
  
  pt->numberOfMaturities=numOfMat;
  pt->tenor=tenorVal;
  pt->value=(double *)malloc(sizeof(double)*pt->numberOfMaturities);
  for (i=0;i<pt->numberOfMaturities;i++)
    {
      pt->value[i]=0.0;
    }

  *ptZb=pt;
  return(1);
}

int freeZeroBond(ZeroBond **ptZb)
{
  ZeroBond  *pt;

  pt=(*ptZb);
  *ptZb=NULL;
  free(pt->value);
  
  return(1);
}


int printZeroBond(ZeroBond *ptZb)
{
  int i;
      
  for(i=0;i<ptZb->numberOfMaturities;i++)
    {
      printf("maturity %lf value %lf \n",(i+1)*ptZb->tenor , ptZb->value[i]);
    }

  printf("\n");

  return(1);
}

