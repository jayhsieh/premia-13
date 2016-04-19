
#include"lmm_header.h"



static double period;

double f1(double t ,double T)
{
  if (t>=T)
    return(0.0);
  else
    return(0.20);
}

double f2(double t ,double T)
{
  
  double val;
  if (t>=T)
    return(0.0);
  else
    val=1./sqrt(0.04+0.00075*t) * (0.01 - 0.05*exp(-0.1*(T-t)));
  return(val);
}


double fd2(double t, double T)
{
  int j;
  int k;
  double val;
  
  if (t>=T)
    return(0.0);
  else
    {
      
      j=(int)(t/period);
      k=(int)(T/period);
      val=0.01-0.05*exp(-0.1*(j-k))/sqrt(0.04+0.00075*j);
      return(val);
    }
}


double fd1(double t ,double T)
{
  if (t>=T)
    return(0.0);
  else
    return(0.20);
}



int mallocVolatility(Volatility **ptVol , int numOfFac )
{
  int i;
  Volatility *pt;
  pt=(Volatility *)malloc(sizeof(Volatility));

  if (numOfFac>2)
    {
      printf("only two factors allowed !!!");
      numOfFac=2;
    }
  pt->numberOfFactors=numOfFac;
  
  pt->vol=(funcVol *)malloc(sizeof(funcVol)*pt->numberOfFactors);
  if (numOfFac==1)
    {
      pt->vol[0]=f1;
    }
  else
    {
      pt->vol[0]=f1;
      pt->vol[1]=f2;
    }
      
  *ptVol=pt;
  return(1);
}


int mallocVolatilityInteger(Volatility **ptVol , int numOfFac, float tenor)
{
  int i;
  Volatility *pt;
  pt=(Volatility *)malloc(sizeof(Volatility));

  period = tenor;

  if (numOfFac>2)
    {
      printf("only two factors allowed !!!");
      numOfFac=2;
    }
  pt->numberOfFactors=numOfFac;
  
  pt->vol=(funcVol *)malloc(sizeof(funcVol)*pt->numberOfFactors);
  if (numOfFac==1)
    {
      pt->vol[0]=fd1;
    }
  else
    {
      pt->vol[0]=fd1;
      pt->vol[1]=fd2;
    }
  *ptVol=pt;
  return(1);
}

double evalVolatility( Volatility* ptVol ,int factorNumber, double t, double T)
// returns the value at time t of the volatility factor number factorNumber with maturity T 
{
  if(factorNumber>ptVol->numberOfFactors)
    {
      printf("not enough factors!\n");
      exit(1);
    }
  else
    {
      return(ptVol->vol[factorNumber](t,T));
    }
}


int freeVolatility(Volatility **ptVol)
{
  Volatility *pt;
  pt=(*ptVol);
  *ptVol=NULL;
  free(pt->vol);
   
  return(1);
}


int copyVolatility(Volatility *ptSrc, Volatility *ptDest)
{
  int i,j;
  
  ptDest->numberOfFactors=ptSrc->numberOfFactors; /*numOfFac;*/

  for (j=0;j<ptSrc->numberOfFactors;j++)
    {
      ptDest->vol[j]=ptSrc->vol[j];
    }
	  
  return(1);
};






