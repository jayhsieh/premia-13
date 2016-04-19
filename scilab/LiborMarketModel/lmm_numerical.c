#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include"lmm_numerical.h"

///
double BSFormula(Swaption *ptSwpt, Libor* ptLib, double evalTime, double blackVol)
{ // pricing a swaption in the black model

  double Sigma;
  double borne=7.;
  double d1;
  double d2;
  double sum;
  int o,s,m;
  double underlying;
  
  
  s=(int)(ptSwpt->swaptionMaturity/ptLib->tenor);
  m=(int)(ptSwpt->swapMaturity/ptLib->tenor);
  o=(int)(evalTime/ptLib->tenor);

  sum=computeZeroCouponSum(ptLib, o , s+1 , m );
  underlying=computeSwapRate(ptLib, o ,  s, m );

  Sigma=blackVol*sqrt(ptSwpt->swaptionMaturity);
  d1= (log(underlying/ptSwpt->strike)+ 0.5*pow(Sigma,2))/Sigma;
  d2=( log(underlying/ptSwpt->strike)- 0.5*pow(Sigma,2))/Sigma ;
  
  if ( (d1<borne) && (d1>-borne))
    {
      return(sum*(underlying*N(d1) -ptSwpt->strike*N(d2) ));
      
    }
  else
    {
      printf(" can not compute swaption price\n");
      return(0.);
    }
  
}

double ps(int n, double *u, double *v)
{
  int l;
  double s=0;
  
  for(l=0;l<n;l++){s+=u[l]*v[l];}

  return s;
}


double maximum(double a,double b)
{
  if (a>b)
    return(a);
  else
    return(b);
}

//


double maxi(double a,double b)
{
  if (a>b)
    return(a);
  else
    return(b);
}

double ppos(double x){
  return(maxi(x,0));
}

int Set_to_Zero(double *ptr,int dim){
  int l;
  for(l=0;l<dim;l++) ptr[l]=0.0;
  return(1);
}

/*******************   Evolution Routines      ******************************/
int evolutionUnderSpotMeasure(RandomGenerator *ptRand, Libor* ptLibOld , Libor* ptLibNew, Volatility *ptVol, double dt , double t)
{
  // computes the evolution of libor rates under the spot measure
  int i,k,l;
  double val=0.0;
  double drift=0.0;
  double vol=0.0;
  double normVolatility=0.0;
  double scalarProductVol=0.0;
  double T_i,T_k;
  double v_i;
  

  for(i=1;i<ptLibNew->numberOfMaturities;i++)
    {
      if (ptLibOld->libor[i]==0.0)
	{
	  ptLibNew->libor[i]=0.0;
	}
      else
	{
	  // compute the drift
	  drift=0.0;
	  for(k=1;k<=i;k++)
	    {
	      scalarProductVol=0.0;
	      for(l=0;l<ptVol->numberOfFactors;l++)
		{
		  //scalarProductVol+= ptVol->vol[i-1][l]*ptVol->vol[k-1][l];
		  T_i=ptLibOld->maturity[i];
		  T_k=ptLibOld->maturity[k];
		  scalarProductVol+= evalVolatility(ptVol,l,t,T_i)*evalVolatility(ptVol,l,t,T_k);
		  
		}
	      
	      drift+= scalarProductVol * ptLibOld->libor[k]* ptLibOld->tenor ;
	      drift/=(1.+ ptLibOld->tenor * ptLibOld->libor[k]) ; 
	    }

	  // compute de square of the volatility and the random choc
	  normVolatility=0.0;
	  vol=0.0;
	  for(l=0;l<ptVol->numberOfFactors;l++)
	    {
	      T_i=ptLibOld->maturity[i];
	      v_i=evalVolatility(ptVol,l,t,T_i);
	      normVolatility += pow(v_i ,2);
	      vol+= v_i* ptRand->val[l];
	    }
	  
	  drift+=  (-0.5*normVolatility);
	  drift*=dt;
	  val=(drift + sqrt(dt)*vol);

	  //update 
	  ptLibNew->libor[i]=ptLibOld->libor[i]*exp(val);
	}
    }
  return(1);
}

int evolutionUnderForwardMeasure(RandomGenerator *ptRand, Libor* ptLibOld , Libor* ptLibNew, Volatility *ptVol, double dt , double t){
// computes the evolution of libor rates under the forward measure 
// corresponding to the swap Maturity: numeraire is B(t,T_e)
// Forward rate L(t;T_{e-1},T_e)=ptLib->Libor[e-1] is a martingale
  return(1);
}

/***************            Numeraire&Evolution routine     ******************************/

//Computation of Numeraire on the k-th MC path
//NumeraireSpot(T_j)=RollOverBond=Prod_{i=0}^{j-1}[1+tenor*L(T_i;i,i+1)]
//Numeraire(T_j)=ZeroCoupBond(T_j,T_numberofmaturities)

void computeNumeraire(char*MeasureName,Libor* ptLib,Swaption* ptSwpt,double *Numeraire,int j,int k,double auxspot){
  int l,s;
  double aux=1.0;
  s=(int)(ptSwpt->swaptionMaturity/ptLib->tenor);
  
  if(strcmp(MeasureName,"Spot")==0){
    if(((s-2)<=j)&&(j<=(ptLib->numberOfMaturities-3))) Numeraire[k*(ptLib->numberOfMaturities-s)+(j+2-s)]=auxspot;
		  }
  else{
    if ((s-1)<=j){
      for(l=j+1;l<ptLib->numberOfMaturities;l++) aux*=(1./(1.+ptLib->tenor*ptLib->libor[l]));
      Numeraire[k*(ptLib->numberOfMaturities-s)+(j+1-s)]=aux;
    }
  }
  return;
}

void Name_To_Measure(char *ErrorMessage, char *name,
		     int (**computeEvol)(RandomGenerator* ptRand,Libor* ptLibOld,Libor* ptLibNew,Volatility* ptVol,double dt,double t))
{
  /*initialization of evolution */
  if (strcmp("Spot",name)==0){
    *computeEvol=evolutionUnderSpotMeasure;
  } else if (strcmp("Fwd",name)==0){
    *computeEvol=evolutionUnderForwardMeasure;
  }  else {
    strcat(ErrorMessage,"Measure Error:");
    strcat(ErrorMessage,"(");	
    strcat(ErrorMessage,name);	
    strcat(ErrorMessage,") is not good. Please, try with a valid Measure Name");
  }
  return;
}





/****  Pricing routines  ****************************/

int computeSwaptionPriceSpot( int numberMonteCarloSim, int numberTimeStep , double dt ,Libor *ptLib , Swaption *ptSwpt, Volatility *ptVol )
{
  int i,l,j;
  double p,price;
  Libor* ptLibOld;
  Libor* ptLibNew;
  Volatility *ptVolOld;
  RandomGenerator *ptRand;
  double tenor;
  double swapRate;
  double t=0.0;
  /******/double err=0.0;
  tenor=ptLib->tenor;
    
   
  mallocLibor(&ptLibNew , ptLib->numberOfMaturities , ptLib->tenor );
  copyLibor(ptLib, ptLibNew );
  
  mallocLibor(&ptLibOld , ptLib->numberOfMaturities , ptLib->tenor );
  copyLibor(ptLib, ptLibOld);
  
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors);
  
  price=0.0;
  for(l=0;l<numberMonteCarloSim;l++)
    {
      t=0.0;
      copyLibor(ptLibOld,ptLib);
      
      p=1./(1.+ ptLib->tenor*ptLib->libor[0]);
      
      for(i=0;i<(int)(ptSwpt->swaptionMaturity/tenor);i++)
	{
	  for(j=0;j<numberTimeStep;j++) 
	    {
	      
	      randomVector(ptRand);
	      evolutionUnderSpotMeasure(ptRand,ptLib,ptLibNew,ptVol,dt,t);
	      copyLibor(ptLibNew,ptLib);
	      t+=dt;
	    }
	  //p*= 1./(1.+ ptLib->tenor*ptLib->libor[i+1]);
	  /****/if (i<((int)(ptSwpt->swaptionMaturity/tenor)-1)) p*= 1./(1.+ ptLib->tenor*ptLib->libor[i+1]);
	  if(i==((int)(ptSwpt->swaptionMaturity/tenor)-1))
	    break;
	  //  putLiborToZero(ptLib,i+1);
	  /****/      putLiborToZero(ptLib,i+1);
	}
	  
      
      swapRate=computeSwapRate(ptLib,(int)(ptSwpt->swaptionMaturity/tenor), (int)(ptSwpt->swaptionMaturity/tenor),(int)(ptSwpt->swapMaturity/tenor));
      
      p*= maxi(swapRate-ptSwpt->strike,0.0)*computeZeroCouponSum(ptLib, (int)(ptSwpt->swaptionMaturity/tenor) , (int)(ptSwpt->swaptionMaturity/tenor) +1  ,(int)(ptSwpt->swapMaturity/tenor) );
      
	
      price+=p;
      /*****/err+=p*p;
    }
  price/=numberMonteCarloSim;
  /***/err/=(double)numberMonteCarloSim;
  /***/err-=(price*price);err/=(double)numberMonteCarloSim;err=sqrt(err);
  printSwaption(ptSwpt);
  printf(" the price of the swaption (simulated under the spot measure) is \n%lf +/- %lf\n",price*10000,err*10000*1.96);


}



int computeSwaptionPriceForward( int numberMonteCarloSim, int numberTimeStep , double dt ,Libor *ptLib , Swaption *ptSwpt, Volatility *ptVol )
{
  int i,l,j;
  double p,price;
  Libor* ptLibOld;
  Libor* ptLibNew;
  Volatility *ptVolOld;
  RandomGenerator *ptRand;
  double tenor;
  double swapRate;
  double zc;
  double t;
  
  tenor=ptLib->tenor;
    
  mallocVolatility(&ptVolOld , ptVol->numberOfFactors );
  copyVolatility(ptVolOld, ptVol);
   
  mallocLibor(&ptLibNew , ptLib->numberOfMaturities , ptLib->tenor );
  copyLibor(ptLibNew, ptLib);
  
  mallocLibor(&ptLibOld , ptLib->numberOfMaturities , ptLib->tenor );
  copyLibor(ptLibOld, ptLib);
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors);

  price=0.0;
  // compute bond value
  zc=1.;
 
  for (i=0;i<(int)(ptSwpt->swaptionMaturity/tenor);i++)
    {
      zc*=1./(1.+tenor*ptLib->libor[i]);
    }
  for(i=1;i<(int)(ptSwpt->swaptionMaturity/tenor);i++)
    {
      putLiborToZero(ptLibOld,i);
    }

  copyLibor(ptLib,ptLibOld);

  for(l=0;l<numberMonteCarloSim;l++)
    {
      
      t=0.0;
      copyLibor(ptLib,ptLibOld);
      copyVolatility(ptVol,ptVolOld);

      
      for(i=0;i<(int)(ptSwpt->swaptionMaturity/tenor);i++)
	{
	  for(j=0;j<numberTimeStep;j++) 
	    {
	      randomVector(ptRand);
	      evolutionUnderSpotMeasure(ptRand,ptLib,ptLibNew,ptVol,dt,t);
	      copyLibor(ptLib,ptLibNew);
	      t+=dt;
	    }
	  
	}
	
      swapRate=computeSwapRate(ptLib,(int)(ptSwpt->swaptionMaturity/tenor), (int)(ptSwpt->swaptionMaturity/tenor),(int)(ptSwpt->swapMaturity/tenor));
      
      price+= maxi(swapRate-ptSwpt->strike,0.0)*computeZeroCouponSum(ptLib, (int)(ptSwpt->swaptionMaturity/tenor) , (int)(ptSwpt->swaptionMaturity/tenor) +1  ,(int)(ptSwpt->swapMaturity/tenor) );
      
    }
  price/=numberMonteCarloSim;

  printSwaption(ptSwpt);
  printf(" the price of the swaption (simulated under the forward measure) is %lf \n",zc*price*10000);


}


int computeZeroBond( int numberMonteCarloSim, int numberTimeStep , double dt ,Libor *ptLib , Volatility *ptVol , ZeroBond * ptZb)
{
  int i,l,j,k;
  double p,price;
  Libor* ptLibOld;
  Libor* ptLibNew;
  Volatility *ptVolOld;
  RandomGenerator *ptRand;
  double tenor;
  double *tmpVal;
  double t;

  tenor=ptLib->tenor;
  
  tmpVal=(double *)malloc(sizeof(double)*ptZb->numberOfMaturities);
  for (i=0;i<ptZb->numberOfMaturities;i++)
    {
      tmpVal[i]=0.0;
    }

  mallocVolatility(&ptVolOld , ptVol->numberOfFactors  );
  copyVolatility(ptVolOld, ptVol);
   
  mallocLibor(&ptLibNew , ptLib->numberOfMaturities , ptLib->tenor );
  copyLibor(ptLibNew, ptLib);
  
  mallocLibor(&ptLibOld , ptLib->numberOfMaturities , ptLib->tenor );
  copyLibor(ptLibOld, ptLib);
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors);

  for(l=0;l<numberMonteCarloSim;l++)
    {
      t=0.0;
      copyLibor(ptLib,ptLibOld);
      copyVolatility(ptVol,ptVolOld);
      
      for(k=0;k<ptZb->numberOfMaturities;k++)
	{
	  tmpVal[k]= 1./(1.+ ptLib->tenor*ptLib->libor[0]);
	}
      
      for(i=0;i< ptZb->numberOfMaturities  ;i++)
	{
	  for(j=0;j<numberTimeStep;j++) 
	    {
	      randomVector(ptRand);
	      evolutionUnderSpotMeasure(ptRand,ptLib,ptLibNew,ptVol,dt,t);
	      copyLibor(ptLib,ptLibNew);
	      t+=dt;
	    }
	  
	  for(k=i+1 ; k<ptZb->numberOfMaturities ; k++)
	    {
	      tmpVal[k]*= 1./(1.+ ptLib->tenor*ptLib->libor[k]);
	    }
      
	  putLiborToZero(ptLib,i+1);
	  
	}
      for(k=0;k<ptZb->numberOfMaturities;k++)
	{
	  ptZb->value[k]+=tmpVal[k];
	}
    }
  
  for(k=0;k<ptZb->numberOfMaturities;k++)
    {
      ptZb->value[k]/=(numberMonteCarloSim-1);
    }
  printf(" the price of the swaption (simulated under the spot measure) is %lf \n",1.);


}

////////////////////////////////////////////////////////////////////////////////////////////////////
int computeTrajectories( int numberTimeStep , double dt ,HistLibor *ptHistLib , Volatility *ptVol )
{
  int i,l,j,u;
  Libor* ptLibOld;
  Libor* ptLibNew;
  Libor* ptLib;
  RandomGenerator *ptRand;
  double t=0.0;
  
   
  mallocLibor(&ptLibNew , ptHistLib->numberOfMaturities , ptHistLib->tenor );
  //copyLibor(ptLib, ptLibNew );
  mallocLibor(&ptLibOld , ptHistLib->numberOfMaturities , ptHistLib->tenor );
  //copyLibor(ptLib, ptLibOld);
  mallocLibor(&ptLib , ptHistLib->numberOfMaturities , ptHistLib->tenor );
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors);
  
  for(l=0;l<ptHistLib->numberOfTrajectories;l++)
    {
      t=0.0;
      copyLibor(ptLibOld,ptLib);
      for(u=0;u<ptLib->numberOfMaturities;u++)
	{
	  ptHistLib->libor[l][0][u]=ptLib->libor[u];
	}
	  
      for(i=1;i<ptHistLib->numberOfMaturities;i++)
	{
	  
	  for(j=0;j<numberTimeStep;j++) 
	    {
	      
	      randomVector(ptRand);
	      evolutionUnderSpotMeasure(ptRand,ptLib,ptLibNew,ptVol,dt,t);
	      copyLibor(ptLibNew,ptLib);
	      t+=dt;
	    }
	
	  for(u=0;u<ptLib->numberOfMaturities;u++)
	    {
	      ptHistLib->libor[l][i][u]=ptLib->libor[u];
	    }
	  
	}
	  
      
    }


  freeRandomGenerator(&ptRand);
  freeLibor(&ptLib);
  freeLibor(&ptLibNew);
  freeLibor(&ptLibOld);

}
