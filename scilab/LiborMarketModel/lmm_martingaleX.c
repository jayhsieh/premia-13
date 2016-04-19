#include"lmm_martingaleX.h"



//****************************************************************************************************
/*
/*   Added functions for interface
/*
/*****************************************************************************************************/

int lmm_caplet_terminalX_pricer(double *caplets_price ,  double* maturities , int nb_mat , int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor)
{
  
  Volatility *ptVol;
  Libor *ptLib;
  MartingaleX *ptX;

  int i;
  double  dt = tenor/nb_time_step;

  mallocLibor(&ptLib , nb_mat  , tenor );
  mallocVolatility(&ptVol , nb_factors );
  mallocMartingaleX( &ptX , ptLib );
 
  initLibor(ptLib);
  computeCapletTerminalX(nb_MC , dt , ptLib , strike , ptVol , caplets_price);
  printf("Caplet prices:\n");
  
  maturities[0]=0.0;
  
  
  for(i=1;i<nb_mat;i++)
    {
      //  printf("Caplet(Ti=%lf,%lf,K=%lf)=%lf \n",ptLib->maturity[i] , ptLib->maturity[i]+tenor , strike , caplets_price[i]);
      maturities[i]=ptLib->maturity[i];
    }
  
  
  freeLibor(&ptLib);
  freeVolatility(&ptVol);
  freeMartingaleX(&ptX);

  
  return(1);


}


int lmm_swaption_payer_terminalX_pricer(double *swaption_price ,  double swaption_maturity , double swap_maturity    , int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor)
{

  Volatility *ptVol;
  Libor *ptLib;
  MartingaleX *ptX, *ptX2;
  Swaption *ptSwpt;
  
  int numMat;
  double  dt=tenor/nb_time_step;
 
  numMat=(int)(swap_maturity/tenor);
  
  mallocLibor(&ptLib , numMat , tenor );
  mallocVolatility(&ptVol , nb_factors );
  mallocMartingaleX(&ptX, ptLib);
  mallocSwaption(&ptSwpt, swaption_maturity ,swap_maturity , 0.0 , strike , tenor );
 
  initLibor(ptLib);
  *swaption_price=computeSwaptionTerminalX(nb_MC , dt , ptLib , ptSwpt , ptVol );
  /*
    printf("Spt(T=%lf,%lf,K=%lf)=%lf \n", ptSwpt->swaptionMaturity , ptSwpt->swapMaturity , strike , *swaption_price );
  */

  return(1);

}



int mallocMartingaleX(MartingaleX **ptX , Libor *ptL)
{
  int i,N;
  double delta;
  MartingaleX *pt;

  N=ptL->numberOfMaturities;
  delta=ptL->tenor;
  
  pt=(MartingaleX *)malloc(sizeof(MartingaleX));
  pt->maturity=(double *)malloc((N+1)*sizeof(double));
  pt->Xvalue_0=(double *)malloc(N*sizeof(double));
  pt->Xvalue_t=(double *)malloc(N*sizeof(double));
   
  pt->NbMaturity=N;
  pt->tenor=delta;
  pt->Xvalue_0[N-1]=ptL->libor[N-1];
  pt->Xvalue_t[N-1]=ptL->libor[N-1];
  
  for(i=0; i<=N; i++){pt->maturity[i]=ptL->maturity[i];}
  for(i=2; i<=N; i++){pt->Xvalue_0[N-i]=pt->Xvalue_0[N-i+1]*ptL->libor[N-i]*(1+delta*ptL->libor[N-i+1])/ptL->libor[N-i+1];}
  for(i=2; i<=N; i++){pt->Xvalue_t[N-i]=pt->Xvalue_0[N-i];}
   
  *ptX=pt;
  return(1);
}
int initMartingaleX(MartingaleX *ptX)
{
  int i;
  
  for(i=0; i<ptX->NbMaturity; i++){ptX->Xvalue_t[i]=ptX->Xvalue_0[i];}
  return 1;
}

int freeMartingaleX(MartingaleX **ptX)
{
  MartingaleX  *pt;
  
  pt=(*ptX);
  free(pt->maturity);
  free(pt->Xvalue_0);
  free(pt->Xvalue_t);  
		
  *ptX=NULL;
  
  return(1);
}



int printMartingaleX(MartingaleX *ptX)
{
  int i,N;
  N=ptX->NbMaturity;
  
  printf("Maturity scale of the underlaying assets: \n");
  for(i=0;i<N;i++){printf("T=%lf : X= %lf -> %lf \n",ptX->maturity[i],ptX->Xvalue_0[i],ptX->Xvalue_t[i]);}
  printf("\n");
  
  return(1);
}

int drawMartingaleX(MartingaleX *ptX, Volatility *ptV, RandomGenerator *ptW, double dt, double t0, double t)
{
  int i,j,k,N, Nfac;
  double delta;
  double *sigma;
  double Di;
  double *ci;
  double T;
  double *Xprev;
  
  Nfac=ptV->numberOfFactors;
  N=ptX->NbMaturity;
  delta=ptX->tenor;
  k=0;
 
  sigma=(double*)malloc(Nfac*sizeof(double));
  Xprev=(double*)malloc(N*sizeof(double));
  ci=(double*)malloc(Nfac*sizeof(double));
  
  while(k*dt<t)
	{
	  Di=1;
	  randomVector(ptW);
	  T=ptX->maturity[N-1];

	  for(j=0;j<N;j++){Xprev[j]=ptX->Xvalue_t[j];}
	  for(j=0;j<Nfac;j++){ci[j]=0;}
	  for(j=0;j<Nfac;j++){sigma[j]=evalVolatility(ptV,j,t0+k*dt,T);}
	  
	  for(i=N-1; i>0; i--)
		{	
		  ptX->Xvalue_t[i]=Xprev[i]*exp(-0.5*ps(Nfac, sigma, sigma)*dt + sqrt(dt)*ps(Nfac, sigma, ptW->val));
		  Di+=Xprev[i]*delta;
		  if(t0+k*dt<T){ for(j=0; j<Nfac; j++){ci[j]=ci[j]+Xprev[i]*delta*evalVolatility(ptV,j,t0+k*dt,T)/Di;}}
		  T=ptX->maturity[i-1];
		  if(t0+k*dt<T){  for(j=0; j<Nfac; j++){sigma[j]=evalVolatility(ptV,j,t0+k*dt,T) + ci[j];}}
		  else { for(j=0; j<Nfac; j++){sigma[j]=ci[j];}}
		}
	  k+=1;
    }
  
  free(sigma);
  free(Xprev);
  free(ci);
  
  return(1);
}

int SetLiborWithMartingaleX(Libor *ptL, MartingaleX *ptX) 
{
  int i,j,N;
  double Di=1;
  
  N=ptX->NbMaturity;

  if(N!=ptL->numberOfMaturities)
	{"ERROR IN SetLiborWithMartingaleX() in file martingaleX.c, the number of maturities are not the same for L and X \n";}

  ptL->libor[N-1]=ptX->Xvalue_t[N-1];
 
  for(i=N-2; i>0; i--)
	{
	  Di+=ptX->Xvalue_t[i+1]*ptX->tenor;
	  ptL->libor[i]=ptX->Xvalue_t[i]/Di;
	}
  
  return(1);
  
}
  
double computeSwaptionTerminalX(int Nbmontecarlo,double dt,Libor *ptLib,Swaption *ptSwpt,Volatility *ptVol)
{
  int N,alpha, beta,i,j,k,l;
  double price;
  MartingaleX *ptX;
  RandomGenerator *ptRand;
  double tenor;
  double sumbond;
  double zc, Xab,SumXa,Bab,BN, SumXi;
  double t, et, p;
  
  tenor=ptLib->tenor;
  t=ptSwpt->swaptionMaturity;
  N=ptLib->numberOfMaturities;
  alpha=(int)(ptSwpt->swaptionMaturity/tenor);
  beta=(int)(ptSwpt->swapMaturity/tenor);
 
  mallocMartingaleX(&ptX, ptLib);
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors);
  et=0;
  price=0.0;
  // compute bond value P(0,SwaptionMat)
  zc=1.;
  initLibor(ptLib);
  
  for (i=0;i<N;i++)
    {
      zc*=1./(1.+tenor*ptLib->libor[i]);
    }
  initMartingaleX(ptX);
  
  for(l=0;l<Nbmontecarlo;l++)
    {
	  initMartingaleX(ptX);
	  drawMartingaleX(ptX,ptVol,ptRand,dt,0,t);	  

	  
	  Xab=0;
	  for(i=alpha;i<beta;i++)
		{
		  Xab+=ptX->Xvalue_t[i]; 
		}
	  SumXi=1;
	  for(j=N-1;j>beta;j--){SumXi+=tenor*ptX->Xvalue_t[j];}
	  sumbond=0;
	 
	  for(i=beta;i>alpha;i--)
		{
		  SumXi+=tenor*ptX->Xvalue_t[i];
		  sumbond+=SumXi;
		}
	  p=tenor*maximum(Xab - ptSwpt->strike*sumbond,0.0);
	  price+= p;
	  et+=p*p;
    }

  price=price/Nbmontecarlo;
  et=sqrt((et/Nbmontecarlo - price*price)/Nbmontecarlo);
  price=price*zc;
 

  //printSwaption(ptSwpt);
  //printf("ecart type: %lf\n",et*10000);
  
  //printf("Swaption price (simulated under the termimal measure thanks to a discrete martingale term): %lf\n",price*10000);
  freeRandomGenerator(&ptRand);
  freeMartingaleX(&ptX);

  return price;
}

int computeCapletTerminalX(int Nbmontecarlo,double dt,Libor *ptLib, double K, Volatility *ptVol, double *Caps)
{
  int i,l,j,k;
  int N;
  double B0, Bk;
  MartingaleX *ptX;
  RandomGenerator *ptRand;
  double tenor;
  
  tenor=ptLib->tenor;
 
  N=ptLib->numberOfMaturities;

  mallocMartingaleX(&ptX, ptLib);
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors); 
  
  // compute bond value P(0,SwaptionMat)
  B0=1.;
  Bk=1.;
 
  for (i=0;i<N;i++)
    {
      B0*=(1.+tenor*ptLib->libor[i]);
    }
  B0=1/B0;
  
  for(j=1; j<N; j++){Caps[j]=0;}
  
  for(l=0;l<Nbmontecarlo;l++)
    {
	  initMartingaleX(ptX);
	  drawMartingaleX(ptX,ptVol,ptRand,dt,0,tenor);
	  for(k=1;k<N-1;k++)
		{
		  drawMartingaleX(ptX,ptVol,ptRand,dt,k*tenor, tenor);	  
		  SetLiborWithMartingaleX(ptLib, ptX);
		  Bk=1;
		  for (i=k+1;i<N;i++)
			{
			  Bk*=(1.+tenor*ptLib->libor[i]);
			}
		  Caps[k]+=maximum(ptLib->libor[k]-K, 0)*Bk;

		}
	  drawMartingaleX(ptX,ptVol,ptRand,dt,(N-1)*tenor,tenor);
	  SetLiborWithMartingaleX(ptLib, ptX);

	  Caps[N-1]+=maximum(ptLib->libor[N-1]-K, 0);
	  
	}  
  for(j=1; j<N; j++){Caps[j]=Caps[j]*B0/Nbmontecarlo;}

  //printdoubles("Caplets.dat",N,Caps);
  freeRandomGenerator(&ptRand);
  freeMartingaleX(&ptX);
 
  return 1;

}

int computeBiasTerminalX(int Nbmontecarlo,double dt,Libor *ptLib1,double K, Volatility *ptVol, double *Bias)
{
  int i,j,k,l,m;
  int N=ptLib1->numberOfMaturities;
  int Nfac=ptVol->numberOfFactors;
  RandomGenerator *ptW;
  MartingaleX *ptX;
  Libor *ptLib2;
  double t,T;
  double **sigma;
  double delta=ptLib1->tenor;
  double Di,Bi,BN;
  double* cap1;
  double* cap2;
  double* ci;
  double* libor2;
  
  int it;

  ci=(double*)malloc(Nfac*sizeof(double));
  libor2=(double*)malloc(N*sizeof(double));
  cap1  =(double*)malloc(N*sizeof(double));
  cap2  =(double*)malloc(N*sizeof(double));
  sigma=(double**)malloc(N*sizeof(double*));
  
  for(i=0;i<N;i++){sigma[i]=(double*)malloc(Nfac*sizeof(double));}

  mallocRandomGenerator(&ptW,Nfac);
  mallocMartingaleX(&ptX, ptLib1);
  mallocLibor(&ptLib2,N, delta );
  copyLibor(ptLib1, ptLib2 );
  initLibor(ptLib2);

  printf("Nfactor=%d\n",Nfac);

  for (i=0;i<N;i++){cap1[i]=0;}
  for (i=0;i<N;i++){cap2[i]=0;}
  
  initLibor(ptLib1);
 
  for(l=0; l<Nbmontecarlo; l++)
	{
	  k=0;
	  initLibor(ptLib1);
	  initMartingaleX(ptX);
	  it=0;
	  m=1;
	  
	  while(k*dt<ptLib1->maturity[N-1])
		{
		  t=k*dt;
		  Di=1;
		  
		  for(j=0;j<Nfac;j++){ci[j]=0;}
		  
		  randomVector(ptW);

		  if(ptLib1->maturity[it]<=t){it=it+1;}
		  
		  for(i=N-1;i>=it;i--)
			{
			  T=ptLib1->maturity[i];
			  for(j=0; j<Nfac; j++){sigma[i][j]=evalVolatility(ptVol,j,t,T);}
			  ptLib1->libor[i]=ptLib1->libor[i]*exp(-0.5*ps(Nfac, sigma[i], sigma[i])*dt + sqrt(dt)*ps(Nfac, sigma[i], ptW->val));
			}

		  for(i=N-2; i>0; i--)
			{
			  Di+=ptX->Xvalue_t[i+1]*delta;

			  T=ptX->maturity[i+1];
			  if(t<T){ for(j=0; j<Nfac; j++){ci[j]=ci[j]+ptX->Xvalue_t[i+1]*delta*evalVolatility(ptVol,j,t,T)/Di;}}
			  T=ptX->maturity[i];
			  if(t<T){  for(j=0; j<Nfac; j++){sigma[i][j]=evalVolatility(ptVol,j,t,T) + ci[j];}}
			  else { for(j=0; j<Nfac; j++){sigma[i][j]=ci[j];}}
			
			  ptX->Xvalue_t[i]=ptX->Xvalue_t[i]*exp(-0.5*ps(Nfac, sigma[i], sigma[i])*dt + sqrt(dt)*ps(Nfac, sigma[i], ptW->val));
			}
		  T=ptX->maturity[N-1];
		  for(j=0; j<Nfac; j++){sigma[N-1][j]=evalVolatility(ptVol,j,t,T);}
		  ptX->Xvalue_t[N-1]=ptX->Xvalue_t[N-1]*exp(-0.5*ps(Nfac,sigma[N-1],sigma[N-1])*dt+sqrt(dt)*ps(Nfac,sigma[N-1],ptW->val));

		  if(t+dt==ptLib2->maturity[m])
			{
			  SetLiborWithMartingaleX(ptLib2,ptX);
			  libor2[m]=ptLib2->libor[m];
			  Bi=1;
			  for (i=m+1;i<N;i++)
				{
				  Bi*=(1.+delta*libor2[i]);
				}
			  //printf("BN(T_%d)=%lf\n",m+1,Bi);
			  cap2[m]+=maximum(libor2[m]-K, 0)*Bi;
			  m+=1;
			}
		  k+=1;
		}
	  
	  for (i=1;i<N;i++){cap1[i]+=maximum(ptLib1->libor[i]-K, 0);}
	}
  initLibor(ptLib1);

  BN=1;
  for(j=0;j<N;j++){BN=BN/(1+delta*ptLib1->libor[j]);}
  for(i=1;i<N;i++)
	{
	  Bi=1;
	  for(j=0;j<i+1;j++){Bi=Bi/(1+delta*ptLib1->libor[j]);}
	  Bias[i]=1-cap2[i]*BN/(cap1[i]*Bi);
	  if(Bias[i]<0){Bias[i]=-Bias[i];}
	 
	  
	} 

  Bias[0]=0;
  printdoubles("capletbias.dat",N,Bias);
  
  freeRandomGenerator(&ptW);
  freeMartingaleX(&ptX);
  freeLibor(&ptLib2);
  free(cap1);
  free(cap2);
  free(ci);
  free(libor2);
  for(i=0;i<N;i++){free(sigma[i]);}
  free(sigma);
  
  return 1;
}


int computeZCBondTerminalX(int Nbmontecarlo,double dt,Libor *ptLib, double maturity, Volatility *ptVol)
{
  int i,j,l;
  int N=ptLib->numberOfMaturities;
  double BN;
  double delta=ptLib->tenor;
  double *Di;
  double *CumulDi;
  double *B0;
  MartingaleX *ptX;
  RandomGenerator *ptW;

  initLibor(ptLib);
  mallocMartingaleX(&ptX, ptLib);
  mallocRandomGenerator(&ptW,ptVol->numberOfFactors);
  Di=(double*)malloc(N*sizeof(double));
  CumulDi=(double*)malloc(N*sizeof(double));
  B0=(double*)malloc((N+1)*sizeof(double));
 
  for(i=0;i<N;i++){CumulDi[i]=0;}
  for(i=0;i<=N;i++){B0[i]=1;}
  for(i=1;i<=N;i++){for(j=0;j<i;j++){B0[i]=B0[i]/(1+delta*ptLib->libor[j]);}}
  B0[0]=1;

  for(l=0;l<Nbmontecarlo;l++)
    {
	  initMartingaleX(ptX);
	 
	  for(i=1;i<N;i++)
		{
		  Di[i]=0;
		  drawMartingaleX(ptX,ptVol,ptW,dt,(i-1)*delta, delta);	  
		  for(j=i; j<N; j++){Di[i]+=ptX->Xvalue_t[j];}
		  Di[i]=1+delta*Di[i];
		}
	  for(i=0;i<N;i++){CumulDi[i]+=Di[i];}
	  
	}
 
   
  for(i=1;i<N;i++){CumulDi[i]=CumulDi[i]/Nbmontecarlo*B0[N];}
  CumulDi[0]=1;
  for(i=0;i<N;i++){B0[i]=(B0[i]-CumulDi[i])/B0[i];}
  B0[N]=0;
  printdoubles("zcbonderror.dat",N+1,B0);

  free(Di);
  free(B0);
  free(CumulDi);
  freeMartingaleX(&ptX);
  freeRandomGenerator(&ptW);
  
return 1;
}

int printdoubles(char *filename, int n, double *D)
{
  int  i,etat;
  FILE* file;

  file=fopen(filename, "w");

  for(i=0;i<n;i++){fprintf(file,"%lf ",D[i]);}
  fprintf(file, "\n");
  etat=fclose(file);

  return 1;
}
