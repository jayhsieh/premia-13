#include"lmm_martingaleV.h"



//****************************************************************************************************
/*
/*   Added functions for interface
/*
/*****************************************************************************************************/

int lmm_caplet_spotV_pricer(double *caplets_price ,  double* maturities , int nb_mat , int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor)
{

  Volatility *ptVol;
  Libor *ptLib;
  MartingaleV *ptV;

  int i;
  double  dt = tenor/nb_time_step;

  
  mallocLibor(&ptLib , nb_mat  , tenor );
  mallocVolatility(&ptVol , nb_factors );
  mallocMartingaleV( &ptV , ptLib );

  initMartingaleV(ptV);  
  initLibor(ptLib);

  computeCapletSpotV(nb_MC , dt , ptLib , strike , ptVol , caplets_price);
  maturities[0]=0.0;

  
  //  printf("Caplet prices:\n");
    for(i=1;i<nb_mat;i++)
    {
      //    printf("Caplet(Ti=%lf,%lf,K=%lf)=%lf \n",ptLib->maturity[i] , ptLib->maturity[i]+tenor , strike , caplets_price[i]);
      maturities[i]=ptLib->maturity[i];
    }
  

  freeLibor(&ptLib);
  freeVolatility(&ptVol);
  freeMartingaleV(&ptV);

  return(1);
}


int lmm_swaption_payer_spotV_pricer(double *swaption_price ,  double swaption_maturity , double swap_maturity    , int nb_MC , int nb_factors , int nb_time_step , double strike , double tenor)
{

  Volatility *ptVol;
  Libor *ptLib;
  MartingaleV *ptV;
  Swaption *ptSwpt;
  
  int numMat;
  double  dt=tenor/nb_time_step;

  numMat=(int)(swap_maturity/tenor);
  
  mallocLibor(&ptLib , numMat , tenor );
  mallocVolatility(&ptVol , nb_factors );
  mallocSwaption(&ptSwpt, swaption_maturity ,swap_maturity , 0.0 , strike , tenor );
  mallocMartingaleV( &ptV , ptLib );

  initLibor(ptLib);
  
  initMartingaleV(ptV);  
   *swaption_price=computeSwaptionSpotV(nb_MC , dt , ptLib , ptSwpt , ptVol );
   /*
     printf("Spt(T=%lf,%lf,K=%lf)=%lf \n", ptSwpt->swaptionMaturity , ptSwpt->swapMaturity , strike , *swaption_price );
   */

  return(1);

}


///////////////////////////////////////////////////////////////////////////////////////////////


int mallocMartingaleV(MartingaleV **ptV , Libor *ptL)
{
  /*Allocate and set the values of the MartingaleV asset with the libor datas under pointer ptL*/
  int i,N;
  double Pi,delta;
  MartingaleV *pt;

  N=ptL->numberOfMaturities;
  delta=ptL->tenor;

  pt=(MartingaleV *)malloc(sizeof(MartingaleV));
  pt->maturity=(double *)malloc((N+1)*sizeof(double));
  pt->Vvalue_0=(double *)malloc((N+1)*sizeof(double));
  pt->Vvalue_t=(double *)malloc((N+1)*sizeof(double));
   
  pt->NbMaturity=N;
  pt->tenor=delta;
  pt->Vvalue_0[0]=ptL->libor[0]; 
  /*this value does not interfere in the formulas(it is the present spot rate)*/
   
  for(i=0; i<=N; i++){pt->maturity[i]=ptL->maturity[i];}

  Pi=1;
  for(i=1; i<N; i++)
	{
	  Pi=Pi*(1+delta*ptL->libor[i]);
	  pt->Vvalue_0[i]=delta*ptL->libor[i]/Pi;
	}
  pt->Vvalue_0[N]=1/Pi;
  /*MartingaleV has one more terminal value than the Libor or MartingaleX*/

  for(i=0; i<=N; i++){pt->Vvalue_t[i]=pt->Vvalue_0[i];}
  /*Value_t is supposed to change along the computations, whereas Value_0 is supposed to remain unchanged*/
  
  *ptV=pt;
  return(1);
}

int initMartingaleV(MartingaleV *ptV)
{
  /* Reinitialize the time t value to the time 0 value*/
  int i;
  
  for(i=0; i<=ptV->NbMaturity; i++){ptV->Vvalue_t[i]=ptV->Vvalue_0[i];}
  return 1;
}

int freeMartingaleV(MartingaleV **ptV)
{
  MartingaleV  *pt;
  
  pt=(*ptV);
  free(pt->maturity);
  free(pt->Vvalue_0);
  free(pt->Vvalue_t);  
		
  *ptV=NULL;
  
  return(1);
}



int printMartingaleV(MartingaleV *ptV)
{
  int i,N;
  N=ptV->NbMaturity;
  
  printf("Maturity scale of the underlaying assets: \n");
  for(i=1;i<=N;i++){printf("T=%lf : V= %lf -> %lf \n",ptV->maturity[i-1],ptV->Vvalue_0[i],ptV->Vvalue_t[i]);}
  printf("\n");
  
  return(1);
}

int SetMartingaleVWithLibor(MartingaleV *ptV , Libor *ptL)
{
  int i,N;
  double Pi,delta;

  N=ptL->numberOfMaturities;
  delta=ptL->tenor;
  if(N!=ptV->NbMaturity)
	{"ERROR IN SetMartingaleVwithLibor() in file martingaleV.c, the number of maturities are not the same for L and V \n";}
  if(N!=ptV->NbMaturity)
	{"ERROR IN SetMartingaleVwithLibor() in file martingaleV.c, the tenor of L and V  are not the same\n";}

  ptV->Vvalue_t[0]=ptL->libor[0];
  /*Value_t is supposed to change along the computations, whereas Value_0 is supposed to remain unchanged*/
  
  Pi=1;
  for(i=1; i<N; i++)
	{
	  Pi=Pi*(1+delta*ptL->libor[i]);
	  ptV->Vvalue_t[i]=delta*ptL->libor[i]/Pi;
	}
  ptV->Vvalue_t[N]=1/Pi;
  
  return(1);
}
int SetLiborWithMartingaleV(Libor *ptL, MartingaleV *ptV) 
{
  int i,j,N;
  double Di=0;
  double delta;
  
  N=ptV->NbMaturity;
  delta=ptV->tenor;
  if(N!=ptL->numberOfMaturities)
	{"ERROR IN SetLiborWithMartingaleV() in file martingaleV.c, the number of maturities are not the same for L and V \n";}
  if(N!=ptL->numberOfMaturities)
	{"ERROR IN SetLiborWithMartingaleV() in file martingaleV.c, the tenor ot L and V are not the same \n";}
  
  for(i=N-1; i>0; i--)
	{
	  Di+=ptV->Vvalue_t[i+1];
	  ptL->libor[i]=ptV->Vvalue_t[i]/(delta*Di);
	}
  /*Value_t is supposed to change along the computations, whereas Value_0 is supposed to remain unchanged*/

  return(1);
  
}

int drawMartingaleV(MartingaleV *ptV, Volatility *ptVol, RandomGenerator *ptW, double dt, double t0,double t)
{
  int i,j,k,N, eta, Nfac;
  double delta, SumVeta;
  double Di;
  double *ci;
  double *sigma;
  double T;
  double *Vprev;
  
  Nfac=ptVol->numberOfFactors;
  N=ptV->NbMaturity;
  eta=0;
  while(t0>ptV->maturity[eta]){eta+=1;}
  delta=ptV->tenor;
  k=0;
 
  Vprev=(double*)malloc((N+1)*sizeof(double));
  ci=(double*)malloc(Nfac*sizeof(double));
  sigma=(double*)malloc(Nfac*sizeof(double));
  for(i=1;i<N+1;i++){Vprev[i]=ptV->Vvalue_t[i];}
  
  while(t0+k*dt<t)
	{
	  /*Diffusion between time t0 and t of the martingaleV value Vvalue_t*/
	  
	  randomVector(ptW);
	  /*normal variable draw at each time step*/

	  if(t0+k*dt>=ptV->maturity[eta]){eta+=1;}

	  Vprev[0]=0;
	  SumVeta=-1;
	  
	  for(j=0;j<Nfac;j++){ci[j]=0;}
	  
	  for(i=1; i<N; i++)
		{
		  /* Compution of each N-1 MartingaleV assets at time step t0+(k+1)*dt */
		  SumVeta+=Vprev[i-1];
		  T=ptV->maturity[i];
		  for(j=0;j<Nfac;j++){ci[j]+=evalVolatility(ptVol,j,t0+k*dt,T)*Vprev[i]/SumVeta;}
		  for(j=0;j<Nfac;j++){sigma[j]=ci[j]+evalVolatility(ptVol,j,t0+k*dt,T);}
		  
		  ptV->Vvalue_t[i]=Vprev[i]*exp(-0.5*ps(Nfac,sigma,sigma)*dt + sqrt(dt)*ps(Nfac,sigma, ptW->val));  
		}
	  ptV->Vvalue_t[N]=Vprev[N]*exp(-0.5*ps(Nfac,ci,ci)*dt + sqrt(dt)*ps(Nfac,ci, ptW->val));
	  /* Compution of the last MartingaleV asset at time step t0+(k+1)*dt */
	   
	  for(i=1; i<=N; i++){ Vprev[i]=ptV->Vvalue_t[i];} 
	  k+=1;
    }

  free(sigma);
  free(Vprev);
  free(ci);
  
  return(1);
}


double computeSwaptionSpotV(int Nbmontecarlo,double dt,Libor *ptLib,Swaption *ptSwpt,Volatility *ptVol)
{
  int N,alpha, beta,i,j,k,l;
  double price;
  MartingaleV *ptV;
  RandomGenerator *ptRand;
  double tenor;
  double sumbond;
  double Vab,SumVa,B1, SumVi;
  double t, et, p;
  
  tenor=ptLib->tenor;
  t=ptSwpt->swaptionMaturity;
  N=ptLib->numberOfMaturities;
  alpha=(int)(ptSwpt->swaptionMaturity/tenor);
  beta=(int)(ptSwpt->swapMaturity/tenor);
 
  mallocMartingaleV(&ptV, ptLib);
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors);
  et=0;
  price=0.0;
  // compute bond value B(0,T1)
  B1=1/(1.+tenor*ptLib->libor[0]);
  
  initMartingaleV(ptV);
  
  for(l=0;l<Nbmontecarlo;l++)
    {
	  initMartingaleV(ptV);
	  drawMartingaleV(ptV,ptVol,ptRand,dt,0,t);
	  /*drawing martingaleV time swaption first payment t*/
	  Vab=0;
	  for(i=alpha;i<beta;i++)
		{
		  Vab+=ptV->Vvalue_t[i]; 
		}
	  SumVi=0;
	  for(j=N;j>beta;j--){SumVi+=ptV->Vvalue_t[j];}
	  sumbond=0;
	 
	  for(i=beta;i>alpha;i--)
		{
		  SumVi+=ptV->Vvalue_t[i];
		  sumbond+=SumVi;
		}
	  p=maximum(Vab - tenor*ptSwpt->strike*sumbond,0.0);
	  price+= p;
	  et+=p*p;
    }

  price=price/Nbmontecarlo;
  et=sqrt((et/Nbmontecarlo - price*price)/Nbmontecarlo); /*et=ecart type for information*/
  price=price*B1;
 

  //printSwaption(ptSwpt);
  //printf("ecart type: %lf\n",et*10000);
  
  // printf("Swaption (simulated under the termimal measure thanks to a discrete martingale term): %lf\n",price*10000);
  freeRandomGenerator(&ptRand);
  freeMartingaleV(&ptV);

  return price;
}

int computeCapletSpotV(int Nbmontecarlo,double dt,Libor *ptLib, double K, Volatility *ptVol, double *Caps)
{
  int i,l,j,k;
  int N;
  double B1, Dk;
  MartingaleV *ptV;
  RandomGenerator *ptRand;
  double tenor;
  
  tenor=ptLib->tenor;
  B1=1/(1.+tenor*ptLib->libor[0]);
  N=ptLib->numberOfMaturities;

  mallocMartingaleV(&ptV, ptLib);
  mallocRandomGenerator(&ptRand,ptVol->numberOfFactors); 
  
  // compute bond value P(0,SwaptionMat)
  Dk=1.;
  
  for(j=0; j<=N; j++){Caps[j]=0;}
  
  for(l=0;l<Nbmontecarlo;l++)
    {
	  initMartingaleV(ptV);
	  drawMartingaleV(ptV,ptVol,ptRand,dt,0,tenor);
	  /*drawing the martingaleV till first maturity T1*/
	  for(k=1;k<N-1;k++)
		{
		  drawMartingaleV(ptV,ptVol,ptRand,dt,k*tenor, (k+1)*tenor);
		  /*Drawing successively the martingaleV asset till maturity T[k+1]*/
		  
		  SetLiborWithMartingaleV(ptLib, ptV);
		  /*thus ptLib are the drawing libors at time maturity T[k+1]*/ 

		  Dk=1;/*actualisation factor*/
		  for (i=1;i<=k;i++)
			{
			  Dk*=(1.+tenor*ptLib->libor[i]);
			}
		  Caps[k]+=maximum(ptLib->libor[k]-K, 0)/Dk;
		  /*computation of the caplets prices from libor L[k] at time T[k+1] (equal to Libor[k] at time T[k])*/
		}
	  drawMartingaleV(ptV,ptVol,ptRand,dt,(N-1)*tenor,(N)*tenor);
	  SetLiborWithMartingaleV(ptLib, ptV);
	  Dk=1;
		  for (i=1;i<N;i++)
			{
			  Dk*=(1.+tenor*ptLib->libor[i]);
			}
	  Caps[N-1]+=maximum(ptLib->libor[N-1]-K, 0)/Dk;
	  /*Last caplet*/
	}  
  for(j=1; j<N; j++){Caps[j]=Caps[j]*B1/Nbmontecarlo;}

  //printdoubles("Caplets.dat",N,Caps);
  freeRandomGenerator(&ptRand);
  freeMartingaleV(&ptV);
 
  return 1;

}
/*
  Non teste -> exclue du code
int computeBiasSpotV(int Nbmontecarlo,double dt,Libor *ptLib1,double K, Volatility *ptVol, double *Bias)
{
  int i,j,k,l,m;
  int N=ptLib1->numberOfMaturities;
  int Nfac=ptVol->numberOfFactors;
  RandomGenerator *ptW;
  MartingaleV *ptV;
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
  mallocMartingaleV(&ptV, ptLib1);
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
	  initMartingaleV(ptV);
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
			  Di+=ptV->Vvalue_t[i+1]*delta;

			  T=ptV->maturity[i+1];
			  if(t<T){ for(j=0; j<Nfac; j++){ci[j]=ci[j]+ptV->Vvalue_t[i+1]*delta*evalVolatility(ptVol,j,t,T)/Di;}}
			  T=ptV->maturity[i];
			  if(t<T){  for(j=0; j<Nfac; j++){sigma[i][j]=evalVolatility(ptVol,j,t,T) + ci[j];}}
			  else { for(j=0; j<Nfac; j++){sigma[i][j]=ci[j];}}
			
			  ptV->Vvalue_t[i]=ptV->Vvalue_t[i]*exp(-0.5*ps(Nfac, sigma[i], sigma[i])*dt + sqrt(dt)*ps(Nfac, sigma[i], ptW->val));
			}
		  T=ptV->maturity[N-1];
		  for(j=0; j<Nfac; j++){sigma[N-1][j]=evalVolatility(ptVol,j,t,T);}
		  ptV->Vvalue_t[N-1]=ptV->Vvalue_t[N-1]*exp(-0.5*ps(Nfac,sigma[N-1],sigma[N-1])*dt+sqrt(dt)*ps(Nfac,sigma[N-1],ptW->val));

		  if(t+dt==ptLib2->maturity[m])
			{
			  SetLiborWithMartingaleV(ptLib2,ptV);
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
  freeMartingaleV(&ptV);
  freeLibor(&ptLib2);
  free(cap1);
  free(cap2);
  free(ci);
  free(libor2);
  for(i=0;i<N;i++){free(sigma[i]);}
  free(sigma);
  
  return 1;
}
*/

int computeZCBondSpotV(int Nbmontecarlo,double dt,Libor *ptLib, double maturity, Volatility *ptVol)
{
  int i,j,l;
  int N=ptLib->numberOfMaturities;
  double BN;
  double delta=ptLib->tenor;
  double *Di;
  double *CumulDi;
  double *B0;
  MartingaleV *ptV;
  RandomGenerator *ptW;

  initLibor(ptLib);
  mallocMartingaleV(&ptV, ptLib);
  mallocRandomGenerator(&ptW,ptVol->numberOfFactors);
  Di=(double*)malloc(N*sizeof(double));
  CumulDi=(double*)malloc(N*sizeof(double));
  B0=(double*)malloc((N+1)*sizeof(double));
 
  for(i=0;i<N;i++){CumulDi[i]=0;}
  for(i=0;i<=N;i++){B0[i]=1;}
  for(i=1;i<=N;i++){for(j=0;j<i;j++){B0[i]=B0[i]/(1+delta*ptLib->libor[j]);}}
  B0[0]=1;


  for(i=0;i<N;i++){CumulDi[i]=0;}
  
  for(l=0;l<Nbmontecarlo;l++)
    {
	  initMartingaleV(ptV);
	 
	  for(i=1;i<N;i++)
		{
		  Di[i]=0;
		  drawMartingaleV(ptV,ptVol,ptW,dt,(i-1)*delta, delta);	  
		  for(j=i; j<=N; j++){Di[i]+=ptV->Vvalue_t[j];}
		}
	  for(i=0;i<N;i++){CumulDi[i]+=Di[i];}
	  
	}
 
   
  for(i=1;i<N;i++){CumulDi[i]=CumulDi[i]/Nbmontecarlo*B0[1];}
  CumulDi[0]=1;
  for(i=0;i<N;i++){B0[i]=(B0[i]-CumulDi[i])/B0[i];}
  B0[N]=0;
  printdoubles("zcbonderror.dat",N+1,B0);

  free(Di);
  free(CumulDi);
  freeMartingaleV(&ptV);
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
