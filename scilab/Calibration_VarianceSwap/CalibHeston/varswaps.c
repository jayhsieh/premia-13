#include <math.h>
#include "varswaps.h"
#include <stdlib.h>

static int initVSOk = 0;
int TypeVSNorme, TypeVSModel;

static VSMarket *VS;

double kappaFunc(double x, int N, double A, double B, double C, double D1, double D2, double D3, double G);
int calcCoef(double kappa, double *A, double *B, double *C, double *D1, double *D2, double *D3, double *G);

//===================================== VS FIT BY EQUATION =============
// Fit kappa, theta and V0 to Variance Swaps solving a system of (nonlinear) equations

int fitEquationVS(int dimx, double *x)
{

  int ifOK, N;
  double eps;
  double V0, kappa,theta;
  double A, B, C, D1, D2, D3, G, det0;
  double dk, kl, kr, kc, fkl, fkr, fkc, sol;
    
	ifOK=0;
	eps=1e-4;
	dk = 0.1;

	  if(initVSOk != 1) {
		printf("Error! Parameters of Variance Swap were not initialized!\n ");	
		exit(-1);
	  }
	N = VS->nbdata;
	if(N == 0) {ifOK=0; return ifOK;}
 
	  //                 FIND KAPPA
	kl=0.0; kr=15.0;
	calcCoef(kl, &A, &B, &C, &D1, &D2, &D3, &G);
	fkl = kappaFunc(kl, N, A, B, C, D1, D2, D3, G); 
	calcCoef(kr, &A, &B, &C, &D1, &D2, &D3, &G);
	fkr = kappaFunc(kr, N, A, B, C, D1, D2, D3, G);

	while( (fkl*fkr>0.0)&&(kl<kr) )
	{
		kl += dk;
		calcCoef(kl, &A, &B, &C, &D1, &D2, &D3, &G);
		fkl = kappaFunc(kl,N, A, B, C, D1, D2, D3, G); 
	} 
	
	while (fabs(kr-kl) > eps) 
	{
	  kc  = (kl+kr)/2.;
	  calcCoef(kc, &A, &B, &C, &D1, &D2, &D3, &G);
	  fkc = kappaFunc(kc, N, A, B, C, D1, D2, D3, G);
	  if(fkc == 0.) {
		sol = kc;
		ifOK =  1;
		break;
	  }
	  
	  if (fkc*fkr < 0. )
		{kl = kc; fkl = fkc;}
	  else
		{kr = kc; fkr = fkc;}
	}
  //
 
                                 // FIND THETA AND V0
  kappa = (kl+kr)/2.;
  det0 = N*B - C*A;
  
  while( (fabs(det0)<1e-06)&&(kappa>0.) ) 
  {
	kappa -= dk;
	calcCoef(kappa, &A, &B, &C, &D1, &D2, &D3, &G);
	det0 = N*B - C*A;
  }
  
  theta = (D1*B - D2*A)/det0;
  V0 = (N*D2 - C*D1)/det0 + theta;

  ifOK = (kappa>0.)&&(theta>0.)&&(theta<3.0)&&(V0>0.)&&(V0<3.0);

printf("Solving equation:\n kappa= %f  theta= %f  V0= %f  ifOK=%d\n", kappa, theta, V0, ifOK);

  if(ifOK)
  {
	  x[0] = V0;
	  x[1] = kappa;
	  x[2] = theta;
  }
	return ifOK;
}
//---------------------------------------------
int calcCoef(double kappa, double *A, double *B, double *C, double *D1, double *D2, double *D3, double *G)
{
	int i;
	double kt, ekt, expkt;
	double T, pobs;

	*A=0;
	*B=0;
	*C=0;
	*D1=0;
	*D2=0;
	*D3=0;
	*G=0;

	if (kappa<0) {return 0;}

	for(i=0; i<VS->nbdata; i++)
	{
		T   = VS->T[i];
		kt = kappa*T;
		ekt = exp(-kt);
		if(kt==0) { expkt = 1;}
		else { expkt = (1.0-ekt)/kt; } 
		pobs = VS->PriceObs[i];
		*A += expkt;
		*B += expkt*ekt;
		*C += ekt;
		*D1 += pobs;
		*D2 += pobs*ekt;
		*D3 += pobs*expkt;
		*G += expkt*expkt;
	}
	return 1;
}
//--------------------------------------------
double kappaFunc(double x, int N, double A, double B, double C, double D1, double D2, double D3, double G)
{
	return (D1*B - D2*A)*A + (N*D2 - C*D1)*G - (N*B - C*A)*D3;
}
//===================================== VS COST FUNC ==================
// Cost function, to be minimized to fit kappa, theta and V0 to Variance Swaps 

double costFunctionVS(int dimx, double *x)
{
  
  int i;
  double fx=0,fxi;
  
  double V0, kappa,theta, sigmav, rho, prix;
  double T, prixobs;
    
	  if(initVSOk != 1) {
		printf("Error! Parameters of Variance Swap were not initialized!\n ");	
		exit(-1);
	  }
	  
	  V0     = x[0];
	  kappa  = x[1];
	  theta  = x[2];
	  sigmav = x[3];
	  rho    = x[4];

	  for(i=0;i<VS->nbdata;i++)
	    {
	      
	      T   = VS->T[i];
	      
	      prix = VS_Heston(T,V0,kappa,theta); 
	      
	      prixobs = VS->PriceObs[i];
	      if (TypeVSNorme == 1)  
		{ 
			fxi = prix-prixobs;
			fxi = fxi*fxi;
		}
		else if(TypeVSNorme == 4)
		{
			fxi = ( prix - prixobs )*log(prix/prixobs);
		}
		else /*if(TypeNorme == 2)*/ 
		{
			 fxi = (prix-prixobs)/prixobs;
			 fxi = fxi*fxi;
		} 
	      
	      VS->PriceMod[i] = prix;
	      
	      fx = fx + fxi*VS->wi[i];
	}
	return fx;
}
//

//
void FreeCostFunctionVS()
{
  FreeVSMarket(VS);
  VS = NULL;  
}

//======================================VS MARKET=============================

void initParamsVS(int nbdata,VSMarket *_VS, int _TypeNorme, int _TypeModel)
{
  int i;
  //
  if(VS == NULL) 
    VS = CreateVSMarket(nbdata);
  //
  for(i=0;i<nbdata;i++)
	{
	  VS->nbdata     = _VS->nbdata;
	  VS->T[i]       = _VS->T[i];
	  VS->wi[i]      = _VS->wi[i];
	  VS->PriceObs[i] = _VS->PriceObs[i];
	}
  //
  TypeVSNorme     = _TypeNorme;
  TypeVSModel     = _TypeModel;
  //
  initVSOk        = 1;
  // 
}

//----------------------------------------------------

VSMarket * CreateVSMarket(int N)
{
  int i;
  
  VSMarket * VS;
  VS          = (VSMarket*)malloc(sizeof(VSMarket));
  VS->nbdata  = N;
  VS->T       = (double*)malloc(sizeof(double)*N);
  VS->wi      = (double*)malloc(sizeof(double)*N);
  VS->PriceObs = (double*)malloc(sizeof(double)*N);
  VS->PriceMod = (double*)malloc(sizeof(double)*N);

  for(i=0;i<N;i++) VS->wi[i] = 1.;
	
  return VS;
  
}

//-------------------------------------------------------

void readVSMarket(VSMarket *VS)
{
  int i;
  double T, price; 
  FILE *fichier;
  //
  fichier = fopen ("VSMarket.txt", "r");
  //
  //fscanf(fichier,"%d \n", &type);
  //TypeVSNorme = type;

  for(i=0;i<VS->nbdata;i++)
	{
	fscanf(fichier,"%lf %lf \n", &T, &price);
	VS->T[i]   = T;
	price /= 100.0;
	VS->PriceObs[i]   = price*price;
	}
  fclose(fichier);
  //
}

//------------------------------------------------------
void CreateSyntheticVSMarket(VSMarket *VS, int nt, double *TS, double *x)
{
  int i;
  //
  if(VS == NULL) 
    VS = CreateVSMarket(nt);
  //
  VS->nbdata     = nt;
  for(i=0; i<nt; i++)
	{
	  VS->T[i]       = TS[i];
	  VS->PriceObs[i] = VS_Heston(TS[i], x[0], x[1], x[2]);
	}
}
//------------------------------------------------------

void PrintVSMarket(VSMarket *VS)
{
  int i;
  
  printf("VSMarket : %d  data\n", VS->nbdata);
  printf("          T,              PriceObs,        \n");
  for(i=0; i<VS->nbdata; i++)
	{
	  printf("         %f         %f\n", VS->T[i], VS->PriceObs[i]);
	}
}

//------------------------------------------------------

void FreeVSMarket(VSMarket *VS)
{
  //
  free(VS->T);
  free(VS->wi);
  free(VS->PriceObs);
  free(VS->PriceMod);
  free(VS);
  //
  VS = NULL;
  //
}

//===================================PRINT RESULTS=================================

void PrintVSPrices()
{
  int i;
  double p1,p2,e1;
  printf("VarianceSwaps : %d  data\n",VS->nbdata);
  printf("           T,            PriceObs,          PriceMod,           Erreur (%%) \n");
  for(i=0; i<VS->nbdata; i++)
  {
	p1 = sqrt(VS->PriceObs[i])*100.0;
	p2 = sqrt(VS->PriceMod[i])*100.0;
	e1 = 100.*fabs(p1-p2)/p1;
	printf("        %f       %f         %f         %f       \n", VS->T[i], p1, p2, e1); 
  }  
}
//
//----------------------------------------------------------
void PrintVSMarketFile()
{
  int i;
  double p1,p2,e1;
  FILE *fres;
  
  fres = fopen ("CalibResVS.txt", "wt");

  fprintf(fres, "VarianceSwaps : %d  data\n",VS->nbdata);
  fprintf(fres, "           T,            PriceObs,          PriceMod,           Erreur (%%) \n");
  for(i=0;i<VS->nbdata;i++)
  {
	p1 = sqrt(VS->PriceObs[i])*100.0;
	p2 = sqrt(VS->PriceMod[i])*100.0;
	e1 = 100.*fabs(p1-p2)/p1;
	fprintf(fres, "        %f       %f         %f         %f       \n", VS->T[i], p1, p2, e1); 
  }  
 fclose(fres);
}

//========================================VARIANCE SWAP===================

double VS_Heston(double T, double V0, double kappa, double theta)
{
	if(kappa==0.0) { return V0;}
	else {	return theta + (V0 - theta)*(1.0 - exp (-kappa*T) )/(kappa*T);}
}
//==========================================VS=HESTON=GRAD==================

void  gradFunctionVS(int dimx,double *x, double *grad)
{
  int i;
  double V0, kappa, theta, price, price2; 
  double T, priceobs;
  
  double kt, ekt;
 
  for(i=0; i<dimx; i++){ grad[i] =0.;}
   
	V0     = x[0];
	kappa  = x[1];
	theta  = x[2];
	
	  for(i=0;i<VS->nbdata;i++)
	    {
		T   = VS->T[i];
		price = VS_Heston(T,V0,kappa,theta);
		priceobs =VS->PriceObs[i];
		if(TypeVSNorme == 1) 
			{price2 = 2.0*(price - priceobs);}
		else if(TypeVSNorme == 4)
		{
			price2 = (price - priceobs)/price + log(price/priceobs);
		}
		else /*(TypeVSNorme == 2)*/
			{price2 = 2.0*(price - priceobs)/(priceobs*priceobs);}
		kt = kappa*T;
		if(kt==0.0) { ekt = 1.0; }
		else { ekt = (1.0 - exp(-kt))/kt; }

		grad[0] += ekt*price2;
		grad[1] += (V0-theta)/kappa * (exp(-kt) - ekt) * price2;
		grad[2] += (1.0 - ekt)*price2;
	    }
}
