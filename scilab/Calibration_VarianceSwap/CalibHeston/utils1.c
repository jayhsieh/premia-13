#include <math.h>
#include <stdlib.h>
#include "utils1.h"
#include "bsvanillas.h"
#include "pnl/pnl_finance.h"

//============================= datamarket ===============================
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
  
  printf("DataMarket : %d  data lines\n",DM->nbdata);
  printf("    TypeOpt,      St0,             T,             K,             PrixObs,          SigmaImpObs(%%)   \n");
  for(i=0;i<DM->nbdata;i++)
	{
	  printf("        %d     %f       %f       %f        %f         %f         %f\n", DM->TypeOpt[i], DM->St0[i], DM->T[i], DM->r[i], DM->K[i], DM->PrixObs[i], 100.*DM->SigmaImpObs[i]);
	}
}

void PrintDataMarketFile(DataMarket *DM)
{
  int i;
  double p1,p2,e1,e2,sig1,sig2;
  FILE *fres;
  
  fres = fopen ("CalibRes.txt", "wt");
  
  fprintf(fres, "DataMarket : %d  data lines\n",DM->nbdata);
  fprintf(fres, "TypeOpt,    T,        K,          PrixObs,      PrixMod,        Erreur (%%)      SigmaImpObs,       SigmaImp       Erreur(%%)\n");
  for(i=0;i<DM->nbdata;i++)
  {
	p1 = DM->PrixObs[i];
	p2 = DM->PrixMod[i];
	sig1 = DM->SigmaImpObs[i];
	sig2 = DM->SigmaImp[i];
	e1 = 100.*fabs(p1-p2)/p1;
	e2 = 100.*(sig1 - sig2)/sig1;
	fprintf(fres, " %d      %f     %f     %e     %e     %f     %f    %f    %f\n", DM->TypeOpt[i],DM->T[i],DM->K[i],p1,p2,e1,sig1,sig2,e2); 
  } 

  fclose(fres);
}

void PrintDataMarket(DataMarket *DM)
{
  int i;
  double p1,p2,e1,e2,sig1,sig2;
  printf("DataMarket : %d  data lines\n",DM->nbdata);
  printf("    TypeOpt,      St0,             T,             K,             PrixObs,          PrixMod,           Erreur (%%)      SigmaImpObs,       SigmaImp       Erreur\n");
  for(i=0;i<DM->nbdata;i++)
  {
	p1 = DM->PrixObs[i];
	p2 = DM->PrixMod[i];
	sig1 = DM->SigmaImpObs[i];
	sig2 = DM->SigmaImp[i];
	e1 = 100.*fabs(p1-p2)/p1;
	e2 = 100.*(sig1 - sig2);
	printf("        %d     %f       %f       %f         %f         %f            %f         %f       %f        %f\n",DM->TypeOpt[i],DM->St0[i],DM->T[i],DM->K[i],p1,p2,e1,sig1,sig2,e2); 
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

//==================================================================================
// Affiche_Sol
//==================================================================================
void affiche_sol(int dimx,double *sol,double *x,double fx,double *grad,double *grad0)
{
  int i;
  //
  // sol
  //printf("sol  = ");
  //for(i=0;i<dimx;i++) printf("%f ",sol[i]);
  printf("\n");
  printf("x*  = ");
  for(i=0;i<dimx;i++) printf("%f, ",x[i]); printf("\n");
  printf("fx = %f\n",fx);
  printf("g  = ");
  for(i=0;i<dimx;i++) printf("%f ",grad[i]);printf("\n");
  /*printf("g0  = ");
  for(i=0;i<dimx;i++) printf("%f ",grad0[i]);printf("\n");
  printf("g/g0  = ");
  for(i=0;i<dimx;i++) printf("%f ",fabs(grad[i]/grad0[i]));
  printf("\n");*/
  //
  //
}
//==================================================================================
// F_Affiche_Sol
//==================================================================================
void faffiche_sol(FILE *ftest,int dimx,double *x,double fx,double *grad)
{
  int i;
  //
  // sol
  
  fprintf(ftest,"\n");
  fprintf(ftest,"x*  = ");
  for(i=0;i<dimx;i++) fprintf(ftest,"%f, ",x[i]); fprintf(ftest,"\n");
  fprintf(ftest,"fx = %f\n",fx);
  fprintf(ftest,"g  = ");
  for(i=0;i<dimx;i++) fprintf(ftest,"%f ",grad[i]);fprintf(ftest,"\n");
  
  fprintf(ftest,"\n");
  //
  //
}
//==================================================================================
// Init_Sol
//==================================================================================
void init_sol(int TypeModel,int *dimx,double *x,int dimx1,int dimx2,int dimx3,double V0,double kappa,double theta,double sigmav,double rho,double lambda, double m0,double v)
{
  
  if(TypeModel==1 || TypeModel==111)
	{
	  *dimx = dimx1;
	  x[0] = V0;
	  x[1] = kappa;
	  x[2] = theta;
	  x[3] = sigmav;
	  x[4] = rho;
	}
}
//==================================================================================
// Init_Bornes
//==================================================================================
void init_bornes(int TypeModel,int dimx,int *nbd,double *xmin,double *xmax)
{
  int i;
  
  for(i=0;i<dimx;i++) {nbd[i] = 2;}
  if(TypeModel==1 || TypeModel==11 || TypeModel==111)
	{
	  xmin[0] = 25.e-4;  xmax[0] = 0.9;  nbd[0] =2;
	  //	   	  	  xmin[0] = V0;  xmax[0] = V0;
	  xmin[1] = 0.01;  xmax[1] = 5.;  nbd[1] = 2;
	  xmin[2] = 25.e-4;  xmax[2] =0.9; nbd[2] = 2;
	  xmin[3] = 0.01;  xmax[3] = 0.9; nbd[3] = 2;
	  xmin[4] = -0.99;    xmax[4] = 0.99;
	}
  
}


//==================================================================================

void readDataMarket(DataMarket *DM)
{
  int i,j,type;
  double T,K,St0,r,d,wi,sigma,prix,delta;
  FILE *fichier;
  //
  fichier = fopen ("DataMarket.txt", "r");
  //
  for(i=0;i<DM->nbdata;i++)
	{
	  //
	  fscanf(fichier,"%d %lf %lf %lf %lf %lf %lf\n",&type,&T,&K,&sigma,&r,&d,&wi);
	  //
	  r = log(1.+r);
	  d = log(1.+d/100.);
	  //
	  St0 = 100.;
	  //
	  DM->TypeOpt[i]   = type;
	  DM->St0[i] = St0;
	  DM->T[i]   = T;
	  DM->K[i]   = St0*K;
	  DM->r[i]   = r;
	  DM->d[i]   = d;
	  DM->wi[i]   = wi;
	  DM->SigmaImpObs[i]   = sigma;
	  //
	  if(type==1)
		{j =  pnl_cf_call_bs(St0,St0*K,T,r,d,sigma,&prix,&delta);}
	  else
		{j =  pnl_cf_put_bs(St0,St0*K,T,r,d,sigma,&prix,&delta);}
	  DM->PrixObs[i]   = prix;
	  //
	}
  //
  fclose(fichier);
  //
}

//==================================================================================
void read_params(int *TypeModel,int *typen,int *nbdata, int *nvsdata, double *xinit)
{
  //
  FILE *fichier;
  //
  fichier = fopen ("in.dat","r");
  //
  fscanf(fichier,"%d \n", TypeModel);
  fscanf(fichier,"%d \n",typen);
  
  fscanf(fichier,"%d \n",nbdata);
  fscanf(fichier,"%d \n",nvsdata);
//  if((*TypeModel==1)||(*TypeModel==111))
//    {
      fscanf(fichier,"%lf %lf %lf %lf %lf\n",&xinit[0],&xinit[1],&xinit[2],&xinit[3],&xinit[4]);
//    }
  
  fclose(fichier);
  //
}
//==================================================================================     
