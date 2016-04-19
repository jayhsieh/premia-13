
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "utils/intg.h"
#include "pnl/pnl_mathtools.h"

#define dimx3 8


#include "utils1.h"
#include "costFunction.h"
#include "gradFunction.h"
#include "varswaps.h"

extern int optim_scilab( int n, double x[]); 


int TypeModel;
int ifChangeVar;
int ifFixVSParam;

void setScheme(int *TypeModel, int *ifSynth, int *ifChangeVar, int *ifFixVSParam);
int dialog(int TypeModel, int ifSynth, int ifChangeVar, int ifFixVSParam, int typen, int nbdata, int nvsdata, double *x);
int changedialog(int *TypeModel, int *ifSynth, int *ifChangeVar, int *ifFixVSParam, int *typen, int *nbdata, int *nvsdata, double *x);


int main()
{

  int i;

  int dimx;
  int ifSynth, nt, nk;
  double *KS, *TS, *rs, *ds;
  int *TypeOptS;

  int TypeNorme;
  double fx,grad[dimx3];
  double x0[dimx3],sol[dimx3];

  double grad0[dimx3];

  int LogNorme=0,MelangeNormes=1;

  char nom[100];

  int nvsdata;
  int nbdata;
  int typen;

  DataMarket *DM;
  VSMarket *VS;
  double xinit[5], xsynth[5];

  ifSynth = 2;
  nt=4; nk=17;

  FILE *ftest;

  //=========INPUT
  printf("readparams\n");

  dimx=5;

  read_params(&TypeModel, &typen, &nbdata, &nvsdata, xinit);

  setScheme(&TypeModel, &ifSynth, &ifChangeVar, &ifFixVSParam);

  if(dialog(TypeModel, ifSynth, ifChangeVar, ifFixVSParam, typen, nbdata, nvsdata, xinit))
    { 
      changedialog(&TypeModel, &ifSynth, &ifChangeVar, &ifFixVSParam, &typen, &nbdata, &nvsdata, xinit) ;
      //	return 0;
    } 

  if(typen==3){ifChangeVar = 1;}

  for(i=0;i<dimx3;i++)  {grad0[i] = 1.0; sol[i]=0.0;}

  if(ifSynth == 1) { nbdata = nt*nk; nvsdata = nt;}

  // INIT VS MARKET AND DATAMARKET
  printf("create vsmarket\n");
  VS = CreateVSMarket(nvsdata);
  printf("create datamarket\n");
  DM = CreateDataMarket(nbdata);
  printf("read data\n");

  if(ifSynth == 2) //READ DATA FROM FILE
    {
      readVSMarket(VS);
      readDataMarket(DM);
    }
  else // SYNTHETIC DATA
    { 
      TypeOptS = (int*)malloc(sizeof(int)*17);
      TS = (double*)malloc(sizeof(double)*3);
      KS = (double*)malloc(sizeof(double)*17);
      rs = (double*)malloc(sizeof(double)*3);
      ds = (double*)malloc(sizeof(double)*3);
      for(i=0; i<nt; i++)
        {
          TS[i] = 0.5 + i*0.5;
          rs[i] = 0.05;
          ds[i] = 0.0;
        }
      for(i=0; i<8; i++)
        {
          KS[i] = 10.0*(i+2);
          TypeOptS[i] = -1;
        }
      for(i=8; i<nk; i++)
        {
          KS[i] = 10.0*(i+2);
          TypeOptS[i] = 1;
        }
      xsynth[0] = 0.195565;
      xsynth[1] = 7.057543;
      xsynth[2] = 0.202996;
      xsynth[3] = 0.35;
      xsynth[4] = -0.4;
      CreateSyntheticVSMarket(VS, nt, TS, xsynth);
      CreateSyntheticDM(DM, nt, TS, rs, ds, 100.0, nk, KS, TypeOptS, dimx, xsynth);
    }
  TypeNorme=typen;
  initParamsVS(nvsdata,VS,TypeNorme,1);

  if(TypeModel==111) // FIT VARIANCE SWAP FIRST
    {

      printf("init params %d %d %d %d %f %f %f %f %f\n", TypeModel, typen,  nbdata, nvsdata, xinit[0], xinit[1], xinit[2], xinit[3], xinit[4]);

      fx = costFunctionVS(dimx,xinit);
      printf("costFunction %f\n", fx);
      gradFunctionVS(dimx,xinit,grad);
      printf("g  = ");
      for(i=0;i<dimx;i++) {printf("%f ",grad[i]);} printf("\n");

      for(i=0;i<dimx;i++)  x0[i] = xinit[i];

      if(!fitEquationVS(dimx, x0))
        { 
          optim_scilab(dimx,x0);
          if(!(x0[0]>0.)) x0[0] = xinit[0];
          if(!(x0[1]>0.)) x0[1] = xinit[1];
          if(!(x0[2]>0.)) x0[2] = xinit[2];
        }

      fx = costFunctionVS(dimx,x0);
      gradFunctionVS(dimx,x0,grad);

      affiche_sol(dimx,sol,x0,fx,grad,grad0);

      printf("\n print VS file\n");
      PrintVSMarketFile();

      for(i=0;i<dimx;i++)  xinit[i] = x0[i];

      if(ifFixVSParam) {TypeModel=11;}
      else {TypeModel=1;}
    }

  printf("init params %d %d %d %f %f %f %f %f\n", TypeModel, typen,  nbdata, xinit[0], xinit[1], xinit[2], xinit[3], xinit[4]);

  initParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
  initParamsGrad(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);

  fx = costFunction(dimx,xinit);
  printf("costFunction %f\n", fx);

  gradFunction(dimx,xinit,grad); 

  printf("xinit  = ");
  for(i=0;i<dimx;i++) printf("%f ",xinit[i]); printf("\n");
  printf("fx = %f\n",fx);
  printf("g  = ");
  for(i=0;i<dimx;i++) printf("%f ",grad[i]);printf("\n"); 
  /*Calibration*/

  if(ifChangeVar)  // CHANGE VARIABLE TO SET FINITE INTERVALS OF VALUES
    {
      printf("\nChange variables\n");
      xinit[0]=tan((xinit[0]-0.5)*M_PI);
      xinit[1]=tan((xinit[1]/30.0-0.5)*M_PI);
      xinit[2]=tan((xinit[2]-0.5)*M_PI);
      xinit[3]=tan((xinit[3]-0.5)*M_PI);
      xinit[4]=tan(xinit[4]*0.5*M_PI);
    }

  for(i=0;i<dimx;i++)  x0[i] = xinit[i];

  optim_scilab(dimx,x0) ;

  // OUTPUT
  if(ifChangeVar) // If variables were changed
    {	  printf("\n change variables back...\n");
      xinit[0] = 0.5 + atan(x0[0])/M_PI;
      xinit[1]  = 30.0*(0.5 + atan(x0[1])/M_PI);
      xinit[2]  = 0.5 + atan(x0[2])/M_PI;
      xinit[3] = 0.5 + atan(x0[3])/M_PI;
      xinit[4]    = 2.0*atan(x0[4])/M_PI;
    }
  else
    {
      for(i=0;i<dimx;i++)  xinit[i]=x0[i];
    }
  fx = costFunction(dimx,xinit);
  gradFunction(dimx,xinit,grad);

  affiche_sol(dimx,sol,xinit,fx,grad,grad0);

  sprintf(nom,"CalibResParams.dat");
  ftest = fopen (nom, "wt");
  faffiche_sol(ftest,dimx,xinit,fx,grad);
  fclose(ftest);
  printf("\nVS results\n");
  fx = costFunctionVS(dimx,xinit);
  gradFunctionVS(dimx,xinit,grad);

  affiche_sol(dimx,sol,xinit,fx,grad,grad0);

  printf("print VS file\n");
  PrintVSMarketFile();
  printf("ready\n");
  FreeVSMarket(VS);
  FreeCostFunctionVS();

  AfficheDataMarketFile();

  FreeDataMarket(DM);

  return 0;
}

void setScheme(int *TypeModel, int *ifSynth, int *ifChangeVar, int *ifFixVSParam)
{
  switch(*TypeModel)
    {
    case 1: *ifSynth = 2;
            *ifChangeVar =0;
            *ifFixVSParam = 0;
            *TypeModel = 1;
            break;
    case 11: *ifSynth = 2;
             *ifChangeVar =0;
             *ifFixVSParam = 1;
             *TypeModel = 111;
             break;
    case 111: *ifSynth = 2;
              *ifChangeVar =0;
              *ifFixVSParam = 0;
              *TypeModel = 111;
              break;
    case 2: *ifSynth = 2;
            *ifChangeVar = 1;
            *ifFixVSParam = 0;
            *TypeModel = 1;
            break;
    case 22: *ifSynth = 2;
             *ifChangeVar = 1;
             *ifFixVSParam = 1;
             *TypeModel = 111;
             break;
    case 222: *ifSynth = 2;
              *ifChangeVar =1;
              *ifFixVSParam = 0;
              *TypeModel = 111;
              break;
    case 3:  *ifSynth = 1;
             *ifChangeVar = 0;
             *ifFixVSParam = 0;
             *TypeModel = 1;
             break;
    case 33: *ifSynth = 1;
             *ifChangeVar = 0;
             *ifFixVSParam = 1;
             *TypeModel = 111;
             break;
    case 333: *ifSynth = 1;
              *ifChangeVar = 0;
              *ifFixVSParam = 0;
              *TypeModel = 111;
              break;
    case 4: *ifSynth = 1;
            *ifChangeVar = 1;
            *ifFixVSParam =  0;
            *TypeModel = 1;
            break;
    case 44: *ifSynth = 1;
             *ifChangeVar = 1;
             *ifFixVSParam = 1;
             *TypeModel = 111;
             break;
    case 444:*ifSynth = 1;
             *ifChangeVar = 1;
             *ifFixVSParam = 0;
             *TypeModel = 111;
             break;
    }
}

int dialog(int TypeModel, int ifSynth, int ifChangeVar, int ifFixVSParam, int typen, int nbdata, int nvsdata, double *x)
{
  int intans;
  //	double doublans;

  printf("Source of the data: ");
  if(ifSynth==1) {printf("Synthetic\n");}
  else {printf("From file\n");}

  printf("Use Variance Swaps?  ");
  if(TypeModel>10) {printf("Yes\n");}
  else {printf("No\n");}

  printf("Fix parameters after Variance Swaps fitting?  ");
  if(ifFixVSParam) {printf("Yes\n");}
  else {printf("No\n");}

  printf("Change parameters to fix its ranges?  ");
  if(ifChangeVar) {printf("Yes\n");}
  else {printf("No\n");}

  printf("Type of cost function:  ");
  switch(typen)
    {
    case 1:
      printf("Absolute errors in prices\n");
      break;
    case 2:
      printf("Relative errors in prices\n");
      break;
    case 3:
      printf("Absolute errors in implied volatilities\n");
      break;
    case 4:
      printf("Log-differences in prices\n");
      break;
    }	

  printf("Number of Options prices:  %d\n", nbdata);

  printf("Number of Variance Swaps prices:  %d\n", nvsdata);

  printf("Initial guess of parameters: \n");
  printf("current variance V0 =  %f\n", x[0]);
  printf("mean reversion Kappa =  %f\n", x[1]);
  printf("long-run variance Theta =  %f\n", x[2]);
  printf("volatility of volatility sigma =  %f\n", x[3]);
  printf("correlation Rho =  %f\n", x[4]);

  printf("\n Change settings? (0 - No, 1 - Yes) ");
  scanf("%d", &intans);

  return intans;
}

int changedialog(int *TypeModel, int *ifSynth, int *ifChangeVar, int *ifFixVSParam, int *typen, int *nbdata, int *nvsdata, double *x)
{
  // /*
  int intans;
  double doublans;

  printf("\nSource of the data? (1 - Synthetic, 2 - From file): ");
  scanf("%d", &intans);
  if(intans==1) {*ifSynth=1;}
  else {*ifSynth=2;}

  printf("\nUse Variance Swaps? (1 - Yes, 0 - No) ");
  scanf("%d", &intans);
  if(intans==1) {*TypeModel = 111;}
  else {*TypeModel = 1;}

  printf("\nFix parameters after Variance Swaps fitting? (1 - Yes, 0 - No) ");
  scanf("%d", &intans);
  if(intans==1) {*ifFixVSParam = 1;}
  else {*ifFixVSParam = 0;}

  printf("\nChange parameters to fix its ranges? (1 - Yes, 0 - No) ");
  scanf("%d", &intans);
  if(intans==1) {*ifChangeVar = 1;}
  else {*ifChangeVar = 0;}

  printf("\nType of cost function: \n ");
  printf("1 - Absolute errors in prices\n");
  printf("2 - Relative errors in prices\n");
  printf("3 - Absolute errors in implied volatilities\n");
  printf("4 - Log-differences in prices\n");

  printf("Your choice?");
  scanf("%d", &intans);
  if( (intans>0)&&(intans<5) ) {*typen = intans;}
  else {*typen = 1;}

  printf("\nNumber of Options prices:" );
  scanf("%d", nbdata);

  printf("\nNumber of Variance Swaps prices: ");
  scanf("%d", nvsdata);

  printf("\nInitial guess of parameters: \n");

  printf("\ncurrent variance V0 = ");
  scanf("%lf", &doublans );
  if(doublans>0) {x[0] = doublans;}

  printf("\n mean reversion Kappa = "); 
  scanf("%lf", &doublans );
  if(doublans>0) {x[1] = doublans;}

  printf("\nlong-run variance Theta = ");
  scanf("%lf", &doublans );
  if(doublans>0) x[2] = doublans;

  printf("\nvolatility of volatility sigma = ");
  scanf("%lf", &doublans );
  if(doublans>0) x[3] = doublans;

  printf("\ncorrelation Rho = ");
  scanf("%lf", &doublans );
  if( fabs(doublans)<1.0) x[4] = doublans;

  // */
  //	printf("Please make changes to file 'in.dat'\n");
  return 1;
}
