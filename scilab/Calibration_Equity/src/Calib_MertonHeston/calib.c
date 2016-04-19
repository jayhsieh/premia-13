#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

#include "optype.h"


#include "parameter.h"
#include "./svj_std/svj.h"
#include "./svj_std/grad_svj.h"
#include "costFunction.h"
#include "gradFunction.h"
#include "paramsbfgs.h"
//#include "./utils/plot.h"
#include "./utils/datamarket.h"
#include "./utils/minix.h"
#include "./utils/utils1.h"
#include "bfgsb.h"

//
#define dimx1 5
#define dimx2 4
#define dimx3 8
//#define dimx  8
#define nbplot 20



FILE *ftest;

int main(int argc , char ** argv)
{
  int i;
  double St0,r,lambda,divid,kappa,Kmin,Kmax,Tmin,Tmax;
  double m0,v;
  double tmp;
  //  PricingMethod mert;
  int dimx;
  
  int nbd[dimx3],nbd0[dimx3],TypeNorme,TypeModel,TypeModelPrix;
  double fx,grad[dimx3],xmin[dimx3],xmax[dimx3],xmin0[dimx3],xmax0[dimx3];
  double x0[dimx3],sol[dimx3];
  //double x1[5]={.04, .2, 0.02, 0.1, -0.1},x2[4]={.04, .9, 0.1, 0.1},x3[8]={.04, .9, 0.02, 0.1, 0.1, .9, 0.1, 0.1};
  //  double xtoto[dimx3]={0.16, 2., 0.09, 0.40, -0.8, 2.763070, 0.225775, 0.833586,};
	//[dimx3]={0.786537, 0.615608, 0.447587, 0.226881, -0.624934, 2.763070, 0.225775, 0.833586,};
	//{0.009294, 0.009294, 0.009188, 0.173803, -0.479010, 0.508971, 0.101312, -0.199892};
  double theta,sigma0,V0,sigmav,rho;
  double normg;
  double grad0[dimx3],minix[dimx3],minixmin[dimx3],minixmax[dimx3],x00[dimx3];
  int dimminix,mininbd[dimx3],optim,TypeNormeMinix;
  int LogNorme=0,MelangeNormes=1;
  char prefix[20],date[20];
  char nom[100];
  int nbdata1=11,nbdata2=5,nbdata=nbdata1*nbdata2;
  int nbetapes;
  int typen,norme2,norme3,lequel=1;
  time_t today;
  struct tm *tb;
  DataMarket *DM;
  double xinit[5];
  int synthetique = 1;
  // Fin declarations
  //
  if(synthetique==1){
    //=============================================================================
    // Pour generer des donnes synthetiques
    //=============================================================================
    // 
    TypeModelPrix = 2;
    //
    St0     = 100.;
    Kmin    = 0.5*St0;
    Kmax    = 1.5*St0;
    Tmin    = .1;
    Tmax    = 1.;
    r       = 5./100;    
    r       = log(1+r);
    divid   = 0.;
    //
    sigma0  = 0.3;
    V0      = sigma0*sigma0;
    // Pour Merton
    lambda  = 1.;
    m0      = -0.1; 
    v       = 0.2;
    // Pour Heston
    kappa   = 2.; 
    theta   = 0.04; 
    sigmav  = 0.3;
    rho     = -0.5;
    //
    init_PrixObs(nbdata1,nbdata2,TypeModelPrix,Kmin,Kmax,Tmin,Tmax,St0,r,divid,V0,kappa,theta,sigmav,rho,lambda,m0,v);
    //
    //=============================================================================
    //
    nbdata = nbdata1*nbdata2;
    //
  }
  //
  //===============================================================================
  // lecture des parametres de la calibration : type de norme, nb etapes, nb data, et xinit.
  // 
  read_params(TypeModelPrix,&typen,&nbetapes,&nbdata,xinit);
  //
  // Creation des donnees Marches
  DM = CreateDataMarket(nbdata);
  // Lecture a partir du fichier
  readDataMarket(DM);
  // Affichage
  AfficheDataMarket(DM);
  //
  //===============================================================================
  // Le fichier de sortie
  //
  if(TypeModelPrix==1) 
    sprintf(prefix,"CalcHeston");
  else if(TypeModelPrix==2) 
    sprintf(prefix,"CalcMerton");
  else
    sprintf(prefix,"CalcHestMert");
  //
  today = time(NULL);
  tb = localtime(&today);   
  convert_date(tb->tm_mday,tb->tm_mon,tb->tm_year,date);
  //  sprintf(nom,"%s_%s_results-tests-1.dat",prefix,date);
  sprintf(nom,"%s_results.dat",prefix);
  ftest = fopen (nom, "wt");
  //
  //===============================================================================
  // Initialisation des parametres pour la fonction costfunction
  TypeModel = TypeModelPrix;
  TypeNorme = 3;
  LogNorme  = 0;
  initParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
  initParamsGrad(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
  //
  //
  /*
  // Initialisation de la solution
  init_sol(TypeModel,&dimx,x,dimx1,dimx2,dimx3,V0,kappa,theta,sigmav,rho,lambda,m0,v);
  //
  for(i=0;i<dimx;i++) sol[i] = x[i];
  //
  */
  //
  //================================================================================
  //
  //
  if(TypeModel==1) 
    dimx = 5;
  else if(TypeModel==2)
    dimx = 4;
  else
    dimx = 8;
  //
  // init des bornes
  init_bornes(TypeModel,dimx,nbd,xmin,xmax);
  for(i=0;i<dimx;i++) {
	xmin0[i] = xmin[i];
	xmax0[i] = xmax[i];
	nbd0[i]  = nbd[i];
  }
  //
  //
  //================================================================================
  //
  // On test xinit
  fx = costFunction(dimx,xinit);
  gradFunction(dimx,xinit,grad); 
  //
  printf("xinit  = ");
  for(i=0;i<dimx;i++) printf("%f ",xinit[i]); printf("\n");
  printf("fx = %f\n",fx);
  printf("g  = ");
  for(i=0;i<dimx;i++) printf("%f ",grad[i]);printf("\n"); 
  //
  AffichePrixDataMarket(DM);
  //
  //================================================================================
  //
  // Preparation de la calibration
  //
  TypeModel = TypeModelPrix;
  //  
  if(typen==33)
    {
      norme2 = 3;
      norme3 = 3;
    }
  else if(typen==22)
    {
      norme2 = 2;
      norme3 = 2;
    }
  else
    {
      norme2 = 2;
      norme3 = 3;
    }
  //
  TypeNormeMinix = norme3;
  TypeNorme = norme3;
  LogNorme = 0;
  //
  //==============================================================================
  //
  //
  //if(TypeModelPrix==3) 	  TypeModel = 2;
  init_bornes(TypeModel,dimx,nbd,xmin,xmax);
  init_bornes(TypeModel,dimx,nbd0,xmin0,xmax0);
  //
  for(i=0;i<dimx;i++)  x0[i] = xinit[i];
  for(i=0;i<dimx;i++)  x00[i] = x0[i];
  //
  fprintf(ftest,"x0  = ");
  for(i=0;i<dimx;i++) fprintf(ftest,"%f, ",x0[i]); fprintf(ftest,"\n");
  //
  ReInitParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
  //
  optim = 0;
  normg = 1.;
  //
  gradFunction(dimx,x0,grad0);
  for(i=0;i<dimx;i++)  if(grad0[i]==0.)  grad0[i] = 1.e-8;
  normg = norme(dimx,grad0);
  //
  //
  while(optim<nbetapes && normg >1.e-6)
    {
      //
      fflush(NULL);  
      optim++;
      printf("============================================\n");
      printf(" etape %d \n",optim);
      printf("============================================\n");
      printf("x0  = ");
      for(i=0;i<dimx;i++) printf("%f, ",x0[i]); printf("\n");
      //
      initbfgsb(0);
      ReInitBornes(dimx,nbd0,xmin0,xmax0,nbd,xmin,xmax);
      // 
      // D'abord la norme 3 
      //
      TypeNorme = norme3;
      ReInitParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
      //
      bfgsb(dimx,x0,nbd,xmin,xmax);
      //
      fx = costFunction(dimx,x0);
      gradFunction(dimx,x0,grad);
      printf("Apres TypeNorme = %d \n",TypeNorme);
      affiche_sol(dimx,sol,x0,fx,grad,grad0);
      //
      // On fait une recherche dans chaque direction.
      //
      //
      for(i=0;i<dimx;i++)
	{
 	  TypeNorme = TypeNormeMinix;
	  ReInitParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
	  //
	  BloqueVariable_i(dimx,i,nbd,x0,xmin,xmax);
	  //
	  dimminix = initminix(dimx,nbd,xmin,xmax,mininbd,minixmin,minixmax);
	  x2minix(x0,minix);
	  initbfgsb(1);
	  printf("Avant minix %f \n",minix[0]);
	  bfgsb(dimminix,minix,mininbd,minixmin,minixmax);
	  minix2x(x0,minix);
	  FreeMinix();
	  ReInitBornes(dimx,nbd0,xmin0,xmax0,nbd,xmin,xmax);
	  printf("Apres minix %f \n",minix[0]);
        }
      //
      // 
      // Apres on bloque certaines variables
      //
      TypeNorme = TypeNormeMinix;
      ReInitParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
      //
      if(TypeModel==3)
	{
	  // on bascule netre Merton et Heston une fois sur deux. 
	  lequel *=-1;
	}
      BloqueVariables(dimx,TypeModel,lequel,nbd,x0,xmin,xmax);
      //
      dimminix = initminix(dimx,nbd,xmin,xmax,mininbd,minixmin,minixmax);
      x2minix(x0,minix);
      initbfgsb(1);
      bfgsb(dimminix,minix,mininbd,minixmin,minixmax);
      minix2x(x0,minix);
      FreeMinix();
      ReInitBornes(dimx,nbd0,xmin0,xmax0,nbd,xmin,xmax);
      //
      //
      tmp = fabs(grad[0]/grad0[0]);
      for(i=0;i<dimx;i++) if(fabs(grad[i]/grad0[i]) > tmp)
	tmp = fabs(grad[i]/grad0[i]);
      for(i=0;i<dimx;i++)  if(fabs(grad[i]/grad0[i])/tmp < 0.1)  
	{
	  printf(" On bloque la variable %d \n",i);
	  xmin[i] = x0[i];
	  xmax[i] = xmin[i];
	  nbd[i] = 2;
	}
      TypeNorme = TypeNormeMinix;
      ReInitParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
      //
      dimminix = initminix(dimx,nbd,xmin,xmax,mininbd,minixmin,minixmax);
      x2minix(x0,minix);
      initbfgsb(1);
      bfgsb(dimminix,minix,mininbd,minixmin,minixmax);
      minix2x(x0,minix);
      FreeMinix();
      ReInitBornes(dimx,nbd0,xmin0,xmax0,nbd,xmin,xmax);
      //
      //
      // Ensuite la norme 2
      //
      initbfgsb(0);
      TypeNorme = norme2;
      ReInitParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
      //
      bfgsb(dimx,x0,nbd,xmin,xmax);
      //
      fx = costFunction(dimx,x0);
      gradFunction(dimx,x0,grad);
      printf("Apres TypeNorme = %d \n",TypeNorme);
      affiche_sol(dimx,sol,x0,fx,grad,grad0);
      //
      // Test d'arret
      //
      TypeNorme = norme3;
      LogNorme = 0;
      ReInitParams(nbdata,DM,TypeNorme,TypeModel,LogNorme,MelangeNormes);
      //
      fx = costFunction(dimx,x0);
      gradFunction(dimx,x0,grad);
      normg = norme(dimx,grad);
      //
      if(fx<1.) LogNorme = 1; else LogNorme = 0;
      //
      affiche_sol(dimx,sol,x0,fx,grad,grad0);
      AffichePrixDataMarket(DM);
      //
      // A VOIR
      for(i=0;i<dimx;i++) 
	{
	  x00[i] = x0[i];
	  grad0[i] = grad[i];
	  if(grad0[i]==0.)  grad0[i] = 1.e-8;
	}     
      //
    }
  //
  faffiche_sol(ftest,dimx,sol,x0,fx,grad,grad0);
  affiche_sol(dimx,sol,x0,fx,grad,grad0);
  fclose(ftest);
  //checksigmaimp();
  FreeDataMarket(DM);
  //
  return 0;
  //
}

