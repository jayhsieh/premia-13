#include "gradFunction.h"
//
static int initOk=0,TypeNorme,TypeModel;
static int LogNorme=0;
static int MelangeNormes=0;
static int sigmaimpok=0,sigmaimpko=0;
static DataMarket *DM;
static int fd = 0;
//
void  gradFunction(int dimx,double *x, double *grad)
{
  //
  int i,n;
  double fx=0,fxi;
  int ibidon;
  double V0,kappa,theta,sigmav,rho,lambda,m0,v,prix,delta;
  double St0,K,T,prixobs,sigma_imp_obs;
  double *gradi;
  NumFunc_1 *p;
  double ri,di;
  //
  p=(NumFunc_1*)malloc(sizeof(NumFunc_1));
  gradi = (double*)malloc(dimx*sizeof(double));
  //
  // initialisation du vecteur du gradient
  for(n=0;n<dimx;n++){gradi[n] =0.;  grad[n] =0.;}
  //
  // au cas ou le gradient est calcule par DF
  if(fd == 1)
    {
      gradFunctionFD(dimx,x,grad);
    }
  else
    {
      switch(TypeModel)
	{
	case -1:
	  // calcul du gradient de x^2 !!
	  for(i=0;i<dimx;i++)
	    {
	      grad[i] = 2. * x[i];
	    }
	  break;
	  // 
	  // calcul du gradient de x^2 par DF (au cas ou on ne sait pas faire !!)
	  // pour tester la DF
	case 0: 
	  gradFunctionFD(dimx,x,grad);
	  break;
	  //
	  // Call/Put Heston
	case 1:
	  if(initOk != 1) {
	    printf("Erreur, avant d'appeler costFunction, il faut initialiser les paramaetres.\n ");	
	    exit(-1);
	  }
	  //
	  V0     = x[0];
	  kappa  = x[1];
	  theta  = x[2];
	  sigmav = x[3];
	  rho    = x[4];
	  //
	  for(i=0;i<DM->nbdata;i++)
	    {
	      //
	      St0 = DM->St0[i];
	      T   = DM->T[i];
	      K   = DM->K[i];
	      p->Par[0].Val.V_DOUBLE = K;
	      ri =DM->r[i];
	      di =DM->d[i];
	      //
	      if(DM->TypeOpt[i] == 1) {
		ibidon=FT_Call_Heston(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,&prix,&delta);
		ibidon=Grad_FT_Call_Heston(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,dimx,gradi);
	      } else {
		ibidon=FT_Put_Heston(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,&prix,&delta);
		ibidon=Grad_FT_Put_Heston(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,dimx,gradi);
	      }
	      //
	      prixobs = DM->PrixObs[i];
	      sigma_imp_obs = DM->SigmaImpObs[i];
	      fxi = dNormeMarket(TypeNorme,DM->TypeOpt[i],St0,T,K,ri,di,prix,prixobs,sigma_imp_obs,dimx,gradi);
	      //
	      fx = fx + fxi*DM->wi[i];
	      for(n=0;n<dimx;n++) grad[n] = grad[n] + gradi[n]*DM->wi[i];
	      //  
	    }
	  break;
	  //
	  // Call/Put Merton
	case 2:
	  if(initOk != 1) {
	    printf("Erreur, avant d'appeler costFunction, il faut initialiser les paramaetres.\n ");	
	    exit(-1);
	  }
	  //
	  V0    = x[0];
	  lambda= x[1];
	  m0    = x[2];
	  v     = x[3];
	  //
	  for(i=0;i<DM->nbdata;i++)
	    {
	      //
	      St0 = DM->St0[i];
	      T   = DM->T[i];
	      K   = DM->K[i];
	      p->Par[0].Val.V_DOUBLE = K;
	      ri =DM->r[i];
	      di =DM->d[i];
	      //
	      if(DM->TypeOpt[i] == 1) {
		ibidon=FT_Call_Merton(St0,p,T,ri,di,V0,lambda,m0,v,&prix,&delta);
		ibidon=Grad_FT_Call_Merton(St0,p,T,ri,di,V0,lambda,m0,v,dimx,gradi);
	      } else {
		ibidon=FT_Put_Merton(St0,p,T,ri,di,V0,lambda,m0,v,&prix,&delta);
		ibidon=Grad_FT_Put_Merton(St0,p,T,ri,di,V0,lambda,m0,v,dimx,gradi);
	      }
	      //
	      prixobs = DM->PrixObs[i];
	      sigma_imp_obs = DM->SigmaImpObs[i];
	      fxi = dNormeMarket(TypeNorme,DM->TypeOpt[i],St0,T,K,ri,di,prix,prixobs,sigma_imp_obs,dimx,gradi);
	      //
	      fx = fx + fxi*DM->wi[i];
	      for(n=0;n<dimx;n++) grad[n] = grad[n] + gradi[n]*DM->wi[i];
	    }
	  break;
	  //
	  // Call/Put Heston+Merton
	case 3:
	  if(initOk != 1) {
	    printf("Erreur, avant d'appeler costFunction, il faut initialiser les paramaetres.\n ");	
	    exit(-1);
	  }
	  //
	  V0     = x[0];
	  kappa  = x[1];
	  theta  = x[2];
	  sigmav = x[3];
	  rho    = x[4];
	  lambda = x[5];
	  m0     = x[6];
	  v      = x[7];
	  //
	  for(i=0;i<DM->nbdata;i++)
	    {
	      //
	      St0 = DM->St0[i];
	      T   = DM->T[i];
	      K   = DM->K[i];
	      p->Par[0].Val.V_DOUBLE = K;
	      ri =DM->r[i];
	      di =DM->d[i];
	      //
	      if(DM->TypeOpt[i] == 1) {
		ibidon=FT_Call_HestMert(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,lambda,m0,v,&prix,&delta);
		ibidon=Grad_FT_Call_HestMert(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,lambda,m0,v,dimx,gradi);
	      } else {
		ibidon=FT_Put_HestMert(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,lambda,m0,v,&prix,&delta);
		ibidon=Grad_FT_Put_HestMert(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,lambda,m0,v,dimx,gradi);
	      }
	      //
	      prixobs = DM->PrixObs[i];
	      sigma_imp_obs = DM->SigmaImpObs[i];
	      fxi = dNormeMarket(TypeNorme,DM->TypeOpt[i],St0,T,K,ri,di,prix,prixobs,sigma_imp_obs,dimx,gradi);
	      //
	      fx = fx + fxi*DM->wi[i];
	      for(n=0;n<dimx;n++) grad[n] = grad[n] + gradi[n]*DM->wi[i];
	    }
	  break;
	}
    }
  //
  free(p);
  free(gradi);
  // dans le cas de la norme logarithmique : ln(fx)
  if(LogNorme==1)
    {
      if(fx>0.) {
	for(n=0;n<dimx;n++) grad[n] = grad[n]/fx ;
      } else {
	for(n=0;n<dimx;n++) grad[n] =  0.;
      }
    }
  //
}
// Calcul du gradient par differences finies
//
void  gradFunctionFD(int dimx,double *x, double *grad)
{
  int i;
  double fx,fxh,dx;
  //
  fx = costFunction(dimx,x);
  //
  for(i=0;i<dimx;i++)
    {
      dx      = x[i] * h0;
      // Pour  le cas x[i] = 0, on prends dx = h0
      if(dx==0.) dx = h0;
      //	  
      x[i]    = x[i] + dx;
      fxh     = costFunction(dimx,x);
      grad[i] = (fxh - fx) / dx;
      x[i]    = x[i] - dx;
    }
  //
}
//
void initParamsGrad(int nbdata,DataMarket *_DM,int _TypeNorme,int _TypeModel,int _LogNorme,int _MelangeNormes)
{
  //
  int i;
  //
  if(DM == NULL)
    DM = CreateDataMarket(nbdata);
  //
  for(i=0;i<nbdata;i++)
    {
      DM->nbdata     = _DM->nbdata;
      DM->TypeOpt[i] = _DM->TypeOpt[i];
      DM->St0[i]     = _DM->St0[i];
      DM->K[i]       = _DM->K[i];
      DM->T[i]       = _DM->T[i];
      DM->r[i]       = _DM->r[i];
      DM->d[i]       = _DM->d[i];
      DM->wi[i]      = _DM->wi[i];
      DM->PrixObs[i] = _DM->PrixObs[i];
      DM->SigmaImpObs[i] = _DM->SigmaImpObs[i];
    }
  //
  TypeNorme     = _TypeNorme;
  TypeModel     = _TypeModel;
  LogNorme      = _LogNorme;
  MelangeNormes = _MelangeNormes;
  //
  initOk        = 1;
  //
}
//
double dNormeMarket(int TypeNorme,int TypeOpt,double St0,double T,double K,double r,double divid,double prix,double prixobs,double sigma_imp_obs,int dimx,double *gradi)
{
  //
  int j1,i,ibidon;
  double fx,sigma_imp,dprice,sigma_imp0,dsigma_imp;
  //
  // Cas 1 : (ecart des prix)^2
  if(TypeNorme == 1)
    {
      fx     = (prix - prixobs);
	  for(i=0;i<dimx;i++)
	    gradi[i] = gradi[i] * 2.* fx;
	  fx     = fx * fx;
    }
  // Cas 2 : (ecart des prix normalise)^2
  else if(TypeNorme == 2)
    {
      fx     = (prix - prixobs)/prixobs;
      for(i=0;i<dimx;i++)
	gradi[i] = gradi[i] * 2.* fx/prixobs;
      fx     = fx * fx;
    }
  // Cas 3 : (ecart des SigmaImplicites)^2
  else if(TypeNorme == 3)
    {
      sigma_imp     = SigmaImplicite(eps,a0,b0,TypeOpt,St0,T,K,r,divid,prix,&j1);
      //
      if(j1==-1) {
	sigmaimpko++;
	if(MelangeNormes==1){  
	  // Si pas de solutions trouvees dans SigmaImplicte. On utilise les
	  // ecarts de prix normalises.
	  fx = dNormeMarket(2,TypeOpt,St0,T,K,r,divid,prix,prixobs,sigma_imp_obs,dimx,gradi);
	} else {
	  fx = 0.;
	  for(i=0;i<dimx;i++)
	    gradi[i] = 0.;
	}
	//
      } else {
	sigmaimpok++;
	if(TypeOpt==1) {
	  ibidon =  dCall_BlackScholes_73(St0,K,T,r,divid,sigma_imp,&dprice);
	} else {
	  ibidon =  dPut_BlackScholes_73(St0,K,T,r,divid,sigma_imp,&dprice);
	}
	//
	//
	if(dprice==0.) {
	  // si dprice == 0 ==> dsigma ~ \infinie 
	  // printf("sigma_imp=%f,*e8=%f\n",sigma_imp,sigma_imp*1.e8);
	  //printf("TypeOpt=%d,St0=%f,K=%f,T=%f,r=%f,divid=%f\n",TypeOpt,St0,K,T,r,divid);
	  dprice = 1.e-3;
	  // 
	}
	// pour empecher dprice d'exploser.
	if(dprice<-1.e8) { dprice=-1.e8;}
	else if(dprice>1.e8) { dprice=1.e8;}		
	//
	dsigma_imp = 1. / dprice ;
	//
	sigma_imp0 = 1.e-2;
	//
	fx     = (sigma_imp - sigma_imp_obs)/sigma_imp0;
	for(i=0;i<dimx;i++)
	  gradi[i] = gradi[i] * 2.* fx * dsigma_imp/sigma_imp0;
	fx     = fx * fx;
	//
      }
      //
    } else {
      printf("Erreur dans NormeMarket, TypeNorme = %d inprevu\n",TypeNorme);  
      exit(-1);
    }
  //
  return fx;
}
//
void FreeGradFunction()
{
  FreeDataMarket(DM);
  DM = NULL;
}
//
