#include "costFunction.h"
// Variables globales pour ce fichier
static int initOk=0,TypeNorme,TypeModel;
static int LogNorme = 0;
static int MelangeNormes = 0;
static int sigmaimpok=0,sigmaimpko=0;
static DataMarket *DM;
//
double costFunction(int dimx, double *x)
{
  //
  int i;
  double fx=0,fxi;
  int ibidon;
  double V0,kappa,theta,sigmav,rho,lambda,m0,v,prix,delta;
  double St0,K,T,prixobs,sigma_imp_obs;
  int nbok,nbko;
  NumFunc_1 *p;
  double ri,di;
  //
  p=(NumFunc_1*)malloc(sizeof(NumFunc_1));
  //
  nbok = 0;
  nbko = 0;
  // 
  switch(TypeModel)
	{
	  // La fonction x^2 pour tester
	case 0:
	  //
	  for(i=0;i<dimx;i++)
	    fx = fx + x[i]*x[i];
	  //
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
	      ri = DM->r[i];
	      di = DM->d[i];
	      //
	      if(DM->TypeOpt[i] == 1) {
		ibidon=FT_Call_Heston(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,&prix,&delta);
	      } else {
		ibidon=FT_Put_Heston(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,&prix,&delta);
	      }
	      //
	      prixobs = DM->PrixObs[i];
	      sigma_imp_obs = DM->SigmaImpObs[i];
	      fxi = NormeMarket(TypeNorme,DM->TypeOpt[i],St0,T,K,ri,di,prix,prixobs,sigma_imp_obs,&DM->SigmaImp[i]);
	      //
	      DM->PrixMod[i] = prix;
	      //
	      fx = fx + fxi*DM->wi[i];
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
		  } else {
			ibidon=FT_Put_Merton(St0,p,T,ri,di,V0,lambda,m0,v,&prix,&delta);
		  }
		  //
		  prixobs = DM->PrixObs[i];
		  sigma_imp_obs = DM->SigmaImpObs[i];
		  fxi = NormeMarket(TypeNorme,DM->TypeOpt[i],St0,T,K,ri,di,prix,prixobs,sigma_imp_obs,&DM->SigmaImp[i]);
		  //
		  fx = fx + fxi*DM->wi[i];
		  //
		  DM->PrixMod[i] = prix;
		  //
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
		  } else {
			ibidon=FT_Put_HestMert(St0,p,T,ri,di,V0,kappa,theta,sigmav,rho,lambda,m0,v,&prix,&delta);
		  }
		  //
		  prixobs = DM->PrixObs[i];
		  sigma_imp_obs = DM->SigmaImpObs[i];
		  fxi = NormeMarket(TypeNorme,DM->TypeOpt[i],St0,T,K,ri,di,prix,prixobs,sigma_imp_obs,&DM->SigmaImp[i]);
		  //
		  fx = fx + fxi*DM->wi[i];
		  //
		  DM->PrixMod[i] = prix;		  
		  //
		 }
	  break;
	  //
	}
  //
  free(p);
  //
  // si on prend la norme : ln(fx) pour ameliorer la recherche lorsque on est proche de 
  // la solution. 
 if(LogNorme==1)
	{
	  if(fx>0.) fx = log(fx);
	}
  //
  return fx;
}
//
void initParams(int nbdata,DataMarket *_DM,int _TypeNorme,int _TypeModel,int _LogNorme,int _MelangeNormes)
{
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
double NormeMarket(int TypeNorme,int TypeOpt,double St0,double T,double K,double r,double divid,double prix,double prixobs,double sigma_imp_obs,double *sigmaimp)
{
  //
  int j1;
  double fx,sigma_imp,sigma_imp0;
  //
  // Cas 1 : (ecart des prix)^2
  if(TypeNorme == 1)
	{
	  fx     = (prix - prixobs);
	  fx     = fx * fx;
	}
  // Cas 2 : (ecart des prix normalisees)^2
  else if(TypeNorme == 2)
	{
	  fx     = (prix - prixobs)/prixobs;
	  fx     = fx * fx;
	}
  // Cas 3 : (ecart des SigmaImplicites)^2
  else if(TypeNorme == 3)
	{
	  sigma_imp     = SigmaImplicite(eps,a0,b0,TypeOpt,St0,T,K,r,divid,prix,&j1);
	  *sigmaimp = sigma_imp;
	  //	  
	  if(j1==-1) {
		sigmaimpko++;
		if(MelangeNormes==1) {
		  // Si pas de solutions trouvees dans SigmaImplicte. On utilise les
		  // ecarts de prix normalises.
		  fx = NormeMarket(2,TypeOpt,St0,T,K,r,divid,prix,prixobs,sigma_imp_obs,&sigma_imp);
		} else {
		  //sinon on ne prends pas ce point
		  fx = 0.;
		}
		//
	  } else {
		sigmaimpok++;
		sigma_imp0 = 1.e-2;
		
		fx     = (sigma_imp - sigma_imp_obs)/sigma_imp0;
		fx     = fx * fx;
		// si on veut faire du x^4
		//fx     = fx*fx;
	  }
	  //
	}
    else {
	  printf("Erreur dans NormeMarket, TypeNorme = %d inprevu\n",TypeNorme);  
	  exit(-1);
	}
  //  
  return fx;
}
//
void checksigmaimp()
{
  printf("SigmImp OK = %d et Sigma Imp KO = %d \n",sigmaimpok,sigmaimpko);
}
//
void FreeCostFunction()
{
  FreeDataMarket(DM);
  DM = NULL;  
}
//
void AffichePrixDataMarket()
{
  int i;
  double p1,p2,e1,e2,sig1,sig2;
  printf("DataMarket : %d  donnees\n",DM->nbdata);
  printf("    TypeOpt,      St0,             T,             K,             PrixObs,          PrixMod,           Erreur (%%)      SigmaImpObs,       SigmaImp       Erreur\n");
  for(i=0;i<DM->nbdata;i++)
    {
      p1 = DM->PrixObs[i];
      p2 = DM->PrixMod[i];
      sig1 = DM->SigmaImp[i];
      sig2 = DM->SigmaImpObs[i];
      e1 = 100.*fabs(p1-p2)/p1;
      e2 = 100.*(sig1 - sig2);
      printf("        %d     %f       %f       %f         %f         %f            %f         %f       %f        %f\n",DM->TypeOpt[i],DM->St0[i],DM->T[i],DM->K[i],p1,p2,e1,sig1,sig2,e2); 
    }  
}
//



