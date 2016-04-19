#ifndef _SVJ_H_
#define _SVJ_H_

#include "optype.h"
#include "math/integral.h"

typedef struct SVJPARAMS
{
  
  /* phi = 1 pour un call et -1 pour un put*/
  double phi;
  
  /* si type_f = 1 alors utiliser la formule type BS avec les 2 integrales
     et type_f = 2 alors utiliser la formule avec une seule integrale.*/
  int type_f;
 
  /* les parametres de l'option*/
  double St0,K,T,r,divid;
 
  /* les parametres pour la volatlite stochastique
     si heston = 0, alors volatilite constante = sqrt(V0)
     sinon modele de Heston.*/
  int heston;
  double kappa,theta,sigmav,rho,V0;
  
  /* les parametres du saut
     si merton = 0, alors pas de saut, sinon sauts d'intensite lambda
     et de type log-normal.*/
  int merton;
  double lambda,m0,v;
  
}
  SVJPARAMS;

int calc_price_svj(SVJPARAMS *svj,double *ptprice, double *ptdelta);

#endif
