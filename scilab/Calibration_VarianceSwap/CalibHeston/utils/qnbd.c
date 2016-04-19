#include "optim.h"
static int optim_zqnbd (int *indqn, opt_simul simul, double *dh, int *n, double *binf,
			double *bsup, double *x, double *f, double *g, double *zero,
			int *napmax, int *itmax, int *indic, int *izig, int *nfac,
			int *imp, int *io, double *epsx, double *epsf, double *epsg,
			double *x1, double *x2, double *g1, double *dir, double *df0,
			int *ig, int *in, int *irel, int *izag, int *iact,
			double *epsrel, int *ieps1, opt_simul_data *optim_data);

static int optim_ajour (int *mode, int *n, int *nc, int *nr, double *h__, double *w, int *indi);
static int optim_calmaj (double *dh, int *n, double *g1, double *sig, double *w, int *ir,
			 int *mk, double *epsmc, int *nfac);
/*!but 
 *    code de minimisation d une fonction reguliere sous contraintes 
 *    de bornes , aux normes modulopt 
 *!origine 
 *    f. bonnans , inria , 1986 
 *!methode 
 *    principe de l algorithme : quasi-newton + projection 
 *    details dans le rapport inria n. 242,1983 
 *    cette version permet de tester plusieurs variantes de l algorithme 
 *    en modifiant certains parametres internes (cf comment dans le code) 
 *    taille memoire necessaire de l ordre de n**2/2 
 *    pour les problemes de grande taille le code gcbd est mieux adapte 
 * 
 *!sous programmes appeles 
 *    zqnbd   optimiseur effectif 
 *    proj    projection 
 *    calmaj  mise a jour du hessien 
 *    ajour mise a jour des facteurs de choleski 
 *    rlbd,satur   recherche lineaire avec bornes 
 * 
 *!liste d'appel 
 * 
 *    subroutine qnbd(indqn,simul,n,x,f,g,imp,io,zero, 
 *   & napmax,itmax,epsf,epsg,epsx,df0,binf,bsup,nfac, 
 *   & trav,ntrav,itrav,nitrav,izs,rzs,dzs) 
 * 
 *    indqn   indicateur de qnbd                                  es 
 *      en entree =1 standard 
 *                =2 dh et indic initialises au debut de trav et itrav 
 *                   ifac,f,g initialises 
 *      en sortie 
 *       si < 0 incapacite de calculer un point meilleur que le point initial 
 *       si = 0 arret demande par l utilisateur 
 *       si > 0 on fournit un point meilleur que le point de depart 
 *       < -10 parametres d entree non convenables 
 *       = -6  arret lors du calcul de la direction de descente et iter=1 
 *       = -5  arret lors du calcul de l approximation du hessien  iter=1 
 *       = -3  anomalie de simul : indic negatif en un point ou 
 *             f et g ont ete precedemment calcules 
 *       = -2  echec de la recherche lineaire a la premiere iteration 
 *       = -1  f non definie au point initial 
 *       =  1  arret sur epsg 
 *       =  2            epsf 
 *       =  3            epsx 
 *       =  4            napmax 
 *       =  5            itmax 
 *       =  6  pente dans la direction opposee au gradient trop petite 
 *       =  7  arret lors du calcul de la direction de descente 
 *       =  8  arret lors du calcul de l approximation du hessien 
 *       = 10  arret par echec de la recherche lineaire , cause non precisee 
 *       = 11  idem avec indsim < 0 
 *       = 12            un pas trop petit proche d un pas trop grand 
 *                       ceci peut resulter d une erreur dans le gradient 
 *       = 13            trop grand nombre d appels dans une recherche lineaire 
 *    simul voir les normes modulopt 
 *    n dim de x                                                 e 
 *    binf,bsup bornes inf,sup,de dim n                          e 
 *    x variables a optimiser (controle)                          es 
 *    f valeur du critere                                         s 
 *    g gradient de f                                             s 
 *    zero  proche zero machine                                             e 
 *    napmax nombre maximum d appels de simul                               e 
 *    itmax nombre maximum d iterations de descente               e 
 *    itrav vect travail dim nitrav=2n , se decompose en indic et izig 
 *    nfac nombre de variables factorisees                  (e si indqn=2)  s 
 *    imp facteur d impression                                              e 
 *        varie de 0 (pas d impressions) a 3 (nombreuses impressions) 
 *    io numero du fichier de resultats                                     e 
 *    epsx vect dim n precision sur x                                       e 
 *    epsf critere arret sur f                                              e 
 *    epsg arret si sup a norm2(g+)/n                                       e 
 *    trav vect travail dim ntrav 
 *    il faut ntrav > n(n+1)/2 +6n 
 *    df0>0 decroissance f prevue     (prendre 1. par defaut)               e 
 *    izs,rzs,dzs : cf normes modulopt                                     es 
 *! 
 *    indications sur les variables internes a qnbd et zqnbd 
 *    izig  sert a la memorisation des contraintes (actif si izag>1) 
 *    si i ne change pas d ens on enleve 1 a izig  (positif) 
 *    sinon on ajoute izag 
 *    factorisation seulement si izig est nul 
 *    dh estimation hessien dim n(n+1)/2 rangee en troismorceaux es 
 *    indic(i) nouvel indice de l indice i 
 *    indic vect dim n ordre de rangement des indices                       es 
 *    pas necessaire de l initialiser si indqn=1 
 * 
 *    parametres de la recherche lineaire 
 *    amd,amf param. du test de wolfe .    (.7,.1) 
 *    napm nombre max d appels dans la rl  (=15) 
 * 
 * 
 */



int optim_qnbd (int *indqn, opt_simul simul, int *n, double *x, double *f, double *g,
		int *imp, int *io, double *zero, int *napmax, int *itmax,
		double *epsf, double *epsg, double *epsx, double *df0,
		double *binf, double *bsup, int *nfac, double *trav, int *ntrav,
		int *itrav, int *nitrav, opt_simul_data *optim_data)
{
  /* Format strings */
  /* Local variables */
  int iact, izag, irel, ieps1;
  int n1, n2, n3, n4, n5, ig, in;
  double epsrel;
  int ni1, ni2;

  /* Parameter adjustments */
  --bsup;
  --binf;
  --epsx;
  --g;
  --x;
  --trav;
  --itrav;

  /* Function Body */
  if (*imp >= 1)
    {
      Sciprintf( "qnbd: start\n");
    }
  /* 
   * 
   *    parametres caracteristiques de l algorithme 
   *    si les parametres sont nuls l algorithme est celui du rr 242 
   *    ig=1 test sur grad(i) pour relach var 
   *    in=1 limite le nombre de factorisations par iter a n/10 
   *    irel=1 test sur decroissance grad pour rel a iter courante 
   *    epsrel taux de decroissance permettant le relachement (cf irit) 
   *    iact blocage variables dans ib (gestion contraintes actives) 
   *    ieps1 =1 eps1 egal a zero 
   *    note eps1 correspond a eps(xk) 
   */
  ig = 0;
  in = 0;
  irel = 1;
  epsrel = .5;
  izag = 0;
  iact = 1;
  ieps1 = 0;
  /* 
   *    decoupage du vecteur trav 
   */
  n1 = *n * (*n + 1) / 2 + 1;
  n2 = n1 + *n;
  n3 = n2 + *n;
  n4 = n3 + *n;
  n5 = n4 + *n - 1;
  if (*ntrav < n5)
    {
      if (*imp > 0)
	{
	  Sciprintf("qnbd: ntrav=%d, devrait valoir %d\n",
		    (*ntrav),
		    n5);
	}
      *indqn = -11;
      return 0;
    }
  ni1 = *n + 1;
  if (*nitrav < *n << 1)
    {
      ni2 = *n << 1;
      if (*imp > 0)
	{
	  Sciprintf("qnbd: nitrav=%d, devrait valoir %d\n",
		    (*nitrav),
		    ni2);
	}
      *indqn = -12;
      return 0;
    }
  optim_zqnbd (indqn,  simul, &trav[1], n, &binf[1], &bsup[1], &x[1], f,
	       &g[1], zero, napmax, itmax, &itrav[1], &itrav[ni1], nfac, imp,
	       io, &epsx[1], epsf, epsg, &trav[n1], &trav[n2], &trav[n3],
	       &trav[n4], df0, &ig, &in, &irel, &izag, &iact, &epsrel, &ieps1,
	       optim_data);
  return 0;
}





static int optim_zqnbd (int *indqn, opt_simul simul, double *dh, int *n, double *binf,
			double *bsup, double *x, double *f, double *g, double *zero,
			int *napmax, int *itmax, int *indic, int *izig, int *nfac,
			int *imp, int *io, double *epsx, double *epsf, double *epsg,
			double *x1, double *x2, double *g1, double *dir, double *df0,
			int *ig, int *in, int *irel, int *izag, int *iact,
			double *epsrel, int *ieps1, opt_simul_data *optim_data)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  int ifac;
  double diff, difg, scal;
  int mode, napm;
  double teta;
  int iter, irit;
  double tmax;
  int nfac1;
  double difg0, difg1;
  int n2fac;
  double scal1;
  int napm1;
  double teta1, zsig1;
  int idfac, i__, j, k;
  double t;
  int nnfac;
  double v, y, epsmc;
  int indrl, iconv;
  double d1, tiers, d2;
  int i1;
  double tproj;
  int n1, n3;
  double t1, cscal1, aa, dd, bi;
  int ic, ii, ij;
  double fn, bs, ep;
  int ip, mk, ir;
  double gr;
  int np, indsim, nm1;
  double amd, amf;
  int ndh, nap, ifp;
  double sig, fpn;
  int nip;
  double cof1, cof2, sig1, eps0, eps1;

  /*    Copyright INRIA 
   */

  /* Parameter adjustments */
  --dh;
  --dir;
  --g1;
  --x2;
  --x1;
  --epsx;
  --izig;
  --indic;
  --g;
  --x;
  --bsup;
  --binf;

  /* Function Body */
  if (*imp < 4)
    {
      goto L3;
    }
  Sciprintf("qnbd: izag=%d,ig=%d,in=%d,irel=%d,iact=%d,epsrel=%11.4f\n",
	    *izag,
	    *ig,
	    *in,
	    *irel,
	    *iact,
	    *epsrel);
  /* 
   */
  if (*ig == 1)
    {
      Sciprintf("test sur gradient pour sortie ib\n");
    }
  if (*in == 1)
    {
      Sciprintf("test sur nombre de defactorisations pour sortie ib\n");
    }
  if (*izag != 0)
    { 
      Sciprintf("memorisation de variables izag=,%d\n", *izag);
    }
  if (*irel == 1)
    {  
      Sciprintf("methode de minimisations incompletes ; epsrel=,%11.4f\n", *epsrel);
    }
  if (*iact == 1)
    {
      Sciprintf("blocage des variables dans ib\n" );
    }
  if (*ieps1 == 1)
    {  Sciprintf("parametre eps1 nul\n" );
    }
  if (*ieps1 == 2)
    {   Sciprintf("parametre eps1 nul\n");
    }
  /* 
   *    cscal1 utilise pour calculer eps(x) = eps1 cf avant 310 
   */
  cscal1 = 1e8;
  if (*ieps1 == 2)
    { 
      Sciprintf("parametre eps1=eps(x) calcule avec cscal1=,%11.4f\n",cscal1);
    }
 L3:
  /* 
   */
  difg0 = 1.;
  difg1 = 0.;
  /* 
   *    eps0 sert a partitionner les variables 
   */
  eps0 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      izig[i__] = 0;
      /* L5: */
      eps0 += epsx[i__];
    }
  eps0 = eps0 * 10. / *n;
  /* 
   *    section 1  mise en forme de dh 
   *    si indqn=1 on init dh a ident puis scal a it=2 
   * 
   */
  optim_proj (n, &binf[1], &bsup[1], &x[1]);
  ndh = *n * (*n + 1) / 2;
  if (*indqn == 1)
    {
      goto L10;
    }
  if (*indqn == 2)
    {
      goto L30;
    }
  /*    erreur 
   */
  if (*imp > 0)
    {
      Sciprintf("qnbd: valeur non admissible de indqn %d\n", *indqn);
    }
  *indqn = -105;
  if (*imp > 0)
    {
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
 L10:
  /*    on initialise dh a l identite puis a l iteration 2 
   *    on met a l echelle 
   */
  *nfac = 0;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L11: */
      indic[i__] = i__;
    }
  i__1 = ndh;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L12: */
      dh[i__] = 0.;
    }
 L30:
  /* 
   *    section 2  mise a jour dh 
   * 
   *    iter nombre d iterations de descente 
   */
  iter = 0;
  scal = 1.;
  nap = 1;
  indsim = 4;
  if (*indqn == 1)
    {
      (*simul) (&indsim, n, &x[1], f, &g[1], optim_data);
    }
  if (indsim <= 0)
    {
      *indqn = -1;
      if (indsim == 0)
	{
	  *indqn = 0;
	}
      if (*imp > 0)
	{
	  Sciprintf("qnbd: indqn=%d\n",*indqn);
	}
      return 0;
    }
  if (*indqn != 1)
    {
      goto L200;
    }
  /*    mise a echelle dh 
   *    df0 decroissance prevue . si mod quad df0=((dh)-1g,g)/2 
   *    et on cherche dh diag de la forme cst/(dx)**2 
   *    d ou cst=som((y(i)*(dx))**2))/(2*df0) 
   */
  cof1 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L80: */
      /*Computing 2nd power 
       */
      d__1 = g[i__] * epsx[i__];
      cof1 += d__1 * d__1;
    }
  cof1 /= *df0 * 2.;
  i1 = -(*n);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i1 = i1 + *n + 2 - i__;
      /* L82: */
      /*Computing 2nd power 
       */
      d__1 = epsx[i__];
      dh[i1] = (cof1 + *zero) / (d__1 * d__1 + *zero);
    }
  iconv = 0;
 L200:
  ++iter;
  if (iter <= *itmax)
    {
      goto L202;
    }
  if (*imp > 0)
    { 
      Sciprintf("qnbd: maximum d iterations atteint\n");
    }
  *indqn = 5;
  if (*imp > 0)
    {
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
 L202:
  if (*imp >= 2)
    { 
      Sciprintf("(/ qnbd: iter=,%d,  f=,%15.7f\n", iter, *f );
    }
  /*    x1,g1 valeurs a l iteration precedente 
   */
  if (iter == 1)
    {
      goto L300;
    }
  cof1 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      x1[i__] = x[i__] - x1[i__];
      g1[i__] = g[i__] - g1[i__];
      /* L201: */
      cof1 += x1[i__] * g1[i__];
    }
  if (cof1 <= *zero)
    {
      goto L250;
    }
  if (iter > 2 || *indqn != 1)
    {
      goto L250;
    }
  /*    mise a l echelle de dh par methode shanno-phua 
   *     dh=(y,y)/(y,s)*id 
   */
  cof2 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L203: */
      /*Computing 2nd power 
       */
      d__1 = g1[i__];
      cof2 += d__1 * d__1;
    }
  cof2 /= cof1;
  if (*imp > 3)
    { 
      Sciprintf("qnbd: facteur d echelle=,%11.4f\n",  cof2);
    }
  dh[1] = cof2;
  i1 = 1;
  i__1 = *nfac;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i1 = i1 + *n + 1 - i__;
      /* L205: */
      dh[i1] = cof2;
    }
  /* 
   *    scal= (y,s)/(y,y) 
   *    scal sert de coeff a g dans le calcul de dir pour i dans ib 
   */
  scal = 1. / cof2;
 L250:
  /* 
   *    mise a jour dh par methode bfgs (majour) si iter ge 2 
   *    dh1=dh +y*yt/(y,s) - dh*s*st*dh/(s,dh*s) 
   *    exprimons ds=x1 et y=g1 dans les nouv variables soit x2 et g1 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i1 = indic[i__];
      x2[i1] = g1[i__];
      /* L251: */
      dir[i1] = x1[i__];
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L252: */
      g1[i__] = x2[i__];
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i1 = indic[i__];
      /* L253: */
      x2[i1] = x1[i__];
    }
  /*    on stocke d abord dh*s dans x2 
   *    calcul des nfac premieres variables,en deux fois 
   */
  if (*nfac == 0)
    {
      goto L2312;
    }
  if (*nfac > 1)
    {
      goto L2300;
    }
  dir[1] *= dh[1];
  goto L2312;
 L2300:
  np = *nfac + 1;
  ii = 1;
  n1 = *nfac - 1;
  i__1 = n1;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      y = dir[i__];
      if (dh[ii] == 0.)
	{
	  goto L2302;
	}
      ij = ii;
      ip = i__ + 1;
      i__2 = *nfac;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ij;
	  /* L2301: */
	  y += dir[j] * dh[ij];
	}
    L2302:
      dir[i__] = y * dh[ii];
      /* L2303: */
      ii = ii + np - i__;
    }
  dir[*nfac] *= dh[ii];
  i__1 = n1;
  for (k = 1; k <= i__1; ++k)
    {
      i__ = *nfac - k;
      ii = ii - np + i__;
      if (dir[i__] == 0.)
	{
	  goto L2311;
	}
      ip = i__ + 1;
      ij = ii;
      y = dir[i__];
      i__2 = *nfac;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ij;
	  /* L2310: */
	  dir[j] += dh[ij] * dir[i__];
	}
    L2311:
      ;
    }
 L2312:
  nfac1 = *nfac + 1;
  n2fac = *nfac * nfac1 / 2;
  nnfac = *n - *nfac;
  k = n2fac;
  if (*nfac == *n)
    {
      goto L268;
    }
  i__1 = *n;
  for (i__ = nfac1; i__ <= i__1; ++i__)
    {
      /* L255: */
      dir[i__] = 0.;
    }
  if (*nfac == 0)
    {
      goto L265;
    }
  i__1 = *nfac;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i__2 = *n;
      for (j = nfac1; j <= i__2; ++j)
	{
	  ++k;
	  if (x2[j] == 0.)
	    {
	      goto L260;
	    }
	  dir[i__] += dh[k] * x2[j];
	L260:
	  ;
	}
    }
  /*    calcul autres comp de dh*s=d en deux fois 
   */
  k = n2fac;
  i__2 = *nfac;
  for (j = 1; j <= i__2; ++j)
    {
      i__1 = *n;
      for (i__ = nfac1; i__ <= i__1; ++i__)
	{
	  ++k;
	  dir[i__] += dh[k] * x2[j];
	  /* L264: */
	}
    }
 L265:
  k = n2fac + *nfac * nnfac;
  i__1 = *n;
  for (j = nfac1; j <= i__1; ++j)
    {
      i__2 = *n;
      for (i__ = j; i__ <= i__2; ++i__)
	{
	  ++k;
	  if (x2[j] == 0.)
	    {
	      goto L266;
	    }
	  dir[i__] += dh[k] * x2[j];
	L266:
	  ;
	}
    }
  if (*nfac == *n - 1)
    {
      goto L268;
    }
  nm1 = *n - 1;
  k = n2fac + *nfac * nnfac;
  i__2 = nm1;
  for (i__ = nfac1; i__ <= i__2; ++i__)
    {
      ++k;
      i1 = i__ + 1;
      i__1 = *n;
      for (j = i1; j <= i__1; ++j)
	{
	  ++k;
	  if (x2[j] == 0.)
	    {
	      goto L267;
	    }
	  dir[i__] += dh[k] * x2[j];
	L267:
	  ;
	}
    }
  /*    calcul de dh*s fini 
   *    calcul sig1 pour 2eme mise a jour 
   */
 L268:
  sig1 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L271: */
      sig1 += dir[i__] * x2[i__];
    }
  if (sig1 > 0.)
    {
      goto L272;
    }
  if (*imp > 2)
    {
      Sciprintf("qnbd: pb (bs,s) negatif=,%11.4f\n", sig1);
    }
  /* 
   *    ****************************************************** 
   */
  *indqn = 8;
  if (iter == 1)
    {
      *indqn = -5;
    }
  if (*imp > 0)
    {
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
 L272:
  sig1 = -1. / sig1;
  /*    truc powell si (y,s) negatif 
   */
  if (cof1 > *zero)
    {
      goto L277;
    }
  if (*imp > 2)
    { 
      Sciprintf("qnbd: emploi truc powell (y,s)=,%11.4f\n",  cof1);
    }
  teta = -1. / sig1;
  teta = teta * .8 / (teta - cof1);
  teta1 = 1. - teta;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L274: */
      g1[i__] = teta * g1[i__] + teta1 * dir[i__];
    }
  cof1 = -.2 / sig1;
 L277:
  /* 
   *     premiere mise a jour de dh 
   */
  sig = 1. / cof1;
  zsig1 = 1. / sig1;
  mk = 0;
  ir = *nfac;
  epsmc = 1e-9;
  optim_calmaj (&dh[1], n, &g1[1], &sig, &x2[1], &ir, &mk, &epsmc, nfac);
  if (ir != *nfac)
    {
      goto L280;
    }
  optim_calmaj (&dh[1], n, &dir[1], &sig1, &x2[1], &ir, &mk, &epsmc, nfac);
  if (ir != *nfac)
    {
      goto L280;
    }
  goto L300;
 L280:
  if (*imp > 0)
    { 
      Sciprintf("qnbd: pb dans appel majour\n");
    }
  *indqn = 8;
  if (iter == 1)
    {
      *indqn = -5;
    }
  if (*imp > 0)
    {
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
 L300:
  /* 
   *    section 3 determination des variables libres et bloquees 
   * 
   *    calcul eps1 
   * 
   */
  scal1 = scal;
  if (*ieps1 == 1)
    {
      scal1 = 0.;
    }
  if (*ieps1 == 2)
    {
      scal1 = scal * cscal1;
    }
  /* L305: */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L310: */
      x1[i__] = x[i__] - scal1 * (d__1 = g[i__], Abs (d__1)) * g[i__];
    }
  optim_proj (n, &binf[1], &bsup[1], &x1[1]);
  eps1 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L320: */
      eps1 += (d__1 = x1[i__] - x[i__], Abs (d__1));
    }
  eps1 = Min (eps0, eps1);
  if (*ieps1 == 1)
    {
      eps1 = 0.;
    }
  if (*ieps1 == 2)
    {
      eps1 *= 1e4;
    }
  if (*imp > 3)
    {
      Sciprintf("qnbd: val de eps1 servant a partitionner les variables,%11.4f\n", eps1);
    }
  /*    nfac nombre de lignes factorisees (nr pour ajour) 
   */
  ifac = 0;
  idfac = 0;
  k = 0;
  /* 
   * 
   */
  gr = 0.;
  if (*ig == 1)
    {
      gr = difg * .2 / *n;
    }
  n3 = *n;
  if (*in == 1)
    {
      n3 = *n / 10;
    }
  /*    si irit=1 on peut relacher des variables 
   */
  irit = 0;
  if (difg1 <= *epsrel * difg0)
    {
      irit = 1;
    }
  if (*irel == 0 || iter == 1)
    {
      irit = 1;
    }
  if (irit * *irel > 0 && *imp > 3)
    { 
      Sciprintf("qnbd: redemarrage ; difg0=%11.4f,epsrel=%11.4f,difg1=%11.4f\n",
		difg0,
		*epsrel,
		difg1);
    }
  /* 
   */
  tiers = .33333333333333331;
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      --izig[k];
      if (izig[k] <= 0)
	{
	  izig[k] = 0;
	}
      bi = binf[k];
      bs = bsup[k];
      ic = indic[k];
      d1 = x[k] - bi;
      d2 = bs - x[k];
      dd = (bs - bi) * tiers;
      ep = Min (eps1, dd);
      if (d1 > ep)
	{
	  goto L324;
	}
      if (g[k] > 0.)
	{
	  goto L330;
	}
      goto L335;
    L324:
      if (d2 > ep)
	{
	  goto L335;
	}
      if (g[k] > 0.)
	{
	  goto L335;
	}
      goto L330;
      /*    on defactorise si necessaire 
       */
    L330:
      if (ic > *nfac)
	{
	  goto L340;
	}
      ++idfac;
      mode = -1;
      if (*imp >= 4)
	{ 
	  Sciprintf("defactorisation de ,%d\n",   k);
	}
      izig[k] += *izag;
      optim_ajour (&mode, n, &k, nfac, &dh[1], &x2[1], &indic[1]);
      if (mode == 0)
	{
	  goto L340;
	}
      if (*imp > 0)
	{ 
	  Sciprintf("qnbd: pb dans ajour. mode=,%d\n",  mode);
	}
      *indqn = 8;
      if (iter == 1)
	{
	  *indqn = -5;
	}
      if (*imp > 0)
	{
	  Sciprintf("qnbd: indqn=%d\n",*indqn);
	}
      return 0;
      /*    on factorise 
       */
    L335:
      if (irit == 0)
	{
	  goto L340;
	}
      if (ic <= *nfac)
	{
	  goto L340;
	}
      if (izig[k] >= 1)
	{
	  goto L340;
	}
      mode = 1;
      if (ifac >= n3 && iter > 1)
	{
	  goto L340;
	}
      if ((d__1 = g[k], Abs (d__1)) <= gr)
	{
	  goto L340;
	}
      ++ifac;
      if (*imp >= 4)
	{ 
	  Sciprintf("on factorise l indice ,%d\n",    k);

	}
      optim_ajour (&mode, n, &k, nfac, &dh[1], &x2[1], &indic[1]);
      if (mode == 0)
	{
	  goto L340;
	}
      if (*imp > 0)
	{
	  Sciprintf("qnbd: pb dans ajour. mode=,%d\n",    mode);
	}
      *indqn = 8;
      if (iter == 1)
	{
	  *indqn = -5;
	}
      if (*imp > 0)
	{
	  Sciprintf("qnbd: indqn=%d\n",*indqn);
	}
      return 0;
    L340:
      ;
    }
  if (*imp >= 2)
    {
      Sciprintf("qnbd: nbre fact,%d, defact,%d, total var factorisees,%d\n",
		ifac,
		idfac,
		*nfac);
    }
  /* 
   *    *********************************************** a voir 
   */
  if (iconv == 1)
    {
      return 0;
    }
  /* 
   *    section 6 resolution systeme lineaire et expression de dir 
   *    on inverse le syst correspondant aux nl premieres composantes 
   *    dans le nouveau syst d indices 
   */
  ir = *nfac;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i1 = indic[i__];
      /* L640: */
      x2[i1] = g[i__];
    }
  /* L641: */
  if (ir < *nfac)
    {
      goto L412;
    }
  if (*nfac > 1)
    {
      goto L400;
    }
  x2[1] /= dh[1];
  goto L412;
 L400:
  i__1 = *nfac;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      ij = i__;
      i1 = i__ - 1;
      v = x2[i__];
      i__2 = i1;
      for (j = 1; j <= i__2; ++j)
	{
	  v -= dh[ij] * x2[j];
	  /* L401: */
	  ij = ij + *nfac - j;
	}
      x2[i__] = v;
      /* L402: */
      x2[i__] = v;
    }
  x2[*nfac] /= dh[ij];
  np = *nfac + 1;
  i__1 = *nfac;
  for (nip = 2; nip <= i__1; ++nip)
    {
      i__ = np - nip;
      ii = ij - nip;
      v = x2[i__] / dh[ii];
      ip = i__ + 1;
      ij = ii;
      i__2 = *nfac;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ii;
	  /* L410: */
	  v -= dh[ii] * x2[j];
	}
      /* L411: */
      x2[i__] = v;
    }
 L412:
  if (ir == *nfac)
    {
      goto L660;
    }
  if (*imp > 0)
    {
      Sciprintf("qnbd: pb num dans mult par inverse\n");
    }
  *indqn = 7;
  if (iter == 1)
    {
      *indqn = -6;
    }
  if (*imp > 0)
    {
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
 L660:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i1 = indic[i__];
      dir[i__] = -g[i__] * scal;
      /* L610: */
      if (i1 <= *nfac)
	{
	  dir[i__] = -x2[i1];
	}
    }
  /* 
   *    gestion contraintes actives (si iact=1) 
   */
  if (*iact != 1)
    {
      goto L675;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (izig[i__] > 0)
	{
	  dir[i__] = 0.;
	}
      if (indic[i__] > *nfac)
	{
	  dir[i__] = 0.;
	}
      /* L670: */
    }
 L675:
  /* 
   *    recherche lineaire 
   *    conservation de x et g . calcul de dir+ et fpn 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      g1[i__] = g[i__];
      /* L700: */
      x1[i__] = x[i__];
    }
  /*    ifp =1 si fpn trop petit. on prend alors d=-g 
   */
  ifp = 0;
  fn = *f;
 L709:
  fpn = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (x[i__] - binf[i__] <= epsx[i__] && dir[i__] < 0.)
	{
	  dir[i__] = 0.;
	}
      if (bsup[i__] - x[i__] <= epsx[i__] && dir[i__] > 0.)
	{
	  dir[i__] = 0.;
	}
      /* L710: */
      fpn += g[i__] * dir[i__];
    }
  if (fpn > 0.)
    {
      if (ifp == 1)
	{
	  if (*imp > 0)
	    { 
	      Sciprintf("qnbd: arret fpn non negatif=,%11.4f\n", fpn);
	    }
	  *indqn = 6;
	  if (iter == 1)
	    {
	      *indqn = -3;
	    }
	  if (*imp > 0)
	    {
	      Sciprintf("qnbd: indqn=%d\n",*indqn);
	    }
	  return 0;
	}
      else
	{
	  ifp = 1;
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      /* L707: */
	      if (izig[i__] > 0)
		{
		  dir[i__] = -scal * g[i__];
		}
	    }
	  irit = 1;
	  goto L709;
	}
    }
  /*    calcul du t initial suivant une idee de fletcher 
   */
  t1 = t;
  if (iter == 1)
    {
      diff = *df0;
    }
  t = diff * -2. / fpn;
  if (t > .3 && t < 3.)
    {
      t = 1.;
    }
  if (eps1 < eps0)
    {
      t = 1.;
    }
  if (*indqn == 2)
    {
      t = 1.;
    }
  if (iter > 1 && t1 > .01 && t1 < 100.)
    {
      t = 1.;
    }
  tmax = 1e10;
  t = Min (t, tmax);
  /*Computing MAX 
   */
  d__1 = t, d__2 = *zero * 10.;
  t = Max (d__1, d__2);
  /*    amd,amf tests sur h'(t) et diff 
   */
  amd = .7;
  amf = .1;
  napm = 15;
  napm1 = nap + napm;
  if (napm1 > *napmax)
    {
      napm1 = *napmax;
    }
  optim_rlbd (&indrl, n,  simul, &x[1], &binf[1], &bsup[1], &fn, &fpn,
	      &t, &tmax, &dir[1], &g[1], &tproj, &amd, &amf, imp, io, zero,
	      &nap, &napm1, &x2[1], optim_data);
  if (indrl >= 10)
    {
      indsim = 4;
      ++nap;
      (*simul) (&indsim, n, &x[1], f, &g[1], optim_data);
      if (indsim <= 0)
	{
	  *indqn = -3;
	  if (indsim == 0)
	    {
	      *indqn = 0;
	    }
	  if (*imp > 0)
	    {
	      Sciprintf("qnbd: indqn=%d\n",*indqn);
	    }
	  return 0;
	}
    }
  if (indrl <= 0)
    {
      *indqn = 10;
      if (indrl == 0)
	{
	  *indqn = 0;
	}
      if (indrl == -3)
	{
	  *indqn = 13;
	}
      if (indrl == -4)
	{
	  *indqn = 12;
	}
      if (indrl <= -1000)
	{
	  *indqn = 11;
	}
      if (*imp > 0)
	{ 
	  Sciprintf("qnbd: indqn=%d\n",*indqn);
	}
      return 0;
    }
  /* 
   */
  /* L753: */
  if (*imp < 6)
    {
      goto L778;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L760: */ 
      Sciprintf("i=%d, x=%11.4f, g=%11.4f, dir=%11.4f\n", i__,
		x[i__],
		g[i__],
		dir[i__]);
    }
  /* 
   */
 L778:
  if (nap < *napmax)
    {
      goto L758;
    }
  *f = fn;
  if (*imp > 0)
    { 
      Sciprintf("qnbd: retour cause max appels simul %d\n", *napmax);
    }
  *indqn = 4;
  if (*imp > 0)
    {
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
 L758:
  /*    section 8 test de convergence 
   * 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if ((d__1 = x[i__] - x1[i__], Abs (d__1)) > epsx[i__])
	{
	  goto L806;
	}
      /* L805: */
    }
  *f = fn;
  if (*imp > 0)
    {

      Sciprintf("qnbd: retour apres convergence de x\n");
    }
  *indqn = 3;
  if (*imp > 0)
    {
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
 L806:
  difg = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      aa = g[i__];
      if (x[i__] - binf[i__] <= epsx[i__])
	{
	  aa = Min (0., aa);
	}
      if (bsup[i__] - x[i__] <= epsx[i__])
	{
	  aa = Max (0., aa);
	}
      /* L810: */
      /*Computing 2nd power 
       */
      d__1 = aa;
      difg += d__1 * d__1;
    }
  difg1 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (indic[i__] > *nfac)
	{
	  goto L820;
	}
      aa = g[i__];
      if (x[i__] - binf[i__] <= epsx[i__])
	{
	  aa = Min (0., aa);
	}
      if (bsup[i__] - x[i__] <= epsx[i__])
	{
	  aa = Max (0., aa);
	}
      /*Computing 2nd power 
       */
      d__1 = aa;
      difg1 += d__1 * d__1;
    L820:
      ;
    }
  difg1 = sqrt (difg1);
  difg = sqrt (difg);
  difg /= sqrt ((double) (*n));
  diff = (d__1 = *f - fn, Abs (d__1));
  *df0 = -diff;
  if (irit == 1)
    {
      difg0 = difg1;
    }
  *f = fn;
  if (*imp >= 2)
    {
      Sciprintf("qnbd: epsg=%11.4f, difg=%11.4f, epsf=%11.4f,diff=%11.4f, nap=%d\n",
		*epsg,
		difg,
		*epsf,
		diff,
		nap);
    }
  if (diff < *epsf)
    {
      *indqn = 2;
      if (*imp > 0)
	{
	  Sciprintf("qnbd: retour cause decroissance f trop petite=%11.4f\n",
		    diff);
	}
      if (*imp > 0)
	{ 
	  Sciprintf("qnbd: indqn=%d\n",*indqn);
	}
      return 0;
    }
  if (difg > *epsg)
    {
      goto L200;
    }
  *indqn = 1;
  if (*imp > 0)
    {
      Sciprintf("qnbd: retour cause gradient projete petit=,%11.4f\n", difg);
    }
  if (*imp > 0)
    { 
      Sciprintf("qnbd: indqn=%d\n",*indqn);
    }
  return 0;
} 









/*    Copyright INRIA 
 * 
 *    mode = +1 factorise la ligne nc (indices de depart) 
 *         = -1 defactorise   ' 
 *    nr nbre de lignes factorisees 
 *    h mat de dim n 
 *    w,d vect de travail 
 *    indi(i) ligne ou est stockee la ligne i de depart 
 * 
 */

static int optim_ajour (int *mode, int *n, int *nc, int *nr, double *h__, double *w, int *indi)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1;

  /* Local variables */
  double a, b, c__;
  int i__, j, k;
  double u, v, h1;
  int nsaut, i1;
  double h2, ai, di;
  int ii, ij, ik, nh, nj, nk, nl, ko;
  double wi;
  int nw;
  double di1;
  int nh1, nr1, nr2, inc;
  double hij;
  int nii, nkk, nrr, inc1;

  /* Parameter adjustments */
  --indi;
  --w;
  --h__;

  /* Function Body */
  inc = indi[*nc];
  nr1 = *nr + 1;
  nr2 = *nr - 1;
  nrr = *n - *nr;
  nii = *n - inc;
  nkk = *nr - inc;
  if (*mode == -1)
    {
      goto L240;
    }
  /* 
   *     addition d'une ligne a l 
   * 
   *         stockage des elements de la colonne inc dans w 
   */
  nsaut = nii + 1;
  nh = inc * (*n + 1) - inc * (inc + 1) / 2;
  nw = *n;
  if (inc == *n)
    {
      goto L20;
    }
  i__1 = nii;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      w[nw] = h__[nh];
      --nw;
      /* L10: */
      --nh;
    }
 L20:
  w[nr1] = h__[nh];
  --nh;
  if (inc == nr1)
    {
      goto L60;
    }
  i__1 = inc - nr1;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      nl = nii + i__ - 1;
      if (nl == 0)
	{
	  goto L35;
	}
      i__2 = nl;
      for (j = 1; j <= i__2; ++j)
	{
	  h__[nh + nsaut] = h__[nh];
	  /* L30: */
	  --nh;
	}
    L35:
      w[nw] = h__[nh];
      --nw;
      --nh;
      /* L40: */
      ++nsaut;
    }
  i__1 = inc - nr1;
  for (j = 1; j <= i__1; ++j)
    {
      h__[nh + nsaut] = h__[nh];
      /* L50: */
      --nh;
    }
  /* 
   */
 L60:
  --nw;
  nsaut = 1;
  if (*nr == 0)
    {
      goto L125;
    }
  if (inc == *n)
    {
      goto L80;
    }
  i__1 = nii;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      h__[nh + nsaut] = h__[nh];
      /* L70: */
      --nh;
    }
 L80:
  if (*nr == 1)
    {
      goto L110;
    }
  i__1 = nr2;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      w[nw] = h__[nh];
      --nw;
      --nh;
      ++nsaut;
      if (*n == nr1)
	{
	  goto L100;
	}
      i__2 = *n - nr1;
      for (j = 1; j <= i__2; ++j)
	{
	  h__[nh + nsaut] = h__[nh];
	  /* L90: */
	  --nh;
	}
    L100:
      ;
    }
 L110:
  w[nw] = h__[nh];
  --nh;
  ++nsaut;
  if (inc == nr1)
    {
      goto L125;
    }
  i__1 = inc - nr1;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      h__[nh + nsaut] = h__[nh];
      /* L120: */
      --nh;
    }
  /*        mise a jour de l 
   */
 L125:
  if (*nr != 0)
    {
      goto L130;
    }
  if (w[1] > 0.)
    {
      goto L220;
    }
  *mode = -1;
  return 0;
 L130:
  if (*nr == 1)
    {
      goto L160;
    }
  i__1 = *nr;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      ij = i__;
      i1 = i__ - 1;
      v = w[i__];
      i__2 = i1;
      for (j = 1; j <= i__2; ++j)
	{
	  v -= h__[ij] * w[j];
	  /* L140: */
	  ij = ij + *nr - j;
	}
      /* L150: */
      w[i__] = v;
    }
 L160:
  ij = 1;
  v = w[nr1];
  i__1 = *nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      wi = w[i__];
      hij = h__[ij];
      /*Computing 2nd power 
       */
      d__1 = wi;
      v -= d__1 * d__1 / hij;
      w[i__] = wi / hij;
      /* L170: */
      ij = ij + nr1 - i__;
    }
  if (v > 0.)
    {
      goto L180;
    }
  *mode = -1;
  return 0;
 L180:
  w[nr1] = v;
  /*         stockage de w dans h 
   */
  nh = *nr * (*nr + 1) / 2;
  nw = nr1;
  nsaut = nw;
  h__[nh + nsaut] = w[nw];
  --nw;
  --nsaut;
  if (*nr == 1)
    {
      goto L220;
    }
  i__1 = nr2;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      h__[nh + nsaut] = w[nw];
      --nw;
      --nsaut;
      i__2 = i__;
      for (j = 1; j <= i__2; ++j)
	{
	  h__[nh + nsaut] = h__[nh];
	  /* L200: */
	  --nh;
	}
      /* L210: */
    }
 L220:
  h__[nr1] = w[1];
  if (*n == nr1)
    {
      goto L233;
    }
  nh1 = *nr * (*n + 1) - *nr * (*nr + 1) / 2 + 1;
  nw = nr1;
  i__1 = *n - nr1;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L230: */
      h__[nh1 + i__] = w[nw + i__];
    }
  /*         mise a jour de indi 
   */
 L233:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ii = indi[i__];
      if (ii <= *nr || ii >= inc)
	{
	  goto L235;
	}
      indi[i__] = ii + 1;
    L235:
      ;
    }
  ++(*nr);
  indi[*nc] = *nr;
  *mode = 0;
  return 0;
  /* 
   *     soustraction d'une ligne a l 
   * 
   *         recherche des composantes de h 
   */
 L240:
  i__1 = *nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ik = i__;
      ij = inc;
      ii = 1;
      ko = Min (ik, inc);
      v = 0.;
      if (ko == 1)
	{
	  goto L252;
	}
      i__2 = ko - 1;
      for (k = 1; k <= i__2; ++k)
	{
	  nk = nr1 - k;
	  v += h__[ij] * h__[ik] * h__[ii];
	  ij = ij + nk - 1;
	  ii += nk;
	  /* L250: */
	  ik = ik + nk - 1;
	}
    L252:
      a = 1.;
      b = 1.;
      if (ko == i__)
	{
	  goto L253;
	}
      a = h__[ik];
    L253:
      if (ko == inc)
	{
	  goto L260;
	}
      b = h__[ij];
    L260:
      w[i__] = v + a * b * h__[ii];
    }
  /*         mise a jour de l 
   */
  if (inc == *nr)
    {
      goto L315;
    }
  inc1 = inc - 1;
  nh = inc1 * nr1 - inc1 * inc / 2 + 2;
  nh1 = nh + nkk;
  di = h__[nh - 1];
  i__1 = nkk;
  for (j = 1; j <= i__1; ++j)
    {
      di1 = h__[nh1];
      ++nh1;
      a = h__[nh];
      ai = a * di;
      /*Computing 2nd power 
       */
      d__1 = a;
      c__ = d__1 * d__1 * di + di1;
      h__[nh] = c__;
      ++nh;
      if (j == nkk)
	{
	  goto L315;
	}
      i__2 = nkk - j;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  h1 = h__[nh];
	  h2 = h__[nh1];
	  u = ai * h1 + h2 * di1;
	  h__[nh] = u / c__;
	  h__[nh1] = -h1 + a * h2;
	  ++nh;
	  ++nh1;
	  /* L300: */
	}
      ++nh;
      di = di * di1 / c__;
      /* L310: */
    }
  /*         condensation de la matrice l 
   */
 L315:
  nh = inc + 1;
  nsaut = 1;
  nj = *nr - 2;
  if (inc == 1)
    {
      ++nj;
    }
  if (*nr == 1)
    {
      goto L440;
    }
  i__1 = nr2;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i__2 = nj;
      for (j = 1; j <= i__2; ++j)
	{
	  h__[nh - nsaut] = h__[nh];
	  /* L425: */
	  ++nh;
	}
      ++nsaut;
      ++nh;
      if (i__ == inc - 1)
	{
	  goto L430;
	}
      --nj;
      if (nj == 0)
	{
	  goto L440;
	}
    L430:
      ;
    }
  /*         mise a jour de la matrice h 
   */
 L440:
  nh = *nr * nr2 / 2 + 1;
  nw = 1;
  nsaut = *nr;
  if (inc == 1)
    {
      goto L470;
    }
  i__1 = inc - 1;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      h__[nh] = w[nw];
      ++nw;
      --nsaut;
      if (*n == *nr)
	{
	  goto L455;
	}
      i__2 = nrr;
      for (j = 1; j <= i__2; ++j)
	{
	  /* L450: */
	  h__[nh + j] = h__[nh + nsaut + j];
	}
    L455:
      nh = nh + nrr + 1;
      /* L460: */
    }
 L470:
  ++nw;
  if (*nr == *n)
    {
      goto L485;
    }
  i__1 = nrr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L480: */
      w[*nr + i__] = h__[nh + nsaut + i__ - 1];
    }
  nsaut += nrr;
 L485:
  if (inc == *nr)
    {
      goto L510;
    }
  i__1 = nkk;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      --nsaut;
      h__[nh] = w[nw];
      ++nw;
      if (*nr == *n)
	{
	  goto L495;
	}
      i__2 = nrr;
      for (j = 1; j <= i__2; ++j)
	{
	  /* L490: */
	  h__[nh + j] = h__[nh + nsaut + j];
	}
    L495:
      nh = nh + nrr + 1;
      /* L500: */
    }
 L510:
  h__[nh] = w[inc];
  if (*nr == *n)
    {
      goto L540;
    }
  i__1 = nrr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L520: */
      h__[nh + i__] = w[*nr + i__];
    }
  /*         mise a jour de indi 
   */
 L540:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ii = indi[i__];
      if (ii <= inc || ii > *nr)
	{
	  goto L550;
	}
      indi[i__] = ii - 1;
    L550:
      ;
    }
  indi[*nc] = *nr;
  --(*nr);
  *mode = 0;
  return 0;
}				/* ajour_ */



/*    Copyright INRIA 
 *    subroutine de qnbd 
 */

static int optim_calmaj (double *dh, int *n, double *g1, double *sig, double *w, int *ir,
			 int *mk, double *epsmc, int *nfac)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int nfac1, n2fac, i__, j, k, nnfac;

  /* Parameter adjustments */
  --w;
  --g1;
  --dh;

  /* Function Body */
  if (*nfac == *n)
    {
      goto L50;
    }
  nfac1 = *nfac + 1;
  nnfac = *n - *nfac;
  n2fac = *nfac * nfac1 / 2;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L10: */
      w[i__] = g1[i__] * *sig;
    }
  k = n2fac;
  if (*nfac == 0)
    {
      goto L25;
    }
  i__1 = *nfac;
  for (j = 1; j <= i__1; ++j)
    {
      i__2 = *n;
      for (i__ = nfac1; i__ <= i__2; ++i__)
	{
	  ++k;
	  dh[k] += g1[i__] * w[j];
	  /* L20: */
	}
    }
 L25:
  k = n2fac + *nfac * nnfac;
  i__2 = *n;
  for (j = nfac1; j <= i__2; ++j)
    {
      i__1 = *n;
      for (i__ = j; i__ <= i__1; ++i__)
	{
	  ++k;
	  dh[k] += g1[i__] * w[j];
	  /* L30: */
	}
    }
 L50:
  *ir = *nfac;
  if (*nfac == 0)
    {
      return 0;
    }
  optim_majour (&dh[1], &g1[1], &w[1], nfac, sig, ir, mk, epsmc);
  return 0;
}	

