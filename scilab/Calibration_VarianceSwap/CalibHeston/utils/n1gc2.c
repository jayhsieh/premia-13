#include "optim.h"

static int optim_n1gc2a (opt_simul simul, opt_prosca prosca, int *n, double *x, double *f,
			 double *g, double *dx, double *df1, double *eps, int *imp,
			 int *io, int *niter, int *nsim, int *info, int *memh,
			 double *d__, double *xx, double *gg, double *tabaux,
			 double *h__, opt_simul_data *optim_data);

static int optim_n1gc2b (int *n, opt_simul simul, opt_prosca prosca, double *xinit, double *f,
			 double *dg, double *alpha, double *d__, double *xfinal,
			 double *gfinal, int *imp, int *io, int *retour, int *ntotap,
			 int *nsim, int *intfor, double *dx, double *eps, opt_simul_data *optim_data);

static int optim_fmuls1 (int *n, double *h__, double *x, double *hx);

static int optim_fmulb1 (int *n, double *h__, double *x, double *hx, double *tabaux,
			 int *nmisaj, opt_prosca prosca, opt_simul_data *optim_data);

/*!but 
 *    minimisation sans contraintes par un algorithme de quasi-Newton 
 *    a memoire limitee 
 *!origine 
 *    c. lemarechal, inria, 1987 
 *    Copyright INRIA 
 *!commentaires 
 *le sous-programme n1gc2 (gradient conjugue a encombrement variable) 
 *est une interface entre le programme appelant et le sous-programme 
 *n1gc2a, minimiseur proprement dit. 
 *nrz est la dimension declaree pour le tableau de travail rz. 
 *rz est subdivise en 4 vecteurs de dimension n 
 *                    et un tableau de dimension memh. 
 *memh est la dimension allouee a la matrice de quasi newton h. 
 *pour l'usage de la liste d'appel : voir la documentation de n1qn1 
 *! 
 *declaration des tableaux 
 *declaration des scalaires 
 * 
 * 
 */

int optim_n1gc2 (opt_simul simul, opt_prosca prosca, int *n, double *x, double *f, double *g,
		 double *dxmin, double *df1, double *epsrel, int *imp, int *io,
		 int *mode, int *niter, int *nsim, double *rz, int *nrz,
		 opt_simul_data *optim_data)
{
  /* Local variables */
  int memh, iaux;
  int id, ig, ih, ix;

  /* Parameter adjustments */
  --g;
  --x;
  --rz;

  /* Function Body */
  if (*imp > 0)
    {
      Sciprintf("n1gc2: dimension du probleme %d\n",*n);
      Sciprintf("\tnrz=%d, niter=%d, nsim=%d, imp=%d\n",*nrz,*niter,*nsim,*imp);
      Sciprintf("\tepsrel=%8.2f, df1=%8.2f, dxmin=%8.2f\n",
		(*epsrel),
		(*df1),
		(*dxmin));
    }
  if (*n <= 0 || *niter <= 0 || *nsim <= 0 || *dxmin <= 0. || *df1 <= 0.
      || *epsrel <= 0. || *epsrel > 1.)
    {
      *mode = 2;
      if (*imp > 0)
	{
	  Sciprintf("n1gc2: appel incoherent\n");
	}
      return 0;
    }
  /* 
   *calculs des pointeurs destines a subdiviser le tableau rz 
   */
  id = 1;
  ix = id + *n;
  ig = ix + *n;
  iaux = ig + *n;
  ih = iaux + *n;
  /* 
   *calcul du nombre de places memoire affectees a h 
   */
  memh = *nrz - (*n << 2);
  /* 
   */
  if (memh <= 0)
    {
      *mode = 3;
      goto L100;
    }
  else
    {
    }
  /* 
   *appel du sous-programme n1gc2a qui effectue la reelle optimisation 
   */
  optim_n1gc2a ( simul, prosca, n, &x[1], f, &g[1], dxmin, df1,
		 epsrel, imp, io, niter, nsim, mode, &memh, &rz[id], &rz[ix],
		 &rz[ig], &rz[iaux], &rz[ih],optim_data);
  /* 
   */
 L100:
  if (*imp > 0)
    {
      if (*mode == 3)
	{
	  Sciprintf( "n1gc2:  rz insuffisamment dimensionne\n");
	}
      else if (*mode == 6)
	{ 
	  Sciprintf( "n1gc2:   fin sur dxmin\n");
	}
      else
	{  
	  Sciprintf("n1gc2: end with norme of g =%15.9f, niter=%d, nsim=%d\n",
		    (*epsrel),
		    (*niter),
		    (*nsim));
	}
    }
  return 0;
}




static int optim_n1gc2a (opt_simul simul, opt_prosca prosca, int *n, double *x, double *f,
			 double *g, double *dx, double *df1, double *eps, int *imp,
			 int *io, int *niter, int *nsim, int *info, int *memh,
			 double *d__, double *xx, double *gg, double *tabaux,
			 double *h__,  opt_simul_data *optim_data)
{

  /* System generated locals */
  int i__1, i__2;
  double d__1;
  /* Local variables */
  int ieta, iter, i__, j, k, m;
  int l;
  double alpha, omega;
  int redem;
  double sigma;
  int termi;
  double normg;
  double normg0;
  int gc;
  double dg;
  int kj, lk, is, iu;
  double mu, gcarre, ggcarr, nu, sscaek, sscalg, uscalg;
  int nmisaj;
  int redfor;
  double dg1;
  int memuti;
  int intfor, iterqn;
  int ntotap, memsup, km1, kp1, retour, nrzuti;
  double eta;
  int inu;
  double aux1, aux2;

  /*    Copyright INRIA 
   * 
   *parametres 
   *declaration des tableaux 
   *declaration des scalaires 
   * 
   * 
   * ************************************************************* 
   * phase i:determination de la methode ( et de m le cas echeant) 
   * ************************************************************* 
   * 
   */
  /* Parameter adjustments */
  --tabaux;
  --gg;
  --xx;
  --d__;
  --g;
  --x;
  --h__;

  /* Function Body */
  memuti = *n * (*n + 1) / 2;
  /* 
   *memsup est aussi la dimension minimale de la matrice h 
   */
  memsup = (*n << 1) + 2;
  /* 
   */
  if (*memh >= memuti)
    {
      gc = FALSE;
      nrzuti = memuti + (*n << 2);
      if (*imp > 1)
	{
	  Sciprintf("methode de quasi-newton. nrz utile=%d\n",nrzuti);
	}
    }
  else if (*memh < memsup)
    {
      *info = 3;
      return 0;
    }
  else
    {
      gc = TRUE;
      /*m est le nombre de mises a jour admissible 
       */
      m = *memh / memsup;
      /*memuti est ici le nombre de places memoire utilisees pour stocker h 
       */
      memuti = m * memsup;
      nrzuti = memuti + (*n << 2);
      if (*imp > 1)
	{
	  Sciprintf("methode du gradient conjugue avec,%d, mises a jour. nrz utile=%d",m,nrzuti);
  	}
    }
  /* 
   * *********************************************** 
   * phase ii:initialisations propres a l'optimiseur 
   * *********************************************** 
   * 
   *initialisation des compteurs 
   */
  iter = 0;
  ntotap = 1;
  /* 
   * ****************************************************************** 
   * phase iii:demarrage a partir de x(0) avec descente suivant -(grad) 
   * ****************************************************************** 
   * 
   */
 L3000:
  i__ = 0;
  nmisaj = 0;
  /* 
   *calcul de la direction de descente 
   */
  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
    {
      d__[j] = -g[j];
      /* L3100: */
    }
  /* 
   */
  (*prosca) (n, &g[1], &d__[1], &dg1,optim_data);
  normg0 = sqrt ((Abs (dg1)));
  if (iter == 1)
    {
      omega = *eps * normg0;
    }
  /* 
************************************************************* 
*phase iv:debut de l'iteration x(i-1) donne x(i) le long de d 
************************************************************* 
*/
 L4000:
  if (iter == *niter)
    {
      *info = 4;
      goto L99999;
    }
  ++iter;
  ++i__;
  /* 
   *determination du type de pas 
   */
  if (gc)
    {
      iterqn = i__ <= m && 2 <= i__;
    }
  /* 
   * ******************************* 
   * phase v:initialisation de alpha 
   * ******************************* 
   * 
   */
  if (iter == 2)
    {
      alpha = *df1 * 2. / (-dg1);
    }
  else if (gc)
    {
      if (i__ == 1)
	{
	  alpha = 1. / normg0;
	}
      else
	{
	  if (iterqn)
	    {
	      alpha = 1.;
	    }
	  else
	    {
	      alpha = alpha * dg / dg1;
	    }
	}
    }
  else
    {
      alpha = 1.;
    }
  /* 
**************************** 
*phase vi:recherche lineaire 
**************************** 
* 
*/
  dg = dg1;
  intfor = (gc && !iterqn) || (!gc && i__ == 1);
  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
    {
      xx[j] = x[j];
      gg[j] = g[j];
      /* L6000: */
    }
  optim_n1gc2b (n,  simul,  prosca, &xx[1], f, &dg, &alpha,
		&d__[1], &x[1], &g[1], imp, io, &retour, &ntotap, nsim,
		&intfor, dx, eps,optim_data);
  /* 
   */
  if (*imp > 3)
    {
      Sciprintf("()\n");
    }
  if (retour == 4 || (retour == 1 && i__ == 1))
    {
      *info = 6;
      return 0;
    }
  else if (retour == 1)
    {
      if (*imp > 1)
	{
	  Sciprintf("n1gc2 %d, iters %d simuls, necessite d'un redemarrage total",iter,ntotap);
	}
      goto L3000;
    }
  else
    {
      /*calcul de (g,g) 
       */
      if (i__ > 1 && gc)
	{
	  ggcarr = gcarre;
	}
      (*prosca) (n, &g[1], &g[1], &gcarre,optim_data);
      normg = sqrt (gcarre);
      if (*imp > 2)
	{ 
	  Sciprintf("n1gc2 %d, iters %d simuls f=%15.9f",iter,ntotap,*f);
	}
      if (retour == 2)
	{
	  *info = 0;
	  goto L99999;
	}
      else if (retour == 3)
	{
	  *info = 5;
	  goto L99999;
	}
    }
  /* 
   * ****************************************************** 
   * phase vii:test d'arret par obtention de la convergence 
   * ****************************************************** 
   * 
   */
  termi = normg < omega;
  if (termi)
    {
      *info = 1;
      goto L99999;
    }
  else
    {
    }
  /* 
   * ******************************************* 
   * phase viii:test x(i) point de redemarrage? 
   ******************************************** 
   * 
   *doit on forcer un redemarrage? 
   */
  redfor = gc && (i__ == 1 || i__ == m + *n);
  if (redfor)
    {
      redem = TRUE;
    }
  else if (gc && !iterqn)
    {
      (*prosca) (n, &g[1], &gg[1], &aux1,optim_data);
      redem = Abs (aux1) > (d__1 = ggcarr * .2, Abs (d__1));
    }
  else
    {
      redem = FALSE;
    }
  /* 
   * ******************** 
   * phase ix:mise a jour 
   * ******************** 
   * 
   *calcul de s stocke dans d et de y stocke dans xx 
   */
  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
    {
      d__[j] = alpha * d__[j];
      xx[j] = g[j] - gg[j];
      /* L9000: */
    }
  if (redem)
    {
      /*cas ou x(i) est un point de redemarrage 
       */
      i__ = 1;
      nmisaj = 1;
      /*sauvegarde de s qui est actuellement dans d 
       *               u=h*y=y 
       *               nu=(y,hy)=(y,y) 
       *               eta=(s,y) 
       *calcul des indices 
       */
      inu = 1;
      ieta = inu + 1;
      iu = ieta;
      is = iu + *n;
      /* 
       */
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
	{
	  h__[iu + j] = xx[j];
	  h__[is + j] = d__[j];
	  /* L9100: */
	}
      (*prosca) (n, &xx[1], &xx[1], &nu,optim_data);
      h__[inu] = nu;
      (*prosca) (n, &d__[1], &xx[1], &eta,optim_data);
      h__[ieta] = eta;
      /*h1 est maintenant definie 
       *calcul de h1*g que l'on range dans xx 
       */
      optim_fmulb1 (n, &h__[1], &g[1], &xx[1], &tabaux[1], &nmisaj,
		    prosca,optim_data);
      /* 
       */
    }
  else if (gc)
    {
      /*cas de gc sans redamarrage 
       *calcul de h*y range dans gg 
       */
      optim_fmulb1 (n, &h__[1], &xx[1], &gg[1], &tabaux[1], &nmisaj,
		    prosca, optim_data);
      /*calculs de  nu, eta, sscalg, uscalg 
       */
      (*prosca) (n, &xx[1], &gg[1], &nu,optim_data);
      (*prosca) (n, &d__[1], &xx[1], &eta,optim_data);
      (*prosca) (n, &d__[1], &g[1], &sscalg,optim_data );
      (*prosca) (n, &gg[1], &g[1], &uscalg,optim_data);
      /*calcul de sigma et de mu 
       */
      sigma = (uscalg - (nu / eta + 1.) * sscalg) / eta;
      mu = sscalg / eta;
      /*calcul de h*g que l'on range dans xx 
       */
      optim_fmulb1 (n, &h__[1], &g[1], &xx[1], &tabaux[1], &nmisaj,
		    prosca, optim_data);
      /*calcul de la nouvelle direction de recherche: 
       *h*g - mu * u - sigma * s 
       */
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
	{
	  xx[j] = xx[j] - mu * gg[j] - sigma * d__[j];
	  /* L9200: */
	}
      /* 
       *cas d'une iteration de type quasi newton 
       */
      if (iterqn)
	{
	  ++nmisaj;
	  /*sauvegarde des termes utiles pour stocker la matrice mise a jour 
	   */
	  inu += memsup;
	  ieta = inu + 1;
	  iu = ieta;
	  is = iu + *n;
	  i__1 = *n;
	  for (j = 1; j <= i__1; ++j)
	    {
	      h__[iu + j] = gg[j];
	      h__[is + j] = d__[j];
	      /* L9300: */
	    }
	  h__[inu] = nu;
	  h__[ieta] = eta;
	}
      /*cas de la methode quasi newton 
       */
    }
  else
    {
      /*calcul de eta=(s,y) 
       */
      (*prosca) (n, &d__[1], &xx[1], &eta,optim_data);
      if (i__ == 1)
	{
	  /*etape initiale calcul de l'approximation initiale de l'inverse de la 
	   *matrice hessienne 
	   *calcul de nu=(y,h0*y)=(y,y) 
	   */
	  (*prosca) (n, &xx[1], &xx[1], &nu,optim_data);
	  /*stockage de cette matrice h=(eta / nu) * i 
	   */
	  kj = 1;
	  aux1 = eta / nu;
	  i__1 = *n;
	  for (k = 1; k <= i__1; ++k)
	    {
	      h__[kj] = aux1;
	      ++kj;
	      kp1 = k + 1;
	      if (*n >= kp1)
		{
		  i__2 = *n;
		  for (j = kp1; j <= i__2; ++j)
		    {
		      h__[kj] = 0.;
		      ++kj;
		      /* L9400: */
		    }
		}
	      gg[k] = aux1 * xx[k];
	      /* L9500: */
	    }
	  nu = eta;
	}
      else
	{
	  optim_fmuls1 (n, &h__[1], &xx[1], &gg[1]);
	  (*prosca) (n, &xx[1], &gg[1], &nu,optim_data);
	}
      /*calcul de la matrice mise a jour (utilisation de la formule bfgs ) 
       *nu, eta et h*y (stocke dans gg) sont connus 
       */
      aux1 = nu / eta + 1.;
      kj = 1;
      i__1 = *n;
      for (k = 1; k <= i__1; ++k)
	{
	  /*calcul du vecteur contenant la keme colonne de h 
	   */
	  lk = k;
	  km1 = k - 1;
	  if (k >= 2)
	    {
	      i__2 = km1;
	      for (l = 1; l <= i__2; ++l)
		{
		  tabaux[l] = h__[lk];
		  lk += *n - l;
		  /* L9610: */
		}
	    }
	  i__2 = *n;
	  for (l = k; l <= i__2; ++l)
	    {
	      tabaux[l] = h__[lk];
	      ++lk;
	      /* L9620: */
	    }
	  /* 
	   */
	  (*prosca) (n, &xx[1], &tabaux[1], &aux2,optim_data);
	  i__2 = *n;
	  for (l = 1; l <= i__2; ++l)
	    {
	      tabaux[l] = 0.;
	      /* L9630: */
	    }
	  tabaux[k] = 1.;
	  (*prosca) (n, &tabaux[1], &d__[1], &sscaek,optim_data);
	  kj = k - *n;
	  i__2 = k;
	  for (j = 1; j <= i__2; ++j)
	    {
	      kj = kj + *n - j + 1;
	      h__[kj] -=
		((aux2 - aux1 * sscaek) * d__[j] + sscaek * gg[j]) / eta;
	      /* L9700: */
	    }
	  /* L9800: */
	}
    }
  /* 
   * ***************************************************** 
   * phase x :calcul de la nouvelle direction de recherche 
   * ***************************************************** 
   * 
   */
  if (gc)
    {
      /*xx contient -d 
       */
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
	{
	  d__[j] = -xx[j];
	  /* L10000: */
	}
      /* 
       */
    }
  else
    {
      /*cas de la methode de quasi newton 
       *la nouvelle direction d egale -(h * g) 
       */
      optim_fmuls1 (n, &h__[1], &g[1], &d__[1]);
      i__1 = *n;
      for (j = 1; j <= i__1; ++j)
	{
	  d__[j] = -d__[j];
	  /* L10100: */
	}
    }
  /* 
   *test:la direction de recherche est elle bien de descente 
   */
  (*prosca) (n, &d__[1], &g[1], &dg1,optim_data);
  if (dg1 >= 0.)
    {
      *info = 7;
      if (*imp > 1)
	{ 
	  Sciprintf("n1gc2a: erreur dans la hessienne dg=%9.2f\n",dg1);
	}
      goto L99999;
    }
  else
    {
      goto L4000;
    }
  /* 
   *retour au programme appelant 
   */
 L99999:
  *niter = iter;
  *nsim = ntotap;
  if (i__ == 0)
    {
      *eps = normg0;
    }
  else
    {
      *eps = normg;
    }
  return 0;
}


static int optim_fmulb1 (int *n, double *h__, double *x, double *hx, double *tabaux,
			 int *nmisaj, opt_prosca prosca, opt_simul_data *optim_data)
{
  /* System generated locals */
  int i__1;

  /* Local variables */
  int ptnu, k;
  double gamma, sigma;
  int compt, is, iu;
  double mu, nu, sscalx, uscalx;
  int memsup;
  double eta;

  /*    Copyright INRIA 
   * 
   *parametres 
   *declarations 
   * 
   *calcul de la longueur d'un bloc 
   */
  /* Parameter adjustments */
  --tabaux;
  --hx;
  --x;
  --h__;

  /* Function Body */
  memsup = (*n << 1) + 2;
  /*calcul de  h0*x=x=x 
   */
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      hx[k] = x[k];
      /* L1000: */
    }
  /* 
   */
  if (*nmisaj == 0)
    {
      return 0;
    }
  else
    {
      ptnu = 1;
      compt = 1;
    }
  /* 
   */
 L2000:
  iu = ptnu + 1;
  is = iu + *n;
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      tabaux[k] = h__[iu + k];
      /* L3000: */
    }
  (*prosca) (n, &tabaux[1], &x[1], &uscalx,optim_data);
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      tabaux[k] = h__[is + k];
      /* L4000: */
    }
  (*prosca) (n, &tabaux[1], &x[1], &sscalx,optim_data);
  nu = h__[ptnu];
  eta = h__[ptnu + 1];
  /*calcul du nouveau terme et addition dans hx 
   */
  if (compt == 1)
    {
      gamma = eta / nu;
      i__1 = *n;
      for (k = 1; k <= i__1; ++k)
	{
	  hx[k] = gamma * hx[k];
	  /* L5000: */
	}
      mu = sscalx / nu;
      sigma = -(sscalx * 2. / eta) + uscalx / nu;
    }
  else
    {
      mu = sscalx / eta;
      sigma = -(nu / eta + 1.) * mu + uscalx / eta;
    }
  /* 
   */
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      hx[k] = hx[k] - mu * h__[iu + k] - sigma * h__[is + k];
      /* L6000: */
    }
  /* 
   */
  ++compt;
  if (compt <= *nmisaj)
    {
      ptnu += memsup;
      goto L2000;
    }
  else
    {
      return 0;
    }
  return 0;
}				/* fmulb1_ */

/* 
 *ce sous-programme effectue le produit  h * x   avec: 
 *n (e) dimension du probleme 
 *h (e) dimension n(n+1)/2. tiangle superieur, coefficients par ligne 
 *x (e) vecteur de dimension n 
 *hx (s) dimension n. resultat du produit 
 * 
 *    Copyright INRIA 
 * 
 *parametre 
 *declarations 
 * 
 */

static int optim_fmuls1 (int *n, double *h__, double *x, double *hx)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int j, k, kj, km1;
  double aux1;

  /* Parameter adjustments */
  --hx;
  --x;
  --h__;

  /* Function Body */
  i__1 = *n;
  for (k = 1; k <= i__1; ++k)
    {
      /*calcul de la keme composante du produit  h* x 
       */
      aux1 = 0.;
      /*h(kj) est le coefficient (k,j) de la matrice symetrique complete 
       */
      kj = k;
      km1 = k - 1;
      /*contribution du triangle inferieur 
       */
      if (km1 >= 1)
	{
	  i__2 = km1;
	  for (j = 1; j <= i__2; ++j)
	    {
	      aux1 += h__[kj] * x[j];
	      kj += *n - j;
	      /* L1000: */
	    }
	}
      /*contribution du triangle superieur 
       */
      i__2 = *n;
      for (j = k; j <= i__2; ++j)
	{
	  aux1 += h__[kj] * x[j];
	  ++kj;
	  /* L2000: */
	}
      /* 
       */
      hx[k] = aux1;
      /* L3000: */
    }
  /* 
   */
  return 0;
}				/* fmuls1_ */



static int optim_n1gc2b (int *n, opt_simul simul, opt_prosca prosca, double *xinit, double *f,
			 double *dg, double *alpha, double *d__, double *xfinal,
			 double *gfinal, int *imp, int *io, int *retour, int *ntotap,
			 int *nsim, int *intfor, double *dx, double *eps, opt_simul_data *optim_data)
{
  int i__1;
  double d__1;
  double bsup;
  int j, indic;
  double delta;
  int depas;
  double finit, ap, dp, at, fp;
  int encadr, accept, rfinie;
  int nappel;
  int maxpas;
  double dal, pas, aux1, aux2;

  /* Parameter adjustments */
  --gfinal;
  --xfinal;
  --d__;
  --xinit;


  /* Function Body */
  depas = FALSE;
  bsup = 0.;
  finit = *f;
  nappel = 0;
  ap = 0.;
  fp = finit;
  dp = *dg;
  if (*imp > 3)
    {
      Sciprintf("n1gc2b,6x,  pas,%10.3f,  dg=,%9.2f",*alpha,*dg);
    }
  /*calcul de la longueur du pas 
   */
  (*prosca) (n, &d__[1], &d__[1], &pas,optim_data);
  pas = sqrt (pas);
  /*test d'erreur dans la recherche lineaire 
   */
 L1000:
  if (*alpha * pas <= *dx)
    {
      if (*imp > 3)
	{ 
	  Sciprintf("n1gc2b fin sur dx\n");
	}
      *retour = 1;
      return 0;
    }
  else if (*ntotap == *nsim)
    {
      *retour = 3;
      return 0;
    }
  else
    {
    }
  /*calcul du nouveau point susceptible d'etre xfinal 
   */
  i__1 = *n;
  for (j = 1; j <= i__1; ++j)
    {
      xfinal[j] = xinit[j] + *alpha * d__[j];
      /* L2000: */
    }
  /*calculs de f et g en ce point 
   */
  indic = 4;
  (*simul) (&indic, n, &xfinal[1], f, &gfinal[1],optim_data);
  ++nappel;
  ++(*ntotap);
  if (indic < 0)
    {
      depas = TRUE;
      if (*imp > 3)
	{ 
	  Sciprintf("n1gc2b, %10.3f, indic=%d\n",*alpha,indic);
	}
      delta = *alpha - ap;
      if (delta <= *dx)
	{
	  *retour = 4;
	  return 0;
	}
      else
	{
	  bsup = *alpha;
	  *alpha = delta * .1 + ap;
	  goto L1000;
	}
    }
  /*calcul de la derivee suivant d au point xfinal 
   */
  (*prosca) (n, &d__[1], &gfinal[1], &dal,optim_data);
  /* 
   */
  if (*imp > 3)
    { 
      Sciprintf("n1gc2b, %10.3f, %11.3f %11.3f\n",*alpha,aux2,dal);
    }
  if (indic == 0)
    {
      *retour = 2;
      return 0;
    }
  maxpas = *f > finit && dal < 0.;
  if (maxpas)
    {
      *alpha /= 3.;
      ap = 0.;
      fp = finit;
      dp = *dg;
      rfinie = FALSE;
      /* 
       */
    }
  else
    {
      /*test:le nouveau point est il acceptable 
       */
      aux1 = finit + *alpha * 1e-4 * *dg;
      aux2 = (d__1 = dal / *dg, Abs (d__1));
      accept = *f <= aux1 && aux2 <= .9;
      if (accept)
	{
	  /*doit on faire une interpolation 
	   */
	  rfinie = nappel > 1 || !(*intfor) || aux2 <= *eps;
	}
      else
	{
	  rfinie = FALSE; 
	}
      /* 
       */
      if (!rfinie)
	{
	  /*la recherche lineaire n'est pas finie. interpolation 
	   */
	  aux1 = dp + dal - (fp - *f) * 3. / (ap - *alpha);
	  aux2 = aux1 * aux1 - dp * dal;
	  if (aux2 <= 0.)
	    {
	      aux2 = 0.;
	    }
	  else
	    {
	      aux2 = sqrt (aux2);
	    }
	  if (dal - dp + aux2 * 2. == 0.)
	    {
	      *retour = 4;
	      return 0;
	    }
	  at =
	    *alpha - (*alpha - ap) * (dal + aux2 - aux1) / (dal - dp +
							    aux2 * 2.);
	  /*test:le minimum a t-il ete encadre 
	   */
	  encadr = dal / dp <= 0.;
	  if (encadr)
	    {
	      /*le minimum a ete encadre 
	       */
	      if ((d__1 = *alpha - ap, Abs (d__1)) <= *dx)
		{
		  *retour = 4;
		  return 0;
		}
	      aux1 = Min (*alpha, ap) * 1.01;
	      aux2 = Max (*alpha, ap) * .99;
	      if (at < aux1 || at > aux2)
		{
		  at = (*alpha + ap) / 2.;
		}
	    }
	  else
	    {
	      /*le minimum n'a pas ete encadre 
	       */
	      aux1 = Min (ap, *alpha) * .99;
	      if (dal <= 0. || at <= 0. || at >= aux1)
		{
		  aux1 = Max (ap, *alpha) * 1.01;
		  if (dal > 0. || at <= aux1)
		    {
		      if (dal <= 0.)
			{
			  at = Max (ap, *alpha) * 2.;
			}
		      else
			{
			  at = Min (ap, *alpha) / 2.;
			}
		    }
		}
	    }
	  if (depas && at >= bsup)
	    {
	      delta = bsup - *alpha;
	      if (delta <= *dx)
		{
		  *retour = 4;
		  return 0;
		}
	      else
		{
		  at = *alpha + delta * .1;
		}
	    }
	  ap = *alpha;
	  fp = *f;
	  dp = dal;
	  *alpha = at;
	}
    }
  if (rfinie)
    {
      *retour = 0;
      return 0;
    }
  else
    {
      goto L1000;
    }
  return 0;
}	

