#include "optim.h"

static int optim_n1qn1a (opt_simul simul, int *n, double *x, double *f, double *g,
			 double *scale, double *acc, int *mode, int *niter, int *nsim,
			 int *iprint, int *lp, double *h__, double *d__, double *w,
			 double *xa, double *ga, double *xb, double *gb,
			 opt_simul_data *optim_data);

/*!but 
 *    minimisation d une fonction reguliere sans contraintes 
 *!origine 
 *    c. lemarechal, inria, 1987 
 *    Copyright INRIA 
 *!methode 
 *    direction de descente calculee par une methode de quasi-newton 
 *    recherche lineaire de type wolfe 
 *!liste d appel 
 *    simul    : point d'entree au module de simulation (cf normes modulopt i) 
 *    n1qn1 appelle toujours simul avec indic = 4 ; le module de 
 *    simulation doit se presenter sous la forme subroutine simul 
 *    (n,x, f, g, izs, rzs, dzs) et e^tre declare en external dans le 
 *    programme appelant n1qn1. 
 *    n (e)    : nombre de variables dont depend f. 
 *    x (e-s)   : vecteur de dimension n ; en entree le point initial ; 
 *                en sortie : le point final calcule par n1qn1. 
 *    f (e-s)   : scalaire ; en entree valeur de f en x (initial), en sortie 
 *                valeur de f en x (final). 
 *    g (e-s)   : vecteur de dimension n : en entree valeur du gradient en x 
 *                (initial), en sortie valeur du gradient en x (final). 
 *    var (e)   : vecteur strictement positif de dimension n. amplitude de la 
 *                modif souhaitee a la premiere iteration sur x(i).une bonne 
 *                valeur est 10% de la difference (en valeur absolue) avec la 
 *                coordonee x(i) optimale 
 *    eps (e-s) : en entree scalaire definit la precision du test d'arret. 
 *     le programme considere que la convergence est obtenue lorque il lui 
 *     est impossible de diminuer f en attribuant a au moins une coordonnee 
 *     x(i) une variation superieure a eps*var(i). 
 *     en sortie, eps contient le carre de la norme du gradient en x (final). 
 *    mode (e)     : definit l approximation initiale du hessien 
 *                 =1 n1qn1 l initialise lui-meme 
 *                 =2 le hessien est fourni dans zm sous forme compressee (zm 
 *                    contient les colonnes de la partie inferieure du hessien) 
 *    niter (e-s)  : en entree nombre maximal d'iterations : en sortie nombre 
 *                   d'iterations reellement effectuees. 
 *    nsim (e-s)  : en entree nombre maximal d'appels a simul (c'est a dire 
 *        avec indic = 4). en sortie le nombre de tels appels reellement faits. 
 *     imp (e)   : contro^le les messages d'impression : 
 *                 0 rien n'est imprime 
 *                 = 1 impressions initiales et finales 
 *                 = 2 une impression par iteration (nombre d'iterations, 
 *                     nombre d'appels a simul, valeur courante de f). 
 *                 >=3 informations supplementaires sur les recherches 
 *                     lineaires ; 
 *                     tres utile pour detecter les erreurs dans le gradient. 
 *     lp (e)    : le numero du canal de sortie, i.e. les impressions 
 *                 commandees par imp sont faites par write (lp, format). 
 *    zm     : memoire de travail pour n1qn1 de   dimension n*(n+13)/2. 
 *    izs,rzs,dzs memoires reservees au simulateur (cf doc) 
 * 
 *! 
 */



int optim_n1qn1 (opt_simul simul, int *n, double *x, double *f, double *g, double *var,
		 double *eps, int *mode, int *niter, int *nsim, int *imp, int *lp,
		 double *zm, opt_simul_data *optim_data)
{
  int nd, nw, nga, ngb, nxa, nxb,ret;
  /* Parameter adjustments */
  --var;
  --g;
  --x;
  --zm;

  if (*imp > 0)
    {
      Sciprintf("n1qn1: optimization without bound constraints\n");
      Sciprintf("dimension=%d, epsg=%g, verbosity level: imp=%d\n",*n,*eps,*imp);
      Sciprintf("maximum number of iterations allowed: iter=%d\n",*niter);
      Sciprintf("maximum number of calls to costf allowed: nap=%d\n",*nsim);
      Sciprintf("------------------------------------------------\n");
    }
  nd = *n * (*n + 1) / 2 + 1;
  nw = nd + *n;
  nxa = nw + *n;
  nga = nxa + *n;
  nxb = nga + *n;
  ngb = nxb + *n;
  ret = optim_n1qn1a ( simul, n, &x[1], f, &g[1], &var[1], eps, mode, niter,
		       nsim, imp, lp, &zm[1], &zm[nd], &zm[nw], &zm[nxa], &zm[nga],
		       &zm[nxb], &zm[ngb],optim_data);
  if (*imp > 0)
    {
      Sciprintf("n1qn1: end of optimization, gradient norm=%g\n",sqrt (*eps));
    }
  return ret;
}


/*    Copyright INRIA 
 *    A (very) few modifs by Bruno (14 March 2005): I have translated some output 
 *    informations in english (but they don't use format instruction 
 *    which is put in the secong arg of write). Also for the linear 
 *    search output informations I divide by the direction vector norm 
 *    to get the "normalized" directionnal derivative. Note that this is 
 *    just for output (the computing code is normaly not modified). 
 *a better information concerning direction 
 *(blas routine) added by Bruno to get 
 * 
 * calcul initial de fonction-gradient 
 * 
 */

static int optim_n1qn1a (opt_simul simul, int *n, double *x, double *f, double *g,
			 double *scale, double *acc, int *mode, int *niter, int *nsim,
			 int *iprint, int *lp, double *h__, double *d__, double *w,
			 double *xa, double *ga, double *xb, double *gb,
			 opt_simul_data *optim_data)
{
  int ret =OK;
  int c__1 = 1;
  double c_b101 = 0.;
  int i__1, i__2, i__3;
  double d__1, d__2, d__3, d__4;

  /* Local variables */
  double fmin, gmin;
  int nfun, isfv;
  double step;
  double c__;
  int i__, j, k;
  double s;
  int indic;
  double v;
  int iecri, i1;
  double stmin, cc, fa, fb, hh;
  int ii, ij, ik, jk, ni, ip, ir, np;
  double stepbd, steplb;
  double gl1, gl2, dga, dgb, dff;
  int ial, nip, itr;

  /* Parameter adjustments */
  --gb;
  --xb;
  --ga;
  --xa;
  --w;
  --d__;
  --scale;
  --g;
  --x;
  --h__;

  /* Function Body */
  indic = 4;
  (*simul) (&indic, n, &x[1], f, &g[1],optim_data);
  if (indic <= 0)
    {
      if (*iprint != 0)
	{
	  if (indic < 0)
	    {
	      Sciprintf("n1qn1 ne peut demarrer (contrainte implicite)\n");
	    }
	  else if (indic == 0)
	    {  
	      Sciprintf("n1qn1 termine par voeu de l'utilisateur\n");
	    }
	}
      if (indic == 0) ret = FAIL;
      *acc = 0.;
      *niter = 1;
      *nsim = 1;
      return ret;
    }

  nfun = 1;
  iecri = 0;
  itr = 0;
  np = *n + 1;
  /*                 initialisation du hessien, en fonction de var 
   */
  if (*mode >= 2)
    {
      goto L60;
    }
 L20:
  c__ = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L30: */
      /*Computing MAX 
       */
      d__2 = c__, d__3 = (d__1 = g[i__] * scale[i__], Abs (d__1));
      c__ = Max (d__2, d__3);
    }
  if (c__ <= 0.)
    {
      c__ = 1.;
    }
  k = *n * np / 2;
  i__1 = k;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L40: */
      h__[i__] = 0.;
    }
  k = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      h__[k] = c__ * .01 / (scale[i__] * scale[i__]);
      /* L50: */
      k = k + np - i__;
    }
  goto L100;
  /*              factorisation du hessien 
   */
 L60:
  if (*mode >= 3)
    {
      goto L80;
    }
  k = *n;
  if (*n > 1)
    {
      goto L300;
    }
  if (h__[1] > 0.)
    {
      goto L305;
    }
  h__[1] = 0.;
  k = 0;
  goto L305;
 L300:
  np = *n + 1;
  ii = 1;
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      hh = h__[ii];
      ni = ii + np - i__;
      if (hh > 0.)
	{
	  goto L301;
	}
      h__[ii] = 0.;
      --k;
      ii = ni + 1;
      goto L304;
    L301:
      ip = ii + 1;
      ii = ni + 1;
      jk = ii;
      i__2 = ni;
      for (ij = ip; ij <= i__2; ++ij)
	{
	  v = h__[ij] / hh;
	  i__3 = ni;
	  for (ik = ij; ik <= i__3; ++ik)
	    {
	      h__[jk] -= h__[ik] * v;
	      /* L302: */
	      ++jk;
	    }
	  /* L303: */
	  h__[ij] = v;
	}
    L304:
      ;
    }
  if (h__[ii] > 0.)
    {
      goto L305;
    }
  h__[ii] = 0.;
  --k;
 L305:
  /* 
   */
  if (k >= *n)
    {
      goto L100;
    }
 L70:
  if (*iprint != 0)
    { 
      Sciprintf("n1qn1 remplace le hessien initial (qui n'est, pas defini positif)/ par une diagonale positive\n");
    }
  goto L20;
  /*               verification que la diagonale est positive 
   */
 L80:
  k = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (h__[k] <= 0.)
	{
	  goto L70;
	}
      /* L90: */
      k = k + np - i__;
    }
  /*               quelques initialisations 
   */
 L100:
  dff = 0.;
 L110:
  fa = *f;
  isfv = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      xa[i__] = x[i__];
      ga[i__] = g[i__];
    }
  /*                  iteration 
   */
 L130:
  ++itr;
  ial = 0;
  if (itr > *niter)
    {
      goto L250;
    }
  ++iecri;
  if (iecri != -(*iprint))
    {
      goto L140;
    }
  iecri = 0;
  indic = 1;
  (*simul) (&indic, n, &x[1], f, &g[1],optim_data);
  /*              calcul de la direction de recherche 
   */
 L140:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L150: */
      d__[i__] = -ga[i__];
    }
  w[1] = d__[1];
  if (*n > 1)
    {
      goto L400;
    }
  d__[1] /= h__[1];
  goto L412;
 L400:
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      ij = i__;
      i1 = i__ - 1;
      v = d__[i__];
      i__2 = i1;
      for (j = 1; j <= i__2; ++j)
	{
	  v -= h__[ij] * d__[j];
	  /* L401: */
	  ij = ij + *n - j;
	}
      w[i__] = v;
      /* L402: */
      d__[i__] = v;
    }
  d__[*n] /= h__[ij];
  np = *n + 1;
  i__1 = *n;
  for (nip = 2; nip <= i__1; ++nip)
    {
      i__ = np - nip;
      ii = ij - nip;
      v = d__[i__] / h__[ii];
      ip = i__ + 1;
      ij = ii;
      i__2 = *n;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ii;
	  /* L410: */
	  v -= h__[ii] * d__[j];
	}
      /* L411: */
      d__[i__] = v;
    }
 L412:
  /*              calcul du pas minimum 
   *              et de la derivee directionnelle initiale 
   */
  c__ = 0.;
  dga = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*Computing MAX 
       */
      d__2 = c__, d__3 = (d__1 = d__[i__] / scale[i__], Abs (d__1));
      c__ = Max (d__2, d__3);
      /* L160: */
      dga += ga[i__] * d__[i__];
    }
  /*              test si la direction est de descente 
   */
  if (dga >= 0.)
    {
      goto L240;
    }
  /*              initialisation du pas 
   */
  stmin = 0.;
  stepbd = 0.;
  steplb = *acc / c__;
  fmin = fa;
  gmin = dga;
  step = 1.;
  if (dff <= 0.)
    {
      /*Computing MIN 
       */
      d__1 = step, d__2 = 1. / c__;
      step = Min (d__1, d__2);
    }
  if (dff > 0.)
    {
      /*Computing MIN 
       */
      d__1 = step, d__2 = (dff + dff) / (-dga);
      step = Min (d__1, d__2);
    }
  if (*iprint >= 2)
    {
      Sciprintf("iter num %d, nb calls=%d, f=%10.4f\n",itr,nfun,fa);
      if (*iprint >= 3)
	{
	  d__1 = dga / C2F(dnrm2) (n, &d__[1], &c__1);
	  Sciprintf("linear search: initial derivative=%10.4f\n",d__1);
	}
    }
  /*             boucle de reherche lineaire 
   */
 L170:
  c__ = stmin + step;
  if (nfun >= *nsim)
    {
      goto L250;
    }
  ++nfun;
  /*             calcul de fonction-gradient 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L180: */
      xb[i__] = xa[i__] + c__ * d__[i__];
    }
  indic = 4;
  (*simul) (&indic, n, &xb[1], &fb, &gb[1],optim_data);
  /*             test sur indic 
   */
  if (indic > 0)
    {
      goto L185;
    }
  if (indic < 0)
    {
      goto L183;
    }
  if (indic == 0) ret = FAIL;
  if (*iprint > 0)
    {
      Sciprintf("n1qn1 termine par voeu de l'utilisateur\n");
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      x[i__] = xb[i__];
      /* L182: */
      g[i__] = gb[i__];
    }
  goto L250;
 L183:
  stepbd = step;
  ial = 1;
  step /= 10.;
  if (*iprint >= 3)
    {
      Sciprintf("step length=%10.4f, indic=%d\n",c__,indic);
    }
  if (stepbd > steplb)
    {
      goto L170;
    }
  if (*iprint != 0 && isfv < 2)
    {
      Sciprintf("n1qn1 bute sur une contrainte implicite\n");
    }
  goto L240;
  /*            stockage si c'est la plus petite valeur 
   */
 L185:
  isfv = Min (2, isfv);
  if (fb > *f)
    {
      goto L220;
    }
  if (fb < *f)
    {
      goto L200;
    }
  gl1 = 0.;
  gl2 = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*Computing 2nd power 
       */
      d__1 = scale[i__] * g[i__];
      gl1 += d__1 * d__1;
      /* L190: */
      /*Computing 2nd power 
       */
      d__1 = scale[i__] * gb[i__];
      gl2 += d__1 * d__1;
    }
  if (gl2 >= gl1)
    {
      goto L220;
    }
 L200:
  isfv = 3;
  *f = fb;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      x[i__] = xb[i__];
      /* L210: */
      g[i__] = gb[i__];
    }
  /*              calcul de la derivee directionnelle 
   */
 L220:
  dgb = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L230: */
      dgb += gb[i__] * d__[i__];
    }
  if (*iprint < 3)
    {
      goto L231;
    }
  s = fb - fa;
  /*a small change (Bruno): to give a better indication about 
   * the directionnal derivative I scale it by || d || 
   */
  d__1 = dgb / C2F(dnrm2) (n, &d__[1], &c__1);
  Sciprintf("step length=%10.4g, df=%10.4g, derivative=%10.4g\n",c__,s,d__1);
  /*              test si la fonction a descendu 
   */
 L231:
  if (fb - fa <= c__ * .1 * dga)
    {
      goto L280;
    }
  ial = 0;
  /*              iteration terminee si le pas est minimum 
   */
  if (step > steplb)
    {
      goto L270;
    }
 L240:
  if (isfv >= 2)
    {
      goto L110;
    }
  /*              ici, tout est termine 
   */
 L250:
  if (*iprint > 0)
    {
      Sciprintf("iter num %d, nb calls=%d, f=%10.4f\n",itr,nfun,*f);
    }
  *acc = 0.;
  i__1 = *n; 
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L260: */
      *acc += g[i__] * g[i__];
    }
  *niter = itr;
  *nsim = nfun;
  return ret;
  /*              interpolation cubique 
   */
 L270:
  stepbd = step;
  c__ = gmin + dgb - (fb - fmin) * 3. / step;
  if (c__ == 0.)
    {
      goto L250;
    }
  cc = Abs (c__) - gmin * (dgb / Abs (c__));
  cc = sqrt ((Abs (c__))) * sqrt ((Max (0., cc)));
  c__ = (c__ - gmin + cc) / (dgb - gmin + cc + cc);
  step *= Max (.1, c__);
  goto L170;
  /*              ceci est un pas de descente 
   */
 L280:
  if (ial == 0)
    {
      goto L285;
    }
  if (stepbd > steplb)
    {
      goto L285;
    }
  if (*iprint != 0 && isfv < 2)
    { 
      Sciprintf("n1qn1: bute sur une contrainte implicite\n");
    }
  goto L240;
 L285:
  stepbd -= step;
  stmin = c__;
  fmin = fb;
  gmin = dgb;
  /*              extrapolation 
   */
  step = stmin * 9.;
  if (stepbd > 0.)
    {
      step = stepbd * .5;
    }
  c__ = dga + dgb * 3. - (fb - fa) * 4. / stmin;
  if (c__ > 0.)
    {
      /*Computing MIN 
       *Computing MAX 
       */
      d__3 = 1., d__4 = -dgb / c__;
      d__1 = step, d__2 = stmin * Max (d__3, d__4);
      step = Min (d__1, d__2);
    }
  if (dgb < dga * .7)
    {
      goto L170;
    }
  /*                recherche lineaire terminee, test de convergence 
   */
  isfv = 4 - isfv;
  if (stmin + step <= steplb)
    {
      goto L240;
    }
  /*                formule de bfgs 
   */
  ir = -(*n);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      xa[i__] = xb[i__];
      xb[i__] = ga[i__];
      d__[i__] = gb[i__] - ga[i__];
      /* L290: */
      ga[i__] = gb[i__];
    }
  d__1 = 1. / dga;
  optim_majour (&h__[1], &xb[1], &w[1], n, &d__1, &ir, &c__1, &c_b101);
  ir = -ir;
  d__1 = 1. / (stmin * (dgb - dga));
  optim_majour (&h__[1], &d__[1], &d__[1], n, &d__1, &ir, &c__1, &c_b101);
  /*                 test du rang de la nouvelle matrice 
   */
  if (ir < *n)
    {
      goto L250;
    }
  /*                 nouvelle iteration 
   */
  dff = fa - fb;
  fa = fb;
  goto L130;
} 

