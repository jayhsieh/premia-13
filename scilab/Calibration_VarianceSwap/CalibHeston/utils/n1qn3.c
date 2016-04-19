#include "optim.h"

static int optim_ddd2 (opt_prosca prosca,opt_ct ctonb,opt_ct ctcab, int *n, int *nm,
		       double *depl, double *aux, int *jmin, int *jmax, double *diag,
		       double *alpha, double *ybar, double *sbar, opt_simul_data *optim_data);
static int optim_n1qn3a (opt_simul simul, opt_prosca prosca,opt_ct ctonb,opt_ct ctcab, int *n,
			 double *x, double *f, double *g, double *dxmin, double *df1,
			 double *epsg, int *impres, int *io, int *mode, int *niter,
			 int *nsim, int *m, double *d__, double *gg, double *diag,
			 double *aux, double *alpha, double *ybar, double *sbar,
			 opt_simul_data *optim_data);

/*
 * 
 *    n1qn3, version 1.0, septembre 1988. 
 *    Jean Charles Gilbert, Claude Lemarechal, INRIA. 
 *    Copyright INRIA 
 *    email: Jean-Charles.Gilbert@inria.fr 
 * 
 *    Utilise les sous-routines: 
 *        n1qn3a 
 *        ddd2 
 *        nlis0 + dcube (XII/88) 
 * 
 *    La sous-routine n1qn3 est une interface entre le programme 
 *    appelant et la sous-routine n1qn3a, le minimiseur proprement dit. 
 * 
 *    Le module prosca est sense realiser le produit scalaire de deux 
 *    vecteurs de Rn; le module ctonb est sense realiser le changement 
 *    bases: base euclidienne -> base orthonormale (pour le produit 
 *    scalaire prosca); le module ctbas fait la transformation inverse: 
 *    base orthonormale -> base euclidienne. 
 * 
 *    dz est la zone de travail pour n1qn3a, de dimension ndz. 
 *    Elle est subdivisee en 
 *        4 vecteurs de dimension n: d,gg,diag,aux 
 *        m scalaires: alpha 
 *        m vecteurs de dimension n: ybar 
 *        m vecteurs de dimension n: sbar 
 * 
 *    m est alors le plus grand entier tel que 
 *        m*(2*n+1)+4*n .le. ndz, 
 *    soit m := (ndz-4*n) / (2*n+1) 
 *    Il faut avoir m >= 1, donc ndz >= 6n+1. 
 * 
 *    A chaque iteration la metrique est formee a partir d'une matrice 
 *    diagonale D qui est mise a jour m fois par la formule de BFGS en 
 *    utilisant les m couples {y,s} les plus recents. La matrice 
 *    diagonale est egale apres la premiere iteration a 
 *        (y,s)/|y|**2 * identite (facteur d'Oren-Spedicato) 
 *    et est elle-meme mise a jour a chaque iteration en utilisant la 
 *    formule de BFGS directe diagonalisee adaptee a l'ellipsoide de 
 *    Rayleigh. Si on note 
 *        D[i]:=(De[i],e[i]), y[i]:=(y,e[i]), s[i]:=(s,e[i]), 
 *    ou les e[i] forment une base orthonormale pour le produit scalaire 
 *    (.,.) que realise prosca, la formule de mise a jour de D s'ecrit: 
 *        D[i] := 1 / ( (Dy,y)/(y,s)/D[i] + y[i]**2/(y,s) 
 *                       - (Dy,y)*(s[i]/D[i])**2/(y,s)/(D**(-1)s,s) ) 
 * 
 */


int optim_n1qn3 (opt_simul simul,opt_prosca prosca, opt_ct ctonb, 
		 opt_ct ctcab, int *n,
		 double *x, double *f, double *g, double *dxmin, double *df1,
		 double *epsg, int *impres, int *io, int *mode, int *niter,
		 int *nsim, double *dz, int *ndz, opt_simul_data *optim_data)
{
  /* Local variables */
  int iaux, m, isbar;
  int iprec, iybar;
  double r1, r2;
  int l1memo, id, ialpha;
  double ps;
  int ntravu, igg;
  /* Parameter adjustments */
  --dz;
  --g;
  --x;
  if (*impres >= 1)
    {
      Sciprintf("n1qn3: entry point\n");
      Sciprintf("\tdimension of the problem n=%d\n",*n);
      Sciprintf("\tabsolute precision on x (dxmin)%9.2g\n",*dxmin);
      Sciprintf("\texpected decrease for f (df1):%9.2g\n",*df1);
      Sciprintf("\trelative precision on g (epsg):%9.2g\n",*epsg);
      Sciprintf("\tmaximal number of iterations (niter):%d\n",*niter);
      Sciprintf("\tmaximal number of simulations (nsim):%d\n",*nsim);
      Sciprintf("\tprinting level (impres):%d\n",*impres);
    }
  if (*n <= 0 || *niter <= 0 || *nsim <= 0 || *dxmin <= 0. || *epsg <= 0.
      || *epsg > 1.)
    {
      *mode = 2;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn1: inconsistent call\n");
	}
      return 0;
    }
  if (*ndz < *n * 6 + 1)
    {
      *mode = 2;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn1: not enough memory allocated\n");
	}
      goto L904;
    }
  /* 
   *---- calcul de m et des pointeurs subdivisant la zone de travail dz 
   * 
   */
  ntravu = *ndz - (*n << 2);
  l1memo = (*n << 1) + 1;
  m = ntravu / l1memo;
  ntravu = m * l1memo + (*n << 2);
  if (*impres >= 1)
    {
      Sciprintf("\tallocated memory (nrz) :%d\n",*ndz);
      Sciprintf("\tused memory: %d\n",ntravu);
      Sciprintf("\tnumber of updates: %d\n",m);
    }
  id = 1;
  igg = id + *n;
  iprec = igg + *n;
  iaux = iprec + *n;
  ialpha = iaux + *n;
  iybar = ialpha + m;
  isbar = iybar + *n * m;
  /* 
   *---- appel du code d"optimisation 
   * 
   */
  optim_n1qn3a (simul, prosca, ctonb, ctcab, n,
		&x[1], f, &g[1], dxmin, df1, epsg, impres, io, mode, niter,
		nsim, &m, &dz[id], &dz[igg], &dz[iprec], &dz[iaux],
		&dz[ialpha], &dz[iybar], &dz[isbar],optim_data);
  /* 
   *---- impressions finales 
   * 
   */
 L904:
  if (*impres >= 1)
    {
      Sciprintf("n1qn3: output mode is %d\n",*mode);
      Sciprintf("\tnumber of iterations: %d\n",*niter);
      Sciprintf("\tnumber of simulations: %d\n",*nsim);
      Sciprintf("\trealized relative precision on g: %9.2g\n",*epsg);
    }
  (*prosca) (n, &x[1], &x[1], &ps,optim_data);
  r1 = sqrt (ps);
  (*prosca) (n, &g[1], &g[1], &ps,optim_data);
  r2 = sqrt (ps);
  if (*impres >= 1)
    {
      Sciprintf("\tnorm of x = %15.8g\n",r1);
      Sciprintf("\tf         = %15.8g\n",*f);
      Sciprintf("\tnorm of g = %15.8g\n",r2);
    }
  return 0;
}


/*
 * 
 *    Code d'optimisation proprement dit. 
 */

static int optim_n1qn3a (opt_simul simul, opt_prosca prosca,opt_ct ctonb,opt_ct ctcab, int *n,
			 double *x, double *f, double *g, double *dxmin, double *df1,
			 double *epsg, int *impres, int *io, int *mode, int *niter,
			 int *nsim, int *m, double *d__, double *gg, double *diag,
			 double *aux, double *alpha, double *ybar, double *sbar,
			 opt_simul_data *optim_data)
{
  double c_b25 = .9;
  double c_b26 = 1e-4;
  /* System generated locals */
  int ybar_dim1, ybar_offset, sbar_dim1, sbar_offset, i__1;
  double d__1, d__2, d__3;
  /* Local variables */
  int jmin, jmax, isim, iter;
  double tmin, tmax;
  int i__;
  double t;
  int indic;
  double preco, gnorm, r1, ff, dk, ps, ys;
  int moderl;
  double precos, dk1, hp0, ps2, den;
  double eps1;

  /* Parameter adjustments */
  --aux;
  --diag;
  --gg;
  --d__;
  --g;
  --x;
  sbar_dim1 = *n;
  sbar_offset = sbar_dim1 + 1;
  sbar -= sbar_offset;
  ybar_dim1 = *n;
  ybar_offset = ybar_dim1 + 1;
  ybar -= ybar_offset;
  --alpha;

  /* Function Body */
  iter = 0;
  isim = 1;
  /* 
   */
  (*prosca) (n, &g[1], &g[1], &ps, optim_data);
  gnorm = sqrt (ps);
  if (*impres >= 1)
    {
      Sciprintf("\tf         = %15.8f\n",*f);
      Sciprintf("\tnorm of g = %15.8f\n",gnorm);
    }
  /* 
   *    ---- mise a l'echelle de la premiere direction de descente 
   * 
   *Computing 2nd power 
   */
  d__1 = gnorm;
  precos = *df1 * 2. / (d__1 * d__1);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      d__[i__] = -g[i__] * precos;
      /* L10: */
    }
  if (*impres >= 5)
    {
      Sciprintf("n1qn3a: descent direction -g: precon = %10.3g\n",precos);
    }
  if (*impres == 3)
    {
    }
  if (*impres == 4)
    {
    }
  /* 
   *    ---- initialisation pour nlis0 
   * 
   */
  tmax = 1e20;
  (*prosca) (n, &d__[1], &g[1], &hp0,optim_data);
  /* 
   *    ---- initialisation pour dd 
   * 
   */
  jmin = 1;
  jmax = 0;
  /* 
   *---- debut de l'iteration. On cherche x(k+1) de la forme x(k) + t*d, 
   *    avec t > 0. On connait d. 
   * 
   *        debut de la boucle: etiquette 100, 
   *        sortie de la boucle: goto 1000. 
   * 
   */
 L100:
  ++iter;
  if (*impres < 0)
    {
      if (iter % (-(*impres)) == 0)
	{
	  indic = 1;
	  (*simul) (&indic, n, &x[1], f, &g[1],optim_data);
	  goto L100;
	}
    }
  if (*impres >= 5)
    {
    }
  if (*impres >= 4)
    { 
    }
  if (*impres >= 3)
    {
      Sciprintf("n1qn3: iter %d, simul %d, f= %15.8g,  h'(0)= %12.5f\n",iter,isim,*f,hp0);
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      gg[i__] = g[i__];
      /* L101: */
    }
  ff = *f;
  /* 
   *    ---- recherche lineaire et nouveau point x(k+1) 
   * 
   */
  if (*impres >= 5)
    {
      Sciprintf("n1qn3: line search\n");
    }
  /* 
   *        ---- calcul de tmin 
   * 
   */
  tmin = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*Computing MAX 
       */
      d__2 = tmin, d__3 = (d__1 = d__[i__], Abs (d__1));
      tmin = Max (d__2, d__3);
      /* L200: */
    }
  tmin = *dxmin / tmin;
  t = 1.;
  r1 = hp0;
  /* 
   */
  optim_nlis0 (n,  simul,  prosca, &x[1], f, &r1, &t, &tmin,
	       &tmax, &d__[1], &g[1], &c_b25, &c_b26, impres, io, &moderl,
	       &isim, nsim, &aux[1],optim_data);
  /* 
   *         ---- nlis0 renvoie les nouvelles valeurs de x, f et g 
   * 
   */
  if (moderl != 0)
    {
      if (moderl < 0)
	{
	  /* 
	   *            ---- calcul impossible 
	   *                 t, g: ou les calculs sont impossible 
	   *                 x, f: ceux du t_gauche (donc f <= ff) 
	   * 
	   */
	  *mode = moderl;
	}
      else if (moderl == 1)
	{
	  /* 
	   *            ---- descente bloquee sur tmax 
	   *                 [sortie rare (!!) d'apres le code de nlis0] 
	   * 
	   */
	  *mode = 3;
	  if (*impres >= 1)
	    {
	      Sciprintf("n1qn3: iteration %d, line search blocked on tmax: decrease the scaling\n",iter);
 	    }
	}
      else if (moderl == 4)
	{
	  /* 
	   *            ---- nsim atteint 
	   *                 x, f: ceux du t_gauche (donc f <= ff) 
	   * 
	   */
	  *mode = 5;
	}
      else if (moderl == 5)
	{
	  /* 
	   *            ---- arret demande par l'utilisateur (indic = 0) 
	   *                 x, f: ceux en sortie du simulateur 
	   * 
	   */
	  *mode = 0;
	}
      else if (moderl == 6)
	{
	  /* 
	   *            ---- arret sur dxmin ou appel incoherent 
	   *                 x, f: ceux du t_gauche (donc f <= ff) 
	   * 
	   */
	  *mode = 6;
	}
      goto L1000;
    }
  /* 
   *    ---- tests d'arret 
   * 
   */
  (*prosca) (n, &g[1], &g[1], &ps,optim_data);
  eps1 = sqrt (ps) / gnorm;
  /* 
   */
  if (*impres >= 5)
    {
      Sciprintf("n1qn3: stopping criterion on g: %12.5g\n",eps1);
    }
  if (eps1 < *epsg)
    {
      *mode = 1;
      goto L1000;
    }
  if (iter == *niter)
    {
      *mode = 4;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn3: iteration %d, maximal number of iterations\n",iter);
 	}
      goto L1000;
    }
  if (isim > *nsim)
    {
      *mode = 5;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn3: iteration %d, %d simulations (maximal number reached)\n",iter,isim);
	}
      goto L1000;
    }
  /* 
   *    ---- mise a jour de la matrice 
   * 
   */
  if (*m > 0)
    {
      ++jmax;
      if (iter > *m)
	{
	  ++jmin;
	  if (jmin > *m)
	    {
	      jmin -= *m;
	    }
	  if (jmax > *m)
	    {
	      jmax -= *m;
	    }
	}
      /* 
       *         ---- y, s et (y,s) 
       * 
       */
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  sbar[i__ + jmax * sbar_dim1] = t * d__[i__];
	  ybar[i__ + jmax * ybar_dim1] = g[i__] - gg[i__];
	  /* L400: */
	}
      if (*impres >= 5)
	{
	  (*prosca) (n, &sbar[jmax * sbar_dim1 + 1],
		     &sbar[jmax * sbar_dim1 + 1], &ps,optim_data);
	  dk1 = sqrt (ps);
	  if (iter > 1)
	    {
	      d__1 = dk1 / dk;
	      Sciprintf("n1qn3: convergence rate, s(k)/s(k-1) = %12.5f\n",d__1);
	    }
	  dk = dk1;
	}
      (*prosca) (n, &ybar[jmax * ybar_dim1 + 1], &sbar[jmax * sbar_dim1 + 1],
		 &ps,optim_data);
      ys = ps;
      if (ys <= 0.)
	{
	  *mode = 7;
	  if (*impres >= 1)
	    {
	      Sciprintf("n1qn3: iteration %d, the scalar product (y,s) = %f12.5 is not positive\n",
			iter,ys);
	    }
	  goto L1000;
	}
      if (*impres >= 5)
	{
	  Sciprintf("n1qn3: matrix update:\n");
	}
      /* 
       *         ---- ybar et sbar 
       * 
       */
      r1 = sqrt (1. / ys);
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  sbar[i__ + jmax * sbar_dim1] = r1 * sbar[i__ + jmax * sbar_dim1];
	  ybar[i__ + jmax * ybar_dim1] = r1 * ybar[i__ + jmax * ybar_dim1];
	  /* L410: */
	}
      /* 
       *         ---- calcul de la diagonale de preconditionnement 
       * 
       */
      (*prosca) (n, &ybar[jmax * ybar_dim1 + 1], &ybar[jmax * ybar_dim1 + 1],
		 &ps,optim_data);
      precos = 1. / ps;
      if (iter == 1)
	{
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      diag[i__] = precos;
	      /* L401: */
	    }
	}
      else
	{
	  /* 
	   *            ---- ajustememt de la diagonale a l'ellipsoide de Rayleigh 
	   * 
	   */
	  (*ctonb) (n, &ybar[jmax * ybar_dim1 + 1], &aux[1],optim_data);
	  r1 = 0.;
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      /*Computing 2nd power 
	       */
	      d__1 = aux[i__];
	      r1 += diag[i__] * (d__1 * d__1);
	      /* L398: */
	    }
	  if (*impres >= 5)
	    {
	      d__1 = 1. / r1;
	      Sciprintf("fitting the ellipsoid: factor %10.3f\n",d__1);
 	    }
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      diag[i__] /= r1;
	      /* L399: */
	    }
	  /* 
	   *            ---- mise a jour diagonale 
	   *                 gg utilise comme vecteur auxiliaire 
	   * 
	   */
	  (*ctonb) (n, &sbar[jmax * sbar_dim1 + 1], &gg[1],optim_data);
	  den = 0.;
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      /*Computing 2nd power 
	       */
	      d__1 = gg[i__];
	      den += d__1 * d__1 / diag[i__];
	      /* L402: */
	    }
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      /*Computing 2nd power 
	       */
	      d__1 = aux[i__];
	      /*Computing 2nd power 
	       */
	      d__2 = gg[i__] / diag[i__];
	      diag[i__] =
		1. / (1. / diag[i__] + d__1 * d__1 - d__2 * d__2 / den);
	      /* L403: */
	    }
	}
      if (*impres >= 5)
	{
	  preco = 0.;
	  i__1 = *n;
	  for (i__ = 1; i__ <= i__1; ++i__)
	    {
	      preco += diag[i__];
	      /* L406: */
	    }
	  preco /= *n;
	  Sciprintf("Oren-Spedicato factor (not used) = %10.3f\n",precos);
	  Sciprintf("diagonal: average value = %10.3f\n",preco);
	}
    }
  /* 
   *    ---- calcul de la nouvelle direction de descente d = - h.g 
   * 
   */
  if (*m == 0)
    {
      /*Computing 2nd power 
       */
      d__1 = eps1 * gnorm;
      preco = (ff - *f) * 2. / (d__1 * d__1);
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  d__[i__] = -g[i__] * preco;
	  /* L500: */
	}
    }
  else
    {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  d__[i__] = -g[i__];
	  /* L510: */
	}
      optim_ddd2 ( prosca,ctonb, ctcab, n, m , &d__[1],
		   &aux[1], &jmin, &jmax, &diag[1], &alpha[1],
		   &ybar[ybar_offset], &sbar[sbar_offset],optim_data);
    }
  /* 
   *         ---- test: la direction d est-elle de descente ? 
   *              hp0 sera utilise par nlis0 
   * 
   */
  (*prosca) (n, &d__[1], &g[1], &hp0,optim_data);
  if (hp0 >= 0.)
    {
      *mode = 7;
      if (*impres >= 1)
	{
	  Sciprintf("n1qn3 iteration %d: the search direction d is not a descent direction: (g,d) = %12.5f",iter,hp0);
	}
      goto L1000;
    }
  if (*impres >= 5)
    {
      (*prosca) (n, &g[1], &g[1], &ps,optim_data);
      ps = sqrt (ps);
      (*prosca) (n, &d__[1], &d__[1], &ps2,optim_data);
      ps2 = sqrt (ps2);
      ps = hp0 / ps / ps2;
      /*Computing MIN 
       */
      d__1 = -ps;
      ps = Min (d__1, 1.);
      ps = acos (ps);
      r1 = ps * 180. / 3.1415927;
      Sciprintf("n1qn3: descent direction d: angle(-g,d) = %f5.1, degrees\n",r1);
    }
  /* 
   *---- on poursuit les iterations 
   * 
   */
  goto L100;
  /* 
   *---- retour 
   * 
   */
 L1000:
  *niter = iter;
  *nsim = isim;
  *epsg = eps1;
  return 0;
}




/*    Copyright INRIA 
 *---- 
 * 
 *    calcule le produit h g ou 
 *        . h est une matrice construite par la formule de bfgs inverse 
 *          a nm memoires a partir de la matrice diagonale diag 
 *          dans un espace hilbertien dont le produit scalaire 
 *          est donne par prosca 
 *          (cf. J. Nocedal, Math. of Comp. 35/151 (1980) 773-782) 
 *        . g est un vecteur de dimension n (en general le gradient) 
 * 
 *    la matrice diag apparait donc comme un preconditionneur diagonal 
 * 
 *    depl = g (en entree), = h g (en sortie) 
 * 
 *    la matrice h est memorisee par les vecteurs des tableaux 
 *    ybar, sbar et les pointeurs jmin, jmax 
 * 
 *    alpha(nm) est une zone de travail 
 * 
 *    izs(1),rzs(1),dzs(1) sont des zones de travail pour prosca 
 * 
 */

static int optim_ddd2 (opt_prosca prosca,opt_ct ctonb,opt_ct ctcab, int *n, int *nm,
		       double *depl, double *aux, int *jmin, int *jmax, double *diag,
		       double *alpha, double *ybar, double *sbar, opt_simul_data *optim_data)
{
  int ybar_dim1, ybar_offset, sbar_dim1, sbar_offset, i__1, i__2;
  int jfin, i__, j;
  double r__;
  int jp;
  double ps;

  /* Parameter adjustments */
  --diag;
  --aux;
  --depl;
  sbar_dim1 = *n;
  sbar_offset = sbar_dim1 + 1;
  sbar -= sbar_offset;
  ybar_dim1 = *n;
  ybar_offset = ybar_dim1 + 1;
  ybar -= ybar_offset;
  --alpha;


  /* Function Body */
  jfin = *jmax;
  if (jfin < *jmin)
    {
      jfin = *jmax + *nm;
    }
  /* 
   *        phase de descente 
   * 
   */
  i__1 = *jmin;
  for (j = jfin; j >= i__1; --j)
    {
      jp = j;
      if (jp > *nm)
	{
	  jp -= *nm;
	}
      (*prosca) (n, &depl[1], &sbar[jp * sbar_dim1 + 1], &ps,optim_data);
      r__ = ps;
      alpha[jp] = r__;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  depl[i__] -= r__ * ybar[i__ + jp * ybar_dim1];
	  /* L20: */
	}
      /* L100: */
    }
  /* 
   *        preconditionnement 
   * 
   */
  (*ctonb) (n, &depl[1], &aux[1], optim_data);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      aux[i__] *= diag[i__];
      /* L150: */
    }
  (*ctcab) (n, &aux[1], &depl[1], optim_data);
  /* 
   *        remontee 
   * 
   */
  i__1 = jfin;
  for (j = *jmin; j <= i__1; ++j)
    {
      jp = j;
      if (jp > *nm)
	{
	  jp -= *nm;
	}
      (*prosca) (n, &depl[1], &ybar[jp * ybar_dim1 + 1], &ps,optim_data);
      r__ = alpha[jp] - ps;
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  depl[i__] += r__ * sbar[i__ + jp * sbar_dim1];
	  /* L120: */
	}
      /* L200: */
    }
  return 0;
}



