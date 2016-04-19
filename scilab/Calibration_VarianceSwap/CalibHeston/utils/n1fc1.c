#include "optim.h"

static int optim_fprf2 (int *, int *, int *, int *, double *, double *, double *,    
			double *, int *, double *, double *, int *, int *, int *,    double *, double *, 
			double *, double *, double *, double *,    double *, double *); 

static int optim_n1fc1o (int *unit, int *job, int *i1, int *i2, int *i3, int *i4,
			 int *i5, double *d1, double *d2, double *d3, double *d4);


static int optim_frdf1 (opt_prosca prosca, int *n, int *ntot, int *ninf, int *kgrad,
			double *al, double *q, double *s, double *poids, double *aps,
			double *anc, int *mm1, double *r__, double *e, int *ic, opt_simul_data *optim_data);

static int optim_fremf2 (opt_prosca prosca, int *iflag, int *n, int *ntot, int *nta, int *mm1,
			 double *p, double *alfa, double *e, double *a, double *r__,opt_simul_data *optim_data);


static int optim_nlis2 (opt_simul simul, opt_prosca prosca, int *n, double *xn, double *fn,
			double *fpn, double *t, double *tmin, double *tmax, double *d__,
			double *d2, double *g, double *gd, double *amd, double *amf,
			int *imp, int *io, int *logic, int *nap, int *napmax, double *x,
			double *tol, double *a, double *tps, double *tnc, double *gg,
			opt_simul_data *optim_data);
/*
  static int optim_fremf1 (opt_prosca prosca, int *iflag, int *n, int *ntot, int *nta, int *mm1,
  double *p, double *alfa, double *e, double *a, double *r__,
  opt_simul_data *optim_data);
*/
static int optim_n1fc1a (opt_simul simul, opt_prosca prosca, int *n, int *mode, double *xn,
			 double *fn, double *g, double *df0, double *eps0, double *dx,
			 int *imp, double *zero, int *io, int *ntot, int *iter,
			 int *nsim, int *memax, double *s, double *gd, double *x,
			 double *sa, double *gg, double *al, double *aps, double *anc,
			 double *poids, double *q, int *jc, int *ic, double *r__,
			 double *a, double *e, double *rr, double *xga, double *y,
			 double *w1, double *w2, opt_simul_data *optim_data);

static int optim_ffinf1 (int *n, int *nv, int *jc, double *xpr, double *p, double *s);

/*
 * Copyright INRIA 
 *
 * minimisation d'une fonction hemiderivable par une methode de faisceau. 
 * la direction est obtenue par projection de l'origine 
 * sur un polyedre genere par un ensemble de gradients deja calcules 
 * et la recherche lineaire donne un pas de descente ou un pas nul. 
 * l'algorithme minimise f a eps0 pres (si convexite) 
 * ou eps0 est une tolerance donnee par l'utilisateur. 
 * 
 * mode 
 *               <=0 mode=indic de simul 
 *               1 fin normale 
 *               2 appel incoherent 
 *               3 reduire l'echelle des x 
 *               4 max iterations 
 *               5 max simulations 
 *               6 impossible d'aller au dela de dx 
 *               7 fprf2 mis en echec 
 *               8 on commence a boucler 
 * imp 
 *               <0 indic=1 toutes les -imp iterations 
 *               0 pas d'impressions 
 *               1 impressions initiales et finales 
 *               2 impressions a chaque convergence 
 *               3 une impression par iteration 
 *               4 informations n1fc1 et nlis2 
 *               >4 debug 
 *                        5 tolerances diverses 
 *                        6 poids 
 *                        >6 fprf2 
 * -------------------------------------------------- 
 */

int optim_n1fc1 (opt_simul simul, opt_prosca prosca, int *n, double *xn, double *fn,
		 double *g, double *dxmin, double *df1, double *epsf,
		 double *zero, int *imp, int *io, int *mode, int *iter, int *nsim,
		 int *memax, int *iz, double *rz, double *dz,opt_simul_data *optim_data)
{
  int c__1 = 1;
  int c__2 = 2;
  int i__1;
  int nanc, nxga, naps, ntot, i__;
  double d1, d2, d3[1], d4[1];
  int i1, i2, i3, i4, i5[1], na, ne, nq, nr, ns, nx, ny, npoids, nw1, nw2,
    ngd, nic, nal, ngg, njc, nsa, ndz, nrr, niz, nrz;

  /*         dimension iz=2*(memax+1) 
   *         dimension rz=5*n+(n+4)*memax 
   *         dimension dz=(memax+9)*memax+8 
   *    Copyright INRIA 
   * 
   */
  /* Parameter adjustments */
  --g;
  --xn;
  --iz;
  --rz;
  --dz;

  /* Function Body */
  if (*n > 0 && *df1 > 0. && *epsf >= 0. && *zero >= 0. && *iter >= 0
      && *nsim >= 0 && *memax >= 1 && *dxmin > 0.)
    {
      goto L10;
    }
  *mode = 2;
  /*    appel incoherent 
   */
  optim_n1fc1o (io, &c__1, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
  goto L999;
 L10:
  ns = 1;
  ngd = ns + *n;
  nx = ngd + *n;
  nsa = nx + *n;
  ngg = nsa + *n;
  nal = ngg + *n;
  naps = nal + *memax;
  nanc = naps + *memax;
  npoids = nanc + *memax;
  nq = npoids + *memax;
  njc = 1;
  nic = njc + *memax + 1;
  nr = 1;
  na = nr + (*memax + 1) * (*memax + 1);
  ne = na + *memax + 1;
  nrr = ne + *memax + 1;
  nxga = nrr + *memax + 1;
  ny = nxga + *memax + 1;
  nw1 = ny + *memax + 1;
  nw2 = nw1 + *memax + 1;
  /* 
   */
  niz = (*memax + 1) << 1;
  nrz = nq + *n * *memax - 1;
  ndz = nw2 + *memax;
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__2, n, memax, &niz, &nrz, &ndz, &d1, &d2, d3, d4);
    }
  i__1 = niz;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L110: */
      iz[i__] = 0;
    }
  i__1 = nrz;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L120: */
      rz[i__] = 0.;
    }
  i__1 = ndz;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L130: */
      dz[i__] = 0.;
    }
  optim_n1fc1a ( simul, prosca, n, mode, &xn[1], fn, &g[1], df1,
		 epsf, dxmin, imp, zero, io, &ntot, iter, nsim, memax, &rz[ns],
		 &rz[ngd], &rz[nx], &rz[nsa], &rz[ngg], &rz[nal], &rz[naps],
		 &rz[nanc], &rz[npoids], &rz[nq], &iz[njc], &iz[nic], &dz[nr],
		 &dz[na], &dz[ne], &dz[nrr], &dz[nxga], &dz[ny], &dz[nw1],
		 &dz[nw2], optim_data);
  iz[1] = ntot;
 L999:
  return 0;
}				


static int optim_n1fc1a (opt_simul simul, opt_prosca prosca, int *n, int *mode, double *xn,
			 double *fn, double *g, double *df0, double *eps0, double *dx,
			 int *imp, double *zero, int *io, int *ntot, int *iter,
			 int *nsim, int *memax, double *s, double *gd, double *x,
			 double *sa, double *gg, double *al, double *aps, double *anc,
			 double *poids, double *q, int *jc, int *ic, double *r__,
			 double *a, double *e, double *rr, double *xga, double *y,
			 double *w1, double *w2, opt_simul_data *optim_data)
{
  int c__3 = 3;
  int c__4 = 4;
  int c__5 = 5;
  int c__6 = 6;
  int c__7 = 7;
  int c__8 = 8;
  int c__9 = 9;
  int c__10 = 10;
  int c__11 = 11;
  int c__12 = 12;
  int c__13 = 13;
  int c__14 = 14;
  int c__15 = 15;
  int c__16 = 16;
  int c__17 = 17;
  int c__18 = 18;
  int c__19 = 19;
  int c__20 = 20;
  /* System generated locals */
  int i__1;
  double d__1, d__2, d__3;
  /* Local variables */
  double alfa, beta, epsm, tmin, tmax, diam2;
  double f;
  int i__, k, j, iflag;
  double u;
  int indic;
  int kgrad;
  double z__;
  int logic, itmax, itimp;
  double d1, d2, d3[1], d4[1];
  int i1, i2, i3, i4, i5[1];
  double ajust, s2, s3;
  double z1, z2;
  int logic2;
  int memax1;
  double fa, df, ap;
  int nk, mm;
  double ro, ps;
  int nv, napmax, nt1;
  double s3n, roa;
  int nta;
  double fpn, tnc, eps;
  int nki;
  double tol, tps, eta2;

  /* Parameter adjustments */
  --gg;
  --sa;
  --x;
  --gd;
  --s;
  --g;
  --xn;
  --poids;
  --anc;
  --aps;
  --al;
  --q;
  --jc;
  --ic;
  --r__;
  --a;
  --e;
  --rr;
  --xga;
  --y;
  --w1;
  --w2;

  /* Function Body */
  itmax = *iter;
  *iter = 0;
  itimp = 0;
  napmax = *nsim;
  *nsim = 1;
  logic = 1;
  logic2 = 0;
  tmax = 1e20;
  eps = *df0;
  epsm = eps;
  df = *df0;
  *mode = 1;
  *ntot = 0;
  iflag = 0;
  /* 
   *         initialisation du faisceau 
   *         calcul du diametre de l'epure et du test d'arret 
   * 
   */
  aps[1] = 0.;
  anc[1] = 0.;
  poids[1] = 0.;
  nta = 0;
  kgrad = 1;
  memax1 = *memax + 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L50: */
      q[i__] = -g[i__];
    }
  (*prosca) (n, &g[1], &g[1], &ps, optim_data);
  if (ps > 0.)
    {
      goto L60;
    }
  *mode = 2;
  if (*imp != 0)
    {
      optim_n1fc1o (io, &c__3, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  goto L900;
 L60:
  diam2 = *df0 * 100. * *df0 / ps;
  eta2 = *eps0 * .01 * *eps0 / diam2;
  ap = *zero * *df0 / diam2;
  if (*imp > 2)
    {
      optim_n1fc1o (io, &c__4, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  /* 
   *             boucle 
   * 
   */
 L100:
  ++(*iter);
  ++itimp;
  if (*iter < itmax)
    {
      goto L110;
    }
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__5, iter, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  *mode = 4;
  goto L900;
 L110:
  ++(*ntot);
  if (logic == 3)
    {
      ro *= sqrt (s2);
    }
  if (itimp != -(*imp))
    {
      goto L200;
    }
  itimp = 0;
  indic = 1;
  (*simul) (&indic, n, &xn[1], &f, &g[1],optim_data);
  /* 
   *        calcul de la direction 
   * 
   */
 L200:
  eps = Min (eps, epsm);
  eps = Max (eps, *eps0);
  optim_fremf2 ( prosca, &iflag, n, ntot, &nta, &memax1, &q[1],
		 &poids[1], &e[1], &a[1], &r__[1], optim_data);
  optim_fprf2 (&iflag, ntot, &nv, io, zero, &s2, &eps, &al[1], imp, &u, &eta2,
	       &memax1, &jc[1], &ic[1], &r__[1], &a[1], &e[1], &rr[1],
	       &xga[1], &y[1], &w1[1], &w2[1]);
  /* 
   *        fin anormale de fprf2 
   * 
   */
  if (iflag == 0)
    {
      goto L250;
    }
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__6, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  *mode = 7;
  goto L900;
 L250:
  nta = *ntot;
  optim_ffinf1 (n, &nv, &jc[1], &xga[1], &q[1], &s[1]);
  u = Max (u, 0.);
  s2 = Max (s2, 0.);
  /* 
   *         tests d'arret (nb. si nr g est interieur a g 
   *                               alors nr g est "nul") 
   * 
   */
  if (nv < *n + 2)
    {
      goto L260;
    }
  /*Computing MAX 
   */
  d__1 = eta2, d__2 = s2 * 10.;
  eta2 = Max (d__1, d__2);
  if (*imp >= 2)
    {
      optim_n1fc1o (io, &c__7, &i1, &i2, &i3, &i4, i5, &eta2, &d2, d3, d4);
    }
 L260:
  if (s2 > eta2)
    {
      goto L300;
    }
  /* 
   *        calcul de la precision 
   */
  z__ = 0.;
  i__1 = nv;
  for (k = 1; k <= i__1; ++k)
    {
      j = jc[k] - 1;
      if (j > 0)
	{
	  z__ += xga[k] * poids[j];
	}
      /* L270: */
    }
  epsm = Min (eps, z__);
  if (*imp >= 2)
    {
      optim_n1fc1o (io, &c__8, iter, nsim, &i3, &i4, i5, fn, &epsm, &s2, d4);
    }
  if (epsm > *eps0)
    {
      goto L280;
    }
  *mode = 1;
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__9, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  goto L900;
  /* 
   *        diminution de epsilon 
   */
 L280:
  /*Computing MAX 
   */
  d__1 = epsm * .1;
  epsm = Max (d__1, *eps0);
  eps = epsm;
  if (logic == 3)
    {
      tol = eps * .01;
    }
  iflag = 2;
  goto L200;
  /* 
   *                suite des iterations 
   *                   impressions 
   * 
   */
 L300:
  if (*imp > 3)
    {
      optim_n1fc1o (io, &c__10, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  if (*imp > 2)
    {
      optim_n1fc1o (io, &c__11, iter, nsim, &nv, &i4, i5, fn, &eps, &s2, &u);
    }
  if (*imp >= 6)
    {
      optim_n1fc1o (io, &c__12, ntot, &i2, &i3, &i4, i5, &d1, &d2, d3,
		    &poids[1]);
    }
  /*               test de non-pivotage 
   */
  if (logic != 3)
    {
      goto L350;
    }
  z__ = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      z1 = s[i__] - sa[i__];
      /* L310: */
      z__ += z1 * z1;
    }
  if (z__ > *zero * 10. * *zero * s2)
    {
      goto L350;
    }
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__13, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  *mode = 8;
  goto L900;
  /* 
   *               recherche lineaire 
   * 
   */
 L350:
  iflag = 3;
  s3 = s2 + u * eps;
  if (logic == 3)
    {
      goto L365;
    }
  ro = df * 2. / s3;
  tol = eps * .01;
  goto L370;
 L365:
  ro /= sqrt (s2);
  /*Computing MAX 
   */
  d__1 = tol * .6, d__2 = *eps0 * .01;
  tol = Max (d__1, d__2);
 L370:
  fa = *fn;
  alfa = .2;
  beta = .1;
  fpn = -s3;
  if (*memax == 1)
    {
      tol = 0.;
    }
  /*                calcul de la resolution minimale, fonction de dx 
   */
  tmin = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L372: */
      /*Computing MAX 
       */
      d__2 = tmin, d__3 = (d__1 = s[i__] / *dx, Abs (d__1));
      tmin = Max (d__2, d__3);
    }
  tmin = 1. / tmin;
  if (*iter == 1)
    {
      roa = ro;
    }
  optim_nlis2 ( simul,  prosca, n, &xn[1], fn, &fpn, &ro, &tmin,
		&tmax, &s[1], &s2, &g[1], &gd[1], &alfa, &beta, imp, io,
		&logic, nsim, &napmax, &x[1], &tol, &ap, &tps, &tnc, &gg[1],optim_data);
  if (logic == 0 || logic == 2 || logic == 3)
    {
      goto L380;
    }
  /*                sortie par anomalie dans nlis2 
   */
  if (*imp <= 0)
    {
      goto L375;
    }
  if (logic == 6 || logic < 0)
    {
      optim_n1fc1o (io, &c__14, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  if (logic == 4)
    {
      optim_n1fc1o (io, &c__15, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  if (logic == 5)
    {
      optim_n1fc1o (io, &c__16, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  if (logic == 1)
    {
      optim_n1fc1o (io, &c__17, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
 L375:
  if (logic == 1)
    {
      *mode = 3;
    }
  if (logic == 4)
    {
      *mode = 5;
    }
  if (logic == 5)
    {
      *mode = 0;
    }
  if (logic == 6)
    {
      *mode = 6;
    }
  if (logic < 0)
    {
      *mode = logic;
    }
  goto L900;
 L380:
  if (logic != 3)
    {
      goto L385;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L382: */
      sa[i__] = s[i__];
    }
 L385:
  if (*iter > 1)
    {
      goto L390;
    }
  /* 
   *             1ere iteration, ajustement de ap, diam et eta 
   */
  if (logic == 0)
    {
      tps = *fn - fa - ro * fpn;
    }
  ap = *zero * *zero * Abs (tps) / (s2 * ro * ro);
  ajust = ro / roa;
  if (logic != 3)
    {
      diam2 = diam2 * ajust * ajust;
    }
  if (logic != 3)
    {
      eta2 /= ajust * ajust;
    }
  if (*imp >= 2)
    {
      optim_n1fc1o (io, &c__18, &i1, &i2, &i3, &i4, i5, &diam2, &eta2, &ap,
		    d4);
    }
 L390:
  mm = *memax - 1;
  if (logic == 2)
    {
      mm = *memax - 2;
    }
  if (*ntot <= mm)
    {
      goto L400;
    }
  /* 
   *     reduction du faisceau pour entrer le nouvel element 
   * 
   */
  optim_frdf1 ( prosca, n, ntot, &mm, &kgrad, &al[1], &q[1], &s[1],
		&poids[1], &aps[1], &anc[1], &memax1, &r__[1], &e[1], &ic[1],
		optim_data);
  iflag = 1;
  nta = *ntot;
  if (*imp >= 2)
    {
      optim_n1fc1o (io, &c__19, iter, nsim, ntot, &i4, i5, fn, &d2, d3, d4);
    }
  /* 
   */
 L400:
  if (*imp >= 5)
    {
      optim_n1fc1o (io, &c__20, &logic, &i2, &i3, &i4, i5, &ro, &tps, &tnc,
		    d4);
    }
  if (logic == 3)
    {
      goto L500;
    }
  /* 
   *                iteration de descente 
   * 
   */
  iflag = Min (iflag, 2);
  df = fa - *fn;
  if (*ntot == 0)
    {
      goto L500;
    }
  /* 
   *              actualisation des poids 
   * 
   */
  s3n = ro * sqrt (s2);
  i__1 = *ntot;
  for (k = 1; k <= i__1; ++k)
    {
      nk = (k - 1) * *n + 1;
      (*prosca) (n, &q[nk], &s[1], &ps,optim_data);
      z1 = (d__1 = aps[k] + (-df + ro * ps), Abs (d__1));
      z2 = anc[k] + s3n;
      /*Computing MAX 
       */
      d__1 = z1, d__2 = ap * z2 * z2;
      poids[k] = Max (d__1, d__2);
      aps[k] = z1;
      anc[k] = z2;
      /* L430: */
    }
  /* 
   *               actualisation de eps 
   * 
   */
  eps = ro * s3;
  kgrad = *ntot + 1;
  /* 
   *      nouvel element du faisceau (pour les trois types de pas) 
   * 
   */
 L500:
  nt1 = *ntot + 1;
  if (logic == 3)
    {
      goto L510;
    }
  aps[nt1] = 0.;
  anc[nt1] = 0.;
  poids[nt1] = 0.;
  goto L520;
 L510:
  aps[nt1] = tps;
  anc[nt1] = sqrt (tnc);
  /*Computing MAX 
   */
  d__1 = tps, d__2 = ap * tnc;
  poids[nt1] = Max (d__1, d__2);
 L520:
  nk = *ntot * *n;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      nki = nk + i__;
      /* L530: */
      q[nki] = -g[i__];
    }
  /* 
   *     traitement pour logic=2 (on ajoute encore un gradient) 
   */
  if (logic != 2)
    {
      goto L550;
    }
  ++(*ntot);
  logic = 3;
  logic2 = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L540: */
      g[i__] = gd[i__];
    }
  goto L390;
 L550:
  logic -= logic2;
  logic2 = 0;
  goto L100;
  /* 
   *               epilogue 
   * 
   */
 L900:
  if (*iter <= 1)
    {
      goto L990;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L910: */
      g[i__] = -s[i__];
    }
 L990:
  return 0;
}				/* n1fc1a_ */

/* 
 *         cette subroutine calcule s = sigma xpr(.)*p(.) 
 *         sachant que les xpr ont ete calcules par fprf2 
 * 
 */

static int optim_ffinf1 (int *n, int *nv, int *jc, double *xpr, double *p, double *s)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int i__, j, k;
  double ps;
  int nij;

  /* Parameter adjustments */
  --s;
  --xpr;
  --jc;
  --p;

  /* Function Body */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ps = 0.;
      i__2 = *nv;
      for (k = 1; k <= i__2; ++k)
	{
	  j = jc[k] - 1;
	  if (j == 0)
	    {
	      goto L910;
	    }
	  nij = (j - 1) * *n + i__;
	  ps += xpr[k] * p[nij];
	L910:
	  ;
	}
      /* L920: */
      s[i__] = ps;
    }
  return 0;
}				/* ffinf1_ */




static int optim_fremf2 (opt_prosca prosca, int *iflag, int *n, int *ntot, int *nta, int *mm1,
			 double *p, double *alfa, double *e, double *a, double *r__,opt_simul_data *optim_data)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int mekk, i__, j, jj, kk, ni, nj;
  double ps;
  int nt1, mej, nij, nta1, nta2;

  /*    Copyright INRIA 
   * 
   *         cette subroutine remplit les donnees pour fprf2 
   *         (produits scalaires et 2 contraintes lineaires) 
   * 
   *            de 1 a ntot +1  si iflag=0 
   *            de nta+1 +1 a ntot +1 sinon 
   * 
   *            (le +1 est du a l'ecart, place en premier) 
   * 
   *         p contient les opposes des gradients a la queue leu leu 
   *         -g(1), -g(2),..., -g(ntot) soit ntot*n coordonnees 
   * 
   */
  /* Parameter adjustments */
  --alfa;
  --a;
  --e;
  --p;
  --r__;

  /* Function Body */
  nt1 = *ntot + 1;
  nta1 = *nta + 1;
  if (*iflag > 0)
    {
      goto L50;
    }
  /* 
   *               remplissage des anciennes donnees 
   *         (produits scalaires, ecart et contrainte d'egalite) 
   * 
   */
  i__1 = *ntot;
  for (j = 1; j <= i__1; ++j)
    {
      jj = (j - 1) * *mm1 + 1;
      /* L10: */
      r__[jj] = 0.;
    }
  a[1] = 1.;
  e[1] = 0.;
  if (nta1 == 1)
    {
      goto L50;
    }
  i__1 = nta1;
  for (j = 2; j <= i__1; ++j)
    {
      e[j] = 1.;
      nj = (j - 2) * *n;
      mej = (j - 1) * *mm1;
      i__2 = j;
      for (i__ = 2; i__ <= i__2; ++i__)
	{
	  ni = (i__ - 2) * *n;
	  /* 
	   *            produit scalaire de g(i-1) avec g(j-1) 
	   *            pour j-1=1,nta et i-1=1,j-1 
	   * 
	   */
	  (*prosca) (n, &p[ni + 1], &p[nj + 1], &ps,optim_data);
	  nij = mej + i__;
	  /*              le produit scalaire ci-dessus va dans r((j-1)*mm1+i) 
	   */
	  r__[nij] = ps;
	  /* L30: */
	}
    }
  /* 
   * 
   */
 L50:
  nta2 = *nta + 2;
  /* 
   *         remplissage des nouvelles donnees 
   * 
   */
  if (nta2 > nt1)
    {
      goto L100;
    }
  i__2 = nt1;
  for (kk = nta2; kk <= i__2; ++kk)
    {
      mekk = (kk - 1) * *mm1;
      e[kk] = 1.;
      r__[mekk + 1] = 0.;
      nj = (kk - 2) * *n;
      i__1 = kk;
      for (i__ = 2; i__ <= i__1; ++i__)
	{
	  ni = (i__ - 2) * *n;
	  /* 
	   *            produit scalaire de g(kk-1) avec g(i-1) 
	   *            pour kk-1=nta+1,ntot et i-1=1,kk-1 
	   * 
	   */
	  (*prosca) (n, &p[ni + 1], &p[nj + 1], &ps,optim_data);
	  nij = mekk + i__;
	  /*              le produit scalaire ci-dessus va dans r((kk-1)*mm1+i) 
	   */
	  /* L70: */
	  r__[nij] = ps;
	}
    }
  /* 
   *         remplissage de la contrainte d'inegalite 
   *              (tout entiere sauf l'ecart) 
   * 
   */
  i__1 = nt1;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      /* L80: */
      a[i__] = alfa[i__ - 1];
    }
 L100:
  return 0;
}				/* fremf2_ */

/* 
 *         cette subroutine remplit les donnees pour fprf2 
 *         (produits scalaires et 2 contraintes lineaires) 
 * 
 *            de 1 a ntot +1  si iflag=0 
 *            de nta+1 +1 a ntot +1 sinon 
 * 
 *            (le +1 est du a l'ecart, place en premier) 
 * 
 *         p contient les opposes des gradients a la queue leu leu 
 *         -g(1), -g(2),..., -g(ntot) soit ntot*n coordonnees 
 * 
 */

int optim_fremf1 (opt_prosca prosca, int *iflag, int *n, int *ntot, int *nta, int *mm1,
		  double *p, double *alfa, double *e, double *a, double *r__,
		  opt_simul_data *optim_data)
{
  int i__1, i__2;
  int mekk, i__, j, jj, kk, ni, nj, nt1, mej, nij, nta1, nta2;

  --alfa;
  --a;
  --e;
  --p;
  --r__;

  nt1 = *ntot + 1;
  nta1 = *nta + 1;
  if (*iflag > 0)
    {
      goto L50;
    }
  /* 
   *               remplissage des anciennes donnees 
   *         (produits scalaires, ecart et contrainte d'egalite) 
   * 
   */
  i__1 = *ntot;
  for (j = 1; j <= i__1; ++j)
    {
      jj = (j - 1) * *mm1 + 1;
      /* L10: */
      r__[jj] = 0.;
    }
  a[1] = 1.;
  e[1] = 0.;
  if (nta1 == 1)
    {
      goto L50;
    }
  i__1 = nta1;
  for (j = 2; j <= i__1; ++j)
    {
      e[j] = 1.;
      nj = (j - 2) * *n;
      mej = (j - 1) * *mm1;
      i__2 = j;
      for (i__ = 2; i__ <= i__2; ++i__)
	{
	  ni = (i__ - 2) * *n;
	  /* 
	   *            produit scalaire de g(i-1) avec g(j-1) 
	   *            pour j-1=1,nta et i-1=1,j-1 
	   *            qui va dans r((j-1)*mm1+i) 
	   * 
	   */
	  nij = mej + i__;
	  (*prosca) (n, &p[ni + 1], &p[nj + 1], &r__[nij], optim_data);
	  /* L30: */
	}
    }
  /* 
   * 
   */
 L50:
  nta2 = *nta + 2;
  /* 
   *         remplissage des nouvelles donnees 
   * 
   */
  if (nta2 > nt1)
    {
      goto L100;
    }
  i__2 = nt1;
  for (kk = nta2; kk <= i__2; ++kk)
    {
      mekk = (kk - 1) * *mm1;
      e[kk] = 1.;
      r__[mekk + 1] = 0.;
      nj = (kk - 2) * *n;
      i__1 = kk;
      for (i__ = 2; i__ <= i__1; ++i__)
	{
	  ni = (i__ - 2) * *n;
	  /* 
	   *            produit scalaire de g(kk-1) avec g(i-1) 
	   *            pour kk-1=nta+1,ntot et i-1=1,kk-1 
	   *            qui va dans r((kk-1)*mm1+i) 
	   * 
	   */
	  nij = mekk + i__;
	  (*prosca) (n, &p[ni + 1], &p[nj + 1], &r__[nij],optim_data);
	}
    }
  /* 
   *         remplissage de la contrainte d'inegalite 
   *              (tout entiere sauf l'ecart) 
   * 
   */
  i__1 = nt1;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      /* L80: */
      a[i__] = alfa[i__ - 1];
    }
 L100:
  return 0;
}				/* fremf1_ */


typedef struct fprf2c_ fprf2c;

struct fprf2c_ 
{
  double u1;
  int nc;
};




/* 
 *        the dimension is mm1*mm1 for r 
 *    *****  on entry  ***** 
 *      iflag=0-1  initialize on one subgradient (mu in) 
 *      iflag=2    "     "     "     "     "     "     " 
 *                 and strive to enter by priority the 
 *                 points of the previous corral at the 
 *                 beginning of the iterations. 
 *      iflag=3    initialize on the previous projection 
 *                 (with its corresponding corral) 
 * *****  on exit  ***** 
 *                 iflag=0    normal end 
 *                       1    old solution is already optimal 
 *                       2    constraints non feasible 
 *                       3    trying to enter a variable 
 *                            that is already in the corral 
 *                       4    starting to loop 
 *     imp > 5    one prints final information 
 *     imp > 6    one prints information at each iteration 
 *     imp > 7    one prints also 
 *               - at each iteration the choleski matrix 
 * ****   begin   **** 
 */

static int optim_fprf2 (int *iflag, int *ntot, int *nv, int *io, double *zero,
			double *s2, double *eps, double *al, int *imp, double *u,
			double *eta, int *mm1, int *jc, int *ic, double *r__, double *a,
			double *e, double *rr, double *xpr, double *y, double *w1,
			double *w2)
{
  static fprf2c fprf2c_1;/* XXXXX static ? */
  int c__21 = 21;
  int c__22 = 22;
  int c__23 = 23;
  int c__24 = 24;
  int c__25 = 25;
  int c__26 = 26;
  int c__27 = 27;
  int c__28 = 28;
  int c__29 = 29;
  int c__30 = 30;
  int c__31 = 31;
  int c__32 = 32;
  int c__34 = 34;

  /* System generated locals */
  int i__1, i__2, i__3, i__4;
  double d__1;

  /* Local variables */
  double gama;
  int mek01, mekk;
  double deps;
  int incr;
  double teta;
  int ment, i__, j, k, l;
  int itmax;
  double dzero, d1, d2, d3[1], d4[1];
  int j0, i2, i3, i4, i5[1], i1, j1, j2, k1, k0;
  double u2, v1, v2;
  int k00, jj, jk, kk;
  double ps, sp;
  int iterpr, nt1;
  double ps1;
  int nv1;
  double ps0, ps2, w1s, w2s;
  int mej, mek;
  double det, ps12, w12s;

  /* Parameter adjustments */
  --al;
  --w2;
  --w1;
  --y;
  --xpr;
  --rr;
  --e;
  --a;
  --ic;
  --jc;
  --r__;

  /* Function Body */
  iterpr = 0;
  nt1 = *ntot + 1;
  itmax = *ntot * 10;
  deps = *eps;
  incr = 0;
  k00 = 1;
  w1s = 0.;
  w2s = 0.;
  w12s = 0.;
  gama = .99;
  dzero = *zero * 10.;
  /*                    initial printouts 
   */
  if (*imp > 7)
    {
      optim_n1fc1o (io, &c__21, &nt1, mm1, &i3, &i4, i5, &deps, &d2, &a[1],
		    &r__[1]);
    }
  /* 
   *                    initial point 
   * 
   */
  /* L100: */
  if (*iflag != 3)
    {
      goto L110;
    }
  if (*imp > 6)
    {
      optim_n1fc1o (io, &c__22, nv, &i2, &i3, &i4, &jc[1], &d1, &d2, d3, d4);
    }
  j0 = nt1;
  ps = fprf2c_1.u1 * (a[nt1] - deps);
  ment = (nt1 - 1) * *mm1;
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      jk = ment + jc[k];
      /* L103: */
      ps += xpr[k] * r__[jk];
    }
  if (ps < *s2)
    {
      goto L107;
    }
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__23, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  *iflag = 1;
  return 0;
 L107:
  ++(*nv);
  ++fprf2c_1.nc;
  jc[*nv] = j0;
  iterpr = 1;
  goto L300;
 L110:
  if (*iflag <= 1)
    {
      goto L140;
    }
  /*       save the corral of previous call 
   */
  i__1 = nt1;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L120: */
      ic[i__] = 0;
    }
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      jk = jc[k];
      /* L130: */
      ic[jk] = 1;
    }
  ic[nt1] = 1;
  /*          initialize with one feasible gradient 
   */
 L140:
  jc[1] = 1;
  *nv = 2;
  fprf2c_1.nc = 1;
  jc[2] = 0;
  i__1 = nt1;
  for (j = 2; j <= i__1; ++j)
    {
      if (a[j] > deps)
	{
	  goto L150;
	}
      jc[2] = j;
    L150:
      ;
    }
  if (jc[2] > 0)
    {
      goto L160;
    }
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__24, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  *iflag = 2;
  return 0;
 L160:
  j = jc[2];
  rr[1] = 1.;
  jj = (j - 1) * *mm1 + j;
  ps = r__[jj] + 1.;
  if (ps > 0.)
    {
      goto L170;
    }
  *iflag = 3;
  return 0;
 L170:
  rr[2] = sqrt (ps);
  r__[2] = a[j];
  i__1 = nt1;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L180: */
      xpr[i__] = 0.;
    }
  xpr[1] = deps - a[j];
  xpr[2] = 1.;
  fprf2c_1.u1 = 0.;
  u2 = -r__[jj];
  if (*imp > 6)
    {
      optim_n1fc1o (io, &c__25, &j, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  /* 
   *                stopping criterion 
   * 
   */
 L200:
  ++iterpr;
  if (*imp > 6)
    {
      optim_n1fc1o (io, &c__26, nv, &i2, &i3, &i4, i5, &d1, &d2, d3, &xpr[1]);
    }
  if (iterpr <= itmax)
    {
      goto L205;
    }
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__27, &i1, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
    }
  *iflag = 4;
  return 0;
 L205:
  *s2 = -deps * fprf2c_1.u1 - u2;
  if (*s2 <= *eta)
    {
      goto L900;
    }
  sp = gama * *s2;
  /*                   first compute all the tests, 
   *           and test with the corral of previous call 
   */
  j0 = 0;
  i__1 = nt1;
  for (j = 2; j <= i__1; ++j)
    {
      ps = fprf2c_1.u1 * (a[j] - deps);
      i__2 = *nv;
      for (k = 1; k <= i__2; ++k)
	{
	  jj = jc[k];
	  if (jj == 1)
	    {
	      goto L210;
	    }
	  j1 = Max (j, jj);
	  j2 = Min (j, jj);
	  jj = (j1 - 1) * *mm1 + j2;
	  ps += xpr[k] * r__[jj];
	L210:
	  ;
	}
      y[j] = ps;
      if (*iflag != 2)
	{
	  goto L220;
	}
      if (ic[j] != 1)
	{
	  goto L220;
	}
      if (ps >= sp)
	{
	  goto L220;
	}
      j0 = j;
      sp = ps;
    L220:
      ;
    }
  if (j0 == 0)
    {
      goto L240;
    }
  if (sp >= gama * *s2)
    {
      goto L240;
    }
  ps1 = (d__1 = fprf2c_1.u1 * (deps - a[j0]), Abs (d__1));
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      j = jc[k];
      if (j == j0)
	{
	  goto L240;
	}
      if (j == 1)
	{
	  goto L230;
	}
      j1 = Max (j0, j);
      j2 = Min (j0, j);
      jj = (j1 - 1) * *mm1 + j2;
      ps1 += xpr[k] * (d__1 =
		       fprf2c_1.u1 * (deps * 2. - a[j]) + y[j] * 2. - r__[jj],
		       Abs (d__1));
    L230:
      ;
    }
  ps1 = ps1 * 1e3 * dzero;
  if (sp > *s2 - ps1)
    {
      goto L240;
    }
  ic[j0] = 0;
  goto L280;
  /*                    now the remaining ones 
   */
 L240:
  j0 = 0;
  sp = gama * *s2;
  i__1 = nt1;
  for (j = 2; j <= i__1; ++j)
    {
      if (*iflag == 2 && ic[j] == 1)
	{
	  goto L260;
	}
      if (y[j] >= sp)
	{
	  goto L260;
	}
      sp = y[j];
      j0 = j;
    L260:
      ;
    }
  if (j0 == 0)
    {
      goto L290;
    }
  ps1 = (d__1 = fprf2c_1.u1 * (deps - a[j0]), Abs (d__1));
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      j = jc[k];
      if (j == 1)
	{
	  goto L270;
	}
      j1 = Max (j0, j);
      j2 = Min (j0, j);
      jj = (j1 - 1) * *mm1 + j2;
      ps1 += xpr[k] * (d__1 =
		       fprf2c_1.u1 * (deps * 2. - a[j]) + y[j] * 2. - r__[jj],
		       Abs (d__1));
    L270:
      ;
    }
  ps1 = ps1 * 1e3 * dzero;
  if (sp > *s2 - ps1)
    {
      goto L290;
    }
 L280:
  ++fprf2c_1.nc;
  ++(*nv);
  jc[*nv] = j0;
  if (*imp > 6)
    {
      optim_n1fc1o (io, &c__28, &j0, &i2, &i3, &i4, i5, s2, &sp, d3, d4);
    }
  goto L300;
  /*        first set of optimality conditions satisfied 
   */
 L290:
  if (fprf2c_1.u1 >= -((double) (*nv)) * dzero)
    {
      goto L900;
    }
  j0 = 1;
  ++(*nv);
  jc[*nv] = 1;
  if (*imp > 6)
    {
      optim_n1fc1o (io, &c__29, &i1, &i2, &i3, &i4, i5, s2, &fprf2c_1.u1, d3,
		    d4);
    }
  /* 
   *              augmenting r 
   * 
   */
 L300:
  nv1 = *nv - 1;
  i__1 = nv1;
  for (k = 1; k <= i__1; ++k)
    {
      if (jc[k] != j0)
	{
	  goto L305;
	}
      if (*imp > 0)
	{
	  optim_n1fc1o (io, &c__30, &j0, &i2, &i3, &i4, i5, &d1, &d2, d3, d4);
	}
      *iflag = 3;
      return 0;
    L305:
      ;
    }
  j = jc[1];
  j1 = Max (j, j0);
  j2 = Min (j, j0);
  jj = (j1 - 1) * *mm1 + j2;
  r__[*nv] = (a[j] * a[j0] + e[j] * e[j0] + r__[jj]) / rr[1];
  ps0 = r__[*nv] * r__[*nv];
  if (nv1 == 1)
    {
      goto L330;
    }
  i__1 = nv1;
  for (k = 2; k <= i__1; ++k)
    {
      j = jc[k];
      j1 = Max (j, j0);
      j2 = Min (j, j0);
      jj = (j1 - 1) * *mm1 + j2;
      ps = a[j] * a[j0] + e[j] * e[j0] + r__[jj];
      k1 = k - 1;
      i__2 = k1;
      for (kk = 1; kk <= i__2; ++kk)
	{
	  j1 = (kk - 1) * *mm1 + k;
	  j2 = (kk - 1) * *mm1 + *nv;
	  /* L310: */
	  ps -= r__[j1] * r__[j2];
	}
      mek = k1 * *mm1 + *nv;
      r__[mek] = ps / rr[k];
      /* L320: */
      ps0 += r__[mek] * r__[mek];
    }
  jj = (j0 - 1) * *mm1 + j0;
  ps0 = a[j0] * a[j0] + e[j0] * e[j0] + r__[jj] - ps0;
  if (ps0 > 0.)
    {
      goto L330;
    }
  *iflag = 3;
  return 0;
 L330:
  rr[*nv] = sqrt (ps0);
  if (iterpr <= 1)
    {
      goto L400;
    }
  incr = 1;
  k00 = *nv;
  /* 
   *         solving the corral-system 
   * 
   */
 L400:
  k = k00;
  if (k > *nv)
    {
      goto L430;
    }
  if (*imp > 7)
    {
      optim_n1fc1o (io, &c__31, nv, mm1, &i3, &i4, i5, &d1, &d2, &rr[1],
		    &r__[1]);
    }
 L410:
  j = jc[k];
  ps1 = a[j];
  ps2 = e[j];
  if (k == 1)
    {
      goto L420;
    }
  k1 = k - 1;
  i__1 = k1;
  for (kk = 1; kk <= i__1; ++kk)
    {
      jj = (kk - 1) * *mm1 + k;
      ps0 = r__[jj];
      ps1 -= ps0 * w1[kk];
      /* L415: */
      ps2 -= ps0 * w2[kk];
    }
 L420:
  ps0 = rr[k];
  w1[k] = ps1 / ps0;
  w2[k] = ps2 / ps0;
  ++k;
  if (k <= *nv)
    {
      goto L410;
    }
  /*               two-two system 
   */
 L430:
  k = 1;
  if (incr == 1)
    {
      k = *nv;
    }
 L440:
  w1s += w1[k] * w1[k];
  w2s += w2[k] * w2[k];
  w12s += w1[k] * w2[k];
  ++k;
  if (k <= *nv)
    {
      goto L440;
    }
  det = w1s * w2s - w12s * w12s;
  ps2 = w2s * deps - w12s;
  ps1 = w1s - w12s * deps;
  /* L450: */
  v1 = ps2 / det;
  v2 = ps1 / det;
  /* L460: */
  fprf2c_1.u1 = deps - v1;
  u2 = 1. - v2;
  if (*nv == fprf2c_1.nc + 1)
    {
      fprf2c_1.u1 = 0.;
    }
  /*                 backward 
   */
  y[*nv] = (v1 * w1[*nv] + v2 * w2[*nv]) / rr[*nv];
  if (*nv == 1)
    {
      goto L500;
    }
  i__1 = *nv;
  for (l = 2; l <= i__1; ++l)
    {
      k = *nv - l + 1;
      k1 = k + 1;
      ps = v1 * w1[k] + v2 * w2[k];
      mek = (k - 1) * *mm1;
      i__2 = *nv;
      for (kk = k1; kk <= i__2; ++kk)
	{
	  mej = mek + kk;
	  /* L470: */
	  ps -= r__[mej] * y[kk];
	}
      /* L480: */
      y[k] = ps / rr[k];
    }
  /* 
   *               test for positivity 
   * 
   */
 L500:
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      if (y[k] <= 0.)
	{
	  goto L550;
	}
      /* L530: */
    }
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      /* L540: */
      xpr[k] = y[k];
    }
  goto L200;
  /*          interpolating between x and y 
   */
 L550:
  teta = 0.;
  k0 = k;
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      if (y[k] >= 0.)
	{
	  goto L560;
	}
      ps = y[k] / (y[k] - xpr[k]);
      if (teta >= ps)
	{
	  goto L560;
	}
      teta = ps;
      k0 = k;
    L560:
      ;
    }
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      ps = teta * xpr[k] + (1. - teta) * y[k];
      if (ps <= dzero)
	{
	  ps = 0.;
	}
      /* L570: */
      xpr[k] = ps;
    }
  if (*imp <= 6)
    {
      goto L600;
    }
  ps1 = 0.;
  ps2 = 0.;
  i__1 = *nv;
  for (k = 1; k <= i__1; ++k)
    {
      i__2 = *nv;
      for (kk = 1; kk <= i__2; ++kk)
	{
	  /*Computing MAX 
	   */
	  i__3 = jc[k], i__4 = jc[kk];
	  j1 = Max (i__3, i__4);
	  /*Computing MIN 
	   */
	  i__3 = jc[k], i__4 = jc[kk];
	  j2 = Min (i__3, i__4);
	  jj = (j1 - 1) * *mm1 + j2;
	  ps1 += xpr[k] * xpr[kk] * r__[jj];
	  ps2 += y[k] * y[kk] * r__[jj];
	  /* L580: */
	}
    }
  /* 
   *                 compressing the corral 
   * 
   */
 L600:
  --(*nv);
  incr = 0;
  k00 = k0;
  w1s = 0.;
  w2s = 0.;
  w12s = 0.;
  l = jc[k0];
  if (l != 1)
    {
      --fprf2c_1.nc;
    }
  if (*imp > 6)
    {
      optim_n1fc1o (io, &c__32, &k0, &l, &i3, &i4, i5, &y[k0], &ps1, &ps2,
		    d4);
    }
  if (k0 > *nv)
    {
      goto L400;
    }
  k1 = k0 - 1;
  i__2 = *nv;
  for (k = k0; k <= i__2; ++k)
    {
      xpr[k] = xpr[k + 1];
      if (k0 == 1)
	{
	  goto L620;
	}
      i__1 = k1;
      for (kk = 1; kk <= i__1; ++kk)
	{
	  mek = (kk - 1) * *mm1 + k;
	  /* L610: */
	  r__[mek] = r__[mek + 1];
	}
    L620:
      jc[k] = jc[k + 1];
    }
  xpr[*nv + 1] = 0.;
 L630:
  mek = (k0 - 1) * *mm1 + k0 + 1;
  ps = r__[mek];
  ps12 = rr[k0 + 1];
  ps0 = sqrt (ps * ps + ps12 * ps12);
  ps /= ps0;
  ps12 /= ps0;
  rr[k0] = ps0;
  if (k0 == *nv)
    {
      goto L400;
    }
  k1 = k0 + 1;
  mek01 = (k0 - 1) * *mm1;
  mek = k0 * *mm1;
  mekk = mek - *mm1;
  i__2 = *nv;
  for (k = k1; k <= i__2; ++k)
    {
      j1 = mekk + k;
      j2 = mek + k;
      r__[j1] = ps * r__[j1 + 1] + ps12 * r__[j2 + 1];
      if (k > k1)
	{
	  r__[j2] = ps2;
	}
      /* L640: */
      ps2 = -ps12 * r__[j1 + 1] + ps * r__[j2 + 1];
    }
  r__[j2 + 1] = ps2;
  ++k0;
  goto L630;
  /* 
   *                     *** finished *** 
   * 
   */
 L900:
  *iflag = 0;
  i__2 = *ntot;
  for (j = 1; j <= i__2; ++j)
    {
      /* L930: */
      al[j] = 0.;
    }
  i__2 = *nv;
  for (k = 1; k <= i__2; ++k)
    {
      j = jc[k] - 1;
      if (j != 0)
	{
	  al[j] = xpr[k];
	}
      /* L940: */
    }
  *u = fprf2c_1.u1;
  if (*imp <= 5)
    {
      return 0;
    }
  optim_n1fc1o (io, &c__34, &fprf2c_1.nc, nv, &i3, &i4, &jc[1], s2, &sp,
		&fprf2c_1.u1, d4);
  return 0;
}				/* fprf2_ */



/*
 *
 *
 *subroutine effectuant une recherche lineaire sur 0 tmax 
 *partant du point xn dans la direction d. 
 *sous l'hypothese d'hemiderivabilite, donne 
 *un pas serieux, bloque, nul ou semi serieux-nul (2 gradients). 
 *necessite fpn < 0 estimant la derivee a l'origine. 
 *appelle simul systematiquement avec indic = 4 
 * 
 * logic 
 *       0          descente serieuse 
 *       1          descente bloquee 
 *       2          pas semiserieux-nul 
 *       3          pas nul, enrichissement du faiseau 
 *       4          nap > napmax 
 *       5          retour a l'utilisateur 
 *       6          non hemi-derivable (au-dela de dx) 
 *       < 0        contrainte implicite active 
 * 
 *       imp 
 *                  =0 pas d'impressions 
 *                  >0 message en cas de fin anormale 
 *                  >3 informations pour chaque essai de t 
 *           ---------------------------------------- 
 *fait appel aux subroutines: 
 *-------simul(indic,n,x,f,g,izs,rzs,dzs) 
 *-------prosca(n,u,v,ps,izs,rzs,dzs) 
 * 
 */

static int optim_nlis2 (opt_simul simul, opt_prosca prosca, int *n, double *xn, double *fn,
			double *fpn, double *t, double *tmin, double *tmax, double *d__,
			double *d2, double *g, double *gd, double *amd, double *amf,
			int *imp, int *io, int *logic, int *nap, int *napmax, double *x,
			double *tol, double *a, double *tps, double *tnc, double *gg,
			opt_simul_data *optim_data)
{
  int c__35 = 35;
  int c__36 = 36;
  int c__37 = 37;
  int c__38 = 38;
  int c__39 = 39;
  int c__40 = 40;
  int c__41 = 41;
  int c__42 = 42;
  int i__1;
  double d__1, d__2, d__3, d__4;
  /* Local variables */
  double tesd, sign, tesf, anum, test, f;
  int i__;
  double p;
  int indic;
  double z__, d1, d3[1], d4[1];
  int i1, i2, i3, i4, i5[1];
  double z1, fa, fd, fg;
  int mc;
  double ta, fp;
  int indica;
  double td;
  int indicd;
  double tg, tc;
  int mp;
  double ps, tp, fpa, den, ffn, fpd, fpg;

  /* Parameter adjustments */
  --gg;
  --x;
  --gd;
  --g;
  --d__;
  --xn;

  /* Function Body */
  tesf = *amf * *fpn;
  tesd = *amd * *fpn;
  td = 0.;
  tg = 0.;
  fg = *fn;
  fpg = *fpn;
  ta = 0.;
  fa = *fn;
  fpa = *fpn;
  indica = 1;
  *logic = 0;
  /*         elimination d'un t initial ridiculement petit 
   */
  if (*t > *tmin)
    {
      goto L20;
    }
  *t = *tmin;
  if (*t <= *tmax)
    {
      goto L20;
    }
  if (*imp > 0)
    {
      optim_n1fc1o (io, &c__35, &i1, &i2, &i3, &i4, i5, &d1, d2, d3, d4);
    }
  *tmin = *tmax;
 L20:
  if (*fn + *t * *fpn < *fn + *t * .9 * *fpn)
    {
      goto L30;
    }
  *t *= 2.;
  goto L20;
  /* 
   */
 L30:
  if (*t < *tmax)
    {
      goto L40;
    }
  *t = *tmax;
  *logic = 1;
 L40:
  if (*imp >= 4)
    {
      optim_n1fc1o (io, &c__36, &i1, &i2, &i3, &i4, i5, fpn, d2, tmin, tmax);
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L50: */
      x[i__] = xn[i__] + *t * d__[i__];
    }
  /* 
   *                          boucle 
   * 
   */
 L100:
  ++(*nap);
  if (*nap <= *napmax)
    {
      goto L150;
    }
  /*               sortie par maximum de simulations 
   */
  *logic = 4;
  if (*imp >= 4)
    {
      optim_n1fc1o (io, &c__37, nap, &i2, &i3, &i4, i5, &d1, d2, d3, d4);
    }
  if (tg == 0.)
    {
      goto L999;
    }
  *fn = fg;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      g[i__] = gg[i__];
      /* L120: */
      xn[i__] += tg * d__[i__];
    }
  goto L999;
 L150:
  indic = 4;
  (*simul) (&indic, n, &x[1], &f, &g[1],optim_data);
  if (indic != 0)
    {
      goto L200;
    }
  /* 
   *               arret demande par l'utilisateur 
   */
  *logic = 5;
  *fn = f;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L170: */
      xn[i__] = x[i__];
    }
  if (*imp >= 4)
    {
      optim_n1fc1o (io, &c__38, &i1, &i2, &i3, &i4, i5, &d1, d2, d3, d4);
    }
  goto L999;
  /* 
   *               les tests elementaires sont faits, on y va 
   *               tout d'abord, ou en sommes nous ? 
   * 
   */
 L200:
  if (indic > 0)
    {
      goto L210;
    }
  td = *t;
  indicd = indic;
  *logic = 0;
  if (*imp >= 4)
    {
      optim_n1fc1o (io, &c__39, &indic, &i2, &i3, &i4, i5, t, d2, d3, d4);
    }
  *t = tg + (td - tg) * .1;
  goto L905;
  /* 
   *               calcul de la derivee directionnelle h'(t) 
   */
 L210:
  (*prosca) (n, &g[1], &d__[1], &fp, optim_data);
  /* 
   *        test de descente (premiere inegalite pour un pas serieux) 
   */
  ffn = f - *fn;
  if (ffn < *t * tesf)
    {
      goto L300;
    }
  td = *t;
  fd = f;
  fpd = fp;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L230: */
      gd[i__] = g[i__];
    }
  indicd = indic;
  *logic = 0;
  if (*imp >= 4)
    {
      optim_n1fc1o (io, &c__40, &i1, &i2, &i3, &i4, i5, t, &ffn, &fp, d4);
    }
  if (tg != 0.)
    {
      goto L500;
    }
  /*               tests pour un pas nul (si tg=0) 
   */
  if (fpd < tesd)
    {
      goto L500;
    }
  *tps = *fn - f + td * fpd;
  *tnc = *d2 * td * td;
  /*Computing MAX 
   */
  d__1 = *a * *tnc;
  p = Max (d__1, *tps);
  if (p > *tol)
    {
      goto L500;
    }
  *logic = 3;
  goto L999;
  /* 
   *                   descente 
   */
 L300:
  if (*imp >= 4)
    {
      optim_n1fc1o (io, &c__41, &i1, &i2, &i3, &i4, i5, t, &ffn, &fp, d4);
    }
  /* 
   *        test de derivee (deuxieme inegalite pour un pas serieux) 
   */
  if (fp < tesd)
    {
      goto L320;
    }
  /* 
   *               sortie, le pas est serieux 
   */
  *logic = 0;
  *fn = f;
  *fpn = fp;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L310: */
      xn[i__] = x[i__];
    }
  goto L999;
  /* 
   */
 L320:
  if (*logic == 0)
    {
      goto L350;
    }
  /* 
   *               sortie par descente bloquee 
   */
  *fn = f;
  *fpn = fp;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L330: */
      xn[i__] = x[i__];
    }
  goto L999;
  /* 
   *               on a une descente 
   */
 L350:
  tg = *t;
  fg = f;
  fpg = fp;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L360: */
      gg[i__] = g[i__];
    }
  /* 
   */
  if (td != 0.)
    {
      goto L500;
    }
  /*               extrapolation 
   */
  ta = *t;
  *t = tg * 9.;
  z__ = *fpn + fp * 3. - ffn * 4. / tg;
  if (z__ > 0.)
    {
      /*Computing MIN 
       *Computing MAX 
       */
      d__3 = 1., d__4 = -fp / z__;
      d__1 = *t, d__2 = tg * Max (d__3, d__4);
      *t = Min (d__1, d__2);
    }
  *t = tg + *t;
  if (*t < *tmax)
    {
      goto L900;
    }
  *logic = 1;
  *t = *tmax;
  goto L900;
  /* 
   *               interpolation 
   * 
   */
 L500:
  if (indica > 0 && indicd > 0)
    {
      goto L510;
    }
  ta = *t;
  *t = tg * .9 + td * .1;
  goto L900;
 L510:
  test = (td - tg) * .1;
  /*               approximation cubique 
   */
  ps = fp + fpa - (fa - f) * 3. / (ta - *t);
  z1 = ps * ps - fp * fpa;
  if (z1 >= 0.)
    {
      goto L520;
    }
  if (fp < 0.)
    {
      tc = td;
    }
  if (fp >= 0.)
    {
      tc = tg;
    }
  goto L600;
 L520:
  z1 = sqrt (z1);
  if (*t - ta < 0.)
    {
      z1 = -z1;
    }
  sign = (*t - ta) / (d__1 = *t - ta, Abs (d__1));
  if ((ps + fp) * sign > 0.)
    {
      goto L550;
    }
  den = ps * 2. + fp + fpa;
  anum = z1 - fp - ps;
  if ((d__1 = (*t - ta) * anum, Abs (d__1)) >= (td - tg) * Abs (den))
    {
      goto L530;
    }
  tc = *t + anum * (ta - *t) / den;
  goto L600;
 L530:
  tc = td;
  goto L600;
 L550:
  tc = *t + fp * (ta - *t) / (ps + fp + z1);
 L600:
  mc = 0;
  if (tc < tg)
    {
      mc = -1;
    }
  if (tc > td)
    {
      mc = 1;
    }
  /*Computing MAX 
   */
  d__1 = tc, d__2 = tg + test;
  tc = Max (d__1, d__2);
  /*Computing MIN 
   */
  d__1 = tc, d__2 = td - test;
  tc = Min (d__1, d__2);
  /*               approximation polyhedrique 
   */
  ps = fpd - fpg;
  if (ps != 0.)
    {
      goto L620;
    }
  tp = (td + tg) * .5;
  goto L650;
 L620:
  tp = (fg - fpg * tg - (fd - fpd * td)) / ps;
 L650:
  mp = 0;
  if (tp < tg)
    {
      mp = -1;
    }
  if (tp > td)
    {
      mp = 1;
    }
  /*Computing MAX 
   */
  d__1 = tp, d__2 = tg + test;
  tp = Max (d__1, d__2);
  /*Computing MIN 
   */
  d__1 = tp, d__2 = td - test;
  tp = Min (d__1, d__2);
  /*               nouveau t par approximation cp complete securisee 
   */
  ta = *t;
  if (mc == 0 && mp == 0)
    {
      *t = Min (tc, tp);
    }
  if (mc == 0 && mp != 0)
    {
      *t = tc;
    }
  if (mc != 0 && mp == 0)
    {
      *t = tp;
    }
  if (mc == 1 && mp == 1)
    {
      *t = td - test;
    }
  if (mc == -1 && mp == -1)
    {
      *t = tg + test;
    }
  if (mc * mp == -1)
    {
      *t = (tg + td) * .5;
    }
  /* 
   *                fin de boucle 
   * 
   */
 L900:
  fa = f;
  fpa = fp;
 L905:
  indica = indic;
  /*                peut-on faire logic=2 ? 
   */
  if (td == 0.)
    {
      goto L920;
    }
  if (indicd < 0)
    {
      goto L920;
    }
  if (td - tg > *tmin * 10.)
    {
      goto L920;
    }
  if (fpd < tesd)
    {
      goto L920;
    }
  *tps = fg - fd + (td - tg) * fpd;
  *tnc = *d2 * (td - tg) * (td - tg);
  /*Computing MAX 
   */
  d__1 = *a * *tnc;
  p = Max (d__1, *tps);
  if (p > *tol)
    {
      goto L920;
    }
  /*              sortie par pas semiserieux-nul 
   */
  *logic = 2;
  *fn = fg;
  *fpn = fpg;
  *t = tg;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      xn[i__] += tg * d__[i__];
      /* L910: */
      g[i__] = gg[i__];
    }
  goto L999;
  /* 
   *               test d'arret sur la proximite de tg et td 
   * 
   */
 L920:
  if (td == 0.)
    {
      goto L990;
    }
  if (td - tg <= *tmin)
    {
      goto L950;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      z__ = xn[i__] + *t * d__[i__];
      if (z__ != x[i__] && z__ != xn[i__])
	{
	  goto L990;
	}
      /* L930: */
    }
  /*               arret sur dx ou de secours 
   */
 L950:
  *logic = 6;
  if (indicd < 0)
    {
      *logic = indicd;
    }
  if (tg == 0.)
    {
      goto L970;
    }
  *fn = fg;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      xn[i__] += tg * d__[i__];
      /* L960: */
      g[i__] = gg[i__];
    }
 L970:
  if (*imp <= 0)
    {
      goto L999;
    }
  if (*logic < 0)
    {
      optim_n1fc1o (io, &c__42, logic, &i2, &i3, &i4, i5, &d1, d2, d3, d4);
    }
  if (*logic == 6)
    {
      optim_n1fc1o (io, &c__42, &i1, &i2, &i3, &i4, i5, &d1, d2, d3, d4);
    }
  goto L999;
  /* 
   *               recopiage de x et boucle 
   */
 L990:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L995: */
      x[i__] = xn[i__] + *t * d__[i__];
    }
  goto L100;
  /* 
   */
 L999:
  return 0;
}				/* nlis2_ */



static int optim_n1fc1o (int *unit, int *job, int *i1, int *i2, int *i3, int *i4,
			 int *i5, double *d1, double *d2, double *d3, double *d4)
{
  int i__1, i__2, i__3, i__4;

  /* Local variables */
  double deps;
  int iter;
  double epsm;
  int nsim;
  double tmin, tmax;
  int ntot;
  double diam1;
  int i__, j, k, l, n;
  double t, u;
  int indic, logic, memax, lunit, j0, k1, k0;
  double s2, u1, ap;
  int nc;
  double fn;
  int jj, kk;
  double fp;
  int  ll, ln, nn;
  double ro, sp;
  int nv;
  int mm1, nt1;
  double yk0, ps1, ps2, ffn;
  int mej;
  int nap;
  double tnc, fpn, eps;
  int ndz, niz;
  double tps;
  int nrz;
  double eta2;

  /*    impression des traces 
   *    Copyright INRIA 
   * 
   */
  /* Parameter adjustments */
  --d4;
  --d3;
  --i5;

  /* Function Body */
  lunit = *unit;
  /* 
   */
  switch (*job)
    {
    case 1:
      goto L11;
    case 2:
      goto L12;
    case 3:
      goto L13;
    case 4:
      goto L14;
    case 5:
      goto L15;
    case 6:
      goto L16;
    case 7:
      goto L17;
    case 8:
      goto L18;
    case 9:
      goto L19;
    case 10:
      goto L20;
    case 11:
      goto L21;
    case 12:
      goto L22;
    case 13:
      goto L23;
    case 14:
      goto L24;
    case 15:
      goto L25;
    case 16:
      goto L26;
    case 17:
      goto L27;
    case 18:
      goto L28;
    case 19:
      goto L29;
    case 20:
      goto L30;
    case 21:
      goto L31;
    case 22:
      goto L32;
    case 23:
      goto L33;
    case 24:
      goto L34;
    case 25:
      goto L35;
    case 26:
      goto L36;
    case 27:
      goto L37;
    case 28:
      goto L38;
    case 29:
      goto L39;
    case 30:
      goto L40;
    case 31:
      goto L41;
    case 32:
      goto L42;
    case 33:
      goto L43;
    case 34:
      goto L44;
    case 35:
      goto L45;
    case 36:
      goto L46;
    case 37:
      goto L47;
    case 38:
      goto L48;
    case 39:
      goto L49;
    case 40:
      goto L50;
    case 41:
      goto L51;
    case 42:
      goto L52;
    case 43:
      goto L53;
    case 44:
      goto L54;
    case 45:
      goto L55;
    case 46:
      goto L56;
    case 47:
      goto L57;
    case 48:
      goto L58;
    case 49:
      goto L59;
    case 50:
      goto L60;
    }
  /* 
   */
 L11:
  Sciprintf( "n1fc1   incorrect call");
  goto L100;
 L12:
  n = *i1;
  memax = *i2;
  niz = *i3;
  nrz = *i4;
  ndz = i5[1];
  Sciprintf("entry in n1fc1 . n=%d memax=%d\n", n, memax );
  Sciprintf("minimal array sizes iz(%d), rz(%d), dz(%d)\n",    niz,    nrz,   ndz );
  goto L100;
 L13:
  Sciprintf( "n1fc1 initial gradient norm is zero\n");
  goto L100;
 L14:
  goto L100;
 L15:
  iter = *i1;
  Sciprintf("n1fc1 end with iter =%d\n", *i1);
  goto L100;
 L16:
  Sciprintf( "n1fc1 Incorrect end of fprf2\n");
  goto L100;
 L17:
  eta2 = *d1;
  Sciprintf("'n1fc1 eta2 assigned to %10.3f\n",eta2);
  goto L100;
 L18:
  iter = *i1;
  nsim = *i2;
  fn = *d1;
  epsm = *d2;
  s2 = d3[1];
  Sciprintf("n1fc1 %d,%d,%16.7f,   convergence a %10.3f pres   (%9.2f)\n",
	    iter,
	    nsim,
	    fn,
	    epsm,
	    s2);
  goto L100;
 L19:
  Sciprintf( "n1fc1   normal end");
  goto L100;
 L20:
  Sciprintf(" ");
  goto L100;
 L21:
  iter = *i1;
  nsim = *i2;
  nv = *i3;
  fn = *d1;
  eps = *d2;
  s2 = d3[1];
  u = d4[1];
  Sciprintf("n1fc1  %d,%d,%14.7f,%10.2f,%10.2f,%10.2f,%d\n",
	    iter,
	    nsim,
	    fn,
	    eps,
	    s2,
	    u,
	    nv);
  goto L100;
 L22:
  ntot = *i1;
  Sciprintf( "n1fc1  ponderation table\n");
  nn = ntot / 7;
  if (nn * 7 < ntot)
    {
      ++nn;
    }
  l = 0;
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* Computing MIN    */
      i__2 = 7, i__3 = ntot - l;
      ln = Min (i__2, i__3);
      i__2 = ln;
      for (j = 1; j <= i__2; ++j)
	{
	  Sciprintf("%10.3f ,", d4[l + j]);
	}
      l += 7;
      /* L2201: */
    }
 L23:
  Sciprintf( "n1fc1  la direction ne pivote plus\n");
  goto L100;
 L24:
  Sciprintf( "n1fc1  end (dxmin reached)\n");
  goto L100;
 L25:
  Sciprintf( "n1fc1  end (nsim reached)\n");
  goto L100;
 L26:
  Sciprintf( "n1fc1  end (indic=0)\n");
  goto L100;
 L27:
  Sciprintf( "n1fc1  warning txmax reached, reduce scale\n");
  goto L100;
 L28:
  diam1 = *d1;
  eta2 = *d2;
  ap = d3[1];
  Sciprintf("n1fc1, diam1=%10.3f, eta2=%10.3f, ap=%10.3f\n",  diam1,  eta2,  ap);
  goto L100;
 L29:
  iter = *i1;
  nsim = *i2;
  ntot = *i3;
  fn = *d1;
  Sciprintf("n1fc1, %d,%d,%16.7f,   faisceau reduit a,%d, gradients\n",  iter,  nsim,  fn,  ntot );
  goto L100;
 L30:
  logic = *i1;
  ro = *d1;
  tps = *d2;
  tnc = d3[1];
  Sciprintf("n1fc1, logic=%d,ro=%10.3f, tps=%10.3f, tnc=%10.3f\n",  logic,  ro,  tps,  tnc  );
  goto L100;
  /*    ================== 
   *    MESSAGES de frepf2 
   *    ================== 
   */
 L31:
  nt1 = *i1;
  mm1 = *i2;
  deps = *d1;
  Sciprintf( "a = ");
  nn = nt1 / 10;
  if (nn * 10 < nt1)
    {
      ++nn;
    }
  l = 0;
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*Computing MIN 
       */
      i__2 = 10, i__3 = nt1 - l;
      ln = Min (i__2, i__3);
      i__2 = ln;
      for (j = 1; j <= i__2; ++j)
	{
	  Sciprintf("%10.3f ", d3[l + j]);
	}
      l += 10;
    }
  Sciprintf("epsilon =%10.3f\n",deps);
  Sciprintf("g,g = ");
  i__1 = nt1;
  for (j = 1; j <= i__1; ++j)
    {
      mej = (j - 1) * mm1;
      nn = j / 10;
      if (nn * 10 < j)
	{
	  ++nn;
	}
      l = 0;
      i__2 = nn;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  /*Computing MIN 
	   */
	  i__3 = 10, i__4 = j - l;
	  ln = Min (i__3, i__4);
	  i__3 = ln;
	  for (jj = 1; jj <= i__3; ++jj)
	    {
	      Sciprintf("%10.3f ,", d4[mej + l + jj]);
	    }
	  l += 10;
	}
      /* L3103: */
    }
  goto L100;
 L32:
  nv = *i1;
  Sciprintf( "initial corral");
  i__1 = nv;
  for (k = 1; k <= i__1; ++k)
    {
      Sciprintf("%d ", i5[k]);
    }
  goto L100;
 L33:
  Sciprintf("error from fprf2. old solution already optimal\n");
  goto L100;
 L34:
  Sciprintf( "epsilon smaller than a\n");
  goto L100;
 L35:
  j = *i1;
  Sciprintf("start with variables 1 and %d\n",j);
  goto L100;
 L36:
  nv = *i1;
  Sciprintf( "x = ");
  nn = nv / 10;
  if (nn * 10 < nv)
    {
      ++nn;
    }
  l = 0;
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*Computing MIN 
       */
      i__2 = 10, i__3 = nv - l;
      ln = Min (i__2, i__3);
      i__2 = ln;
      for (j = 1; j <= i__2; ++j)
	{
	  Sciprintf("%10.3f ", d4[l + j]);
	}
      l += 10;
      /* L3601: */
    }
  goto L100;
 L37:
  Sciprintf( "fprf2 is apparently looping\n");
  goto L100;
 L38:
  j0 = *i1;
  s2 = *d1;
  sp = *d2;
  Sciprintf("(s,s)=,%12.4f,  variable,%d, (,%12.4f,) coming in.\n",  s2,  j0,  sp  );
  goto L100;
 L39:
  s2 = *d1;
  u1 = *d2;
  Sciprintf("(s,s)=,%12.4f,  u1=,%12.3f,  variable 1 coming in.\n",  s2,  u1  );
  goto L100;
 L40:
  Sciprintf("duplicate variable %d\n",j0);
  goto L100;
 L41:
  nv = *i1;
  mm1 = *i2;
  /*    d3=rr,d4=r 
   */
  Sciprintf("cholesky %11.3f\n",  d3[1]);
  if (nv >= 2)
    {
      i__1 = nv;
      for (ll = 2; ll <= i__1; ++ll)
	{
	  k1 = ll - 1;
	  nn = k1 / 10;
	  if (nn * 10 < k1)
	    {
	      ++nn;
	    }
	  l = 0;
	  if (nn > 1)
	    {
	      i__2 = nn - 1;
	      for (i__ = 1; i__ <= i__2; ++i__)
		{
		  /*Computing MIN 
		   */
		  i__3 = 10, i__4 = k1 - l;
		  ln = Min (i__3, i__4);
		  i__3 = nn;
		  for (kk = 1; kk <= i__3; ++kk)
		    {
		      Sciprintf("%10.3f ", d4[(l + kk - 1) * mm1 + ll]);
		    }
		  l += 10;
		  /* L4102: */
		}
	    }
	  i__2 = nn;
	  for (kk = 1; kk <= i__2; ++kk)
	    {
	      Sciprintf("%10.3f ", d4[(l + kk - 1) * mm1 + ll]);
	    }
	  Sciprintf("%10.3f ", d3[ll]);
	}
    }
  goto L100;
 L42:
  k0 = *i1;
  l = *i2;
  yk0 = *d1;
  ps1 = *d2;
  ps2 = d3[1];
  Sciprintf("variable %d  (%d) =%11.3f  going out. feasible (s,s)=%11.4f, unfeasible=%11.4f\n",
	    k0,
	    l,
	    yk0,
	    ps1,
	    ps2  );
  goto L100;
 L43:
  goto L100;
 L44:
  nc = *i1;
  nv = *i2;
  /*    jc=i5 
   */
  s2 = *d1;
  sp = *d2;
  u1 = d3[1];
  Sciprintf("finished with %d, gradients ,%d, variables./ (s,s)=%11.4f, test=%11.4f/ cost of ,%11.4f",
	    nc,
	    nv,
	    s2,
	    sp,
	    u1 );

  nn = nv / 20;
  if (nn * 10 < nv)
    {
      ++nn;
    }
  l = 0;
  i__1 = nn;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /*Computing MIN 
       */
      i__2 = 20, i__3 = nv - l;
      ln = Min (i__2, i__3);
      i__2 = ln;
      for (k = 1; k <= i__2; ++k)
	{
	  Sciprintf("%d ", i5[l + k]);
	}
      l += 20;
    }
  goto L100;
  /*    ================ 
   *    MESSAGE DE NLIS2 
   *    ================ 
   */
 L45:
  Sciprintf("nlis2,tmin force a tmax\n");
  goto L100;
 L46:
  fpn = *d1;
  tmin = d3[1];
  tmax = d4[1];
  Sciprintf( " ");
  Sciprintf("nlis2  fpn=%10.3f, d2=%9.2f,  tmin=%9.2f, tmax=%9.2f\n",
	    fpn,
	    (*d2),
	    tmin,
	    tmax);
  goto L100;
 L47:
  Sciprintf( " ");
  Sciprintf("nlis2, %d, simulations atteintes\n",  nap);
  goto L100;
 L48:
  Sciprintf( "Stop required by user\n");
  goto L100;
 L49:
  indic = *i1;
  t = *d1;
  Sciprintf("nlis2,i,%10.3f, indic=%d\n",  t,  indic  );
  goto L100;
 L50:
  t = *d1;
  ffn = *d2;
  fp = d3[1];
  Sciprintf("nlis2,i,%10.3f,%11.3f,%11.3f\n",  t,  ffn,  fp  );
  goto L100;
 L51:
  Sciprintf("nlis2,%13.3f,%11.3f,%11.3f,i\n",  t,  ffn,  fp  );
  goto L100;
 L52:
  logic = *i1;
  Sciprintf("nlis2 contrainte implicite %d active\n",  logic  );
  goto L100;
 L53:
  logic = *i1;
  Sciprintf( "nlis2  end (tmin reached)\n");
  goto L100;
 L54:
  goto L100;
 L55:
  goto L100;
 L56:
  goto L100;
 L57:
  goto L100;
 L58:
  goto L100;
 L59:
  goto L100;
 L60:
  goto L100;
  /* 
   */
 L100:
  return 0;
}





static int optim_frdf1 (opt_prosca prosca, int *n, int *ntot, int *ninf, int *kgrad,
			double *al, double *q, double *s, double *poids, double *aps,
			double *anc, int *mm1, double *r__, double *e, int *ic, opt_simul_data *optim_data)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int i__, j, k;
  double z__, z1, z2;
  int nj, nn;
  double ps;
  int nt1, njk;

  /*    Copyright INRIA 
   * 
   *             this subroutine reduces a nonconvex bundle 
   *             of size ntot in rn 
   *             to a size no greater than ninf 
   * 
   */
  /* Parameter adjustments */
  --s;
  --anc;
  --aps;
  --poids;
  --al;
  --q;
  --ic;
  --e;
  --r__;

  /* Function Body */
  if (*ntot <= *ninf)
    {
      goto L900;
    }
  if (*ninf > 0)
    {
      goto L100;
    }
  /* 
   *          pure gradient method 
   * 
   */
  *ntot = 0;
  *kgrad = 0;
  goto L900;
  /* 
   *         reduction to the corral 
   */
 L100:
  nt1 = 0;
  i__1 = *ntot;
  for (j = 1; j <= i__1; ++j)
    {
      if (al[j] == 0. && poids[j] != 0.)
	{
	  goto L150;
	}
      ++nt1;
      ic[nt1] = j;
      if (j == nt1)
	{
	  goto L130;
	}
      nj = *n * (j - 1);
      nn = *n * (nt1 - 1);
      i__2 = *n;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  ++nn;
	  ++nj;
	  /* L110: */
	  q[nn] = q[nj];
	}
      al[nt1] = al[j];
      poids[nt1] = poids[j];
      aps[nt1] = aps[j];
      anc[nt1] = anc[j];
      e[nt1 + 1] = e[j + 1];
    L130:
      if (poids[j] == 0.)
	{
	  *kgrad = nt1;
	}
      nn = nt1 * *mm1 + 1;
      nj = j * *mm1 + 1;
      i__2 = nt1;
      for (k = 1; k <= i__2; ++k)
	{
	  njk = nj + ic[k];
	  ++nn;
	  /* L140: */
	  r__[nn] = r__[njk];
	}
    L150:
      ;
    }
  *ntot = nt1;
  if (*ntot <= *ninf)
    {
      goto L900;
    }
  /* 
   *         corral still too large 
   *             save the near 
   * 
   */
  (*prosca) (n, &s[1], &s[1], &ps,optim_data);
  e[2] = 1.;
  z__ = 0.;
  z1 = 0.;
  z2 = 0.;
  i__1 = *ntot;
  for (k = 1; k <= i__1; ++k)
    {
      z1 += al[k] * aps[k];
      z2 += al[k] * anc[k];
      /* L370: */
      z__ += al[k] * poids[k];
    }
  aps[1] = z1;
  anc[1] = z2;
  poids[1] = z__;
  r__[*mm1 + 2] = ps;
  if (*ninf > 1)
    {
      goto L400;
    }
  *ntot = 1;
  *kgrad = 0;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L380: */
      q[i__] = s[i__];
    }
  goto L900;
  /*               save the gradient 
   */
 L400:
  nn = (*kgrad - 1) * *n;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      nj = *n + i__;
      ++nn;
      q[nj] = q[nn];
      /* L470: */
      q[i__] = s[i__];
    }
  (*prosca) (n, &q[*n + 1], &s[1], &ps,optim_data);
  e[3] = 1.;
  r__[(*mm1 << 1) + 2] = ps;
  (*prosca) (n, &q[*n + 1], &q[*n + 1], &ps,optim_data);
  r__[(*mm1 << 1) + 3] = ps;
  aps[2] = 0.;
  anc[2] = 0.;
  poids[2] = 0.;
  *kgrad = 2;
  *ntot = 2;
 L900:
  return 0;
}

