#include "optim.h"
static int optim_fretc1 (int *mode, int *n, int *nc, int *nr, double *h__, double *w,
			 int *indi, int *indic2);

static int optim_fmani1 (int *mode, int *n, double *d__, double *w, int *indi);

static int optim_fmc11e (double *a, int *n, double *z__, double *w, int *ir);

static int optim_fmc11a (double *a, int *n, double *z__, double *sig, double *w, int *ir,
			 int *mk, double *eps);

static int optim_fmlag1 (int *n, int *nr, double *a, double *z__, double *w);
static int optim_fmc11b (double *a, int *n, int *ir);

static int optim_fmc11z (double *a, int *n, int *nr, double *z__, 
			 double *sig, double *w,
			 int *ir, int *mk, double *eps);

static int optim_fcomp1 (int *indic2, int *ibloc, int *indi, double *h__, double *g,
			 double *d__, double *w, double *w1, int *n, int *nr, int *ncs,
			 double *dga, double *delta, double *prop, double *acc,
			 double *scale);

static int optim_n2qn1a (opt_simul simul, int *n, double *x, double *f, double *ga,
			 double *scale, double *acc, double *df1, int *mode, int *niter,
			 int *nsim, int *iprint, int *lp, double *h__, double *d__,
			 double *w, double *w1, double *g, double *binf, double *bsup,
			 int *indi, int *ibloc, int *iz, opt_simul_data *optim_data);

static int optim_fajc1 (int *n, int *nc, int *nr, double *h__, double *w, int *indi);

int optim_n2qn1 (opt_simul simul, int *n, double *x, double *f, double *g,
		 double *dxmin, double *df1, double *epsabs, int *imp, int *io,
		 int *mode, int *iter, int *nsim, double *binf, double *bsup,
		 int *iz, double *rz, opt_simul_data *optim_data)
{

  int i__1;
  /* Local variables */
  int i__;
  double s;
  int nindi, nd, ni, nw, nibloc, nga, nww, nww1;

  /* Parameter adjustments */
  --bsup;
  --binf;
  --dxmin;
  --g;
  --x;
  --iz;
  --rz;

  /* Function Body */
  /* L1000: */
  /* L1001: */
  /* L1002: */
  if (*imp == 0)
    {
      goto L10;
    }
  nw = *n * (*n + 9) / 2;
  ni = (*n << 1) + 1;
  Sciprintf("n2qn1: start dimension du probleme %d, mode=%d, niter=%d, nsim=%d,\n\t imp=%d df1=%9.2f, epsabs=%9.2f,\n\t dimensions minimales iz(%d) rz(%d)\n",
	    (*n),
	    (*mode),
	    (*iter),
	    (*nsim),
	    (*imp),
	    (*df1),
	    (*epsabs),
	    ni,
	    nw);
 L10:
  if (*n > 0 && (*df1 > 0. || *mode != 1) && *epsabs >= 0. && *mode >= 1
      && *mode <= 4 && *iter > 0 && *nsim > 0)
    {
      goto L100;
    }
  Sciprintf( "n2qn1: appel incoherent\n");
  *mode = 2;
  goto L999;
 L100:
  nd = *n * (*n + 1) / 2 + 1;
  nww = nd + *n;
  nww1 = nww + *n;
  nga = nww1 + *n;
  nindi = 1;
  nibloc = nindi + *n;
  ni = nibloc + *n;
  /* 
   *     calcul du test d arret 
   * 
   */
  s = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L110: */
      s += dxmin[i__] * dxmin[i__];
    }
  *epsabs *= sqrt (s / (double) (*n));
  optim_n2qn1a ( simul, n, &x[1], f, &g[1], &dxmin[1], epsabs, df1,
		 mode, iter, nsim, imp, io, &rz[1], &rz[nd], &rz[nww],
		 &rz[nww1], &rz[nga], &binf[1], &bsup[1], &iz[nindi],
		 &iz[nibloc], &iz[ni], optim_data);
  if (*imp > 0)
    {
      Sciprintf("n2qn1: end. norme du gradient projete %10.3f, mode=%d, niter=%d, nsim=%d\n",
		(*epsabs),
		(*mode),
		(*iter),
		(*nsim));
    }
 L999:
  return 0;
}				/* n2qn1_ */



static int optim_n2qn1a (opt_simul simul, int *n, double *x, double *f, double *ga,
			 double *scale, double *acc, double *df1, int *mode, int *niter,
			 int *nsim, int *iprint, int *lp, double *h__, double *d__,
			 double *w, double *w1, double *g, double *binf, double *bsup,
			 int *indi, int *ibloc, int *iz, opt_simul_data *optim_data)
{
  int c__1 = 1;
  double c_b109 = 0.;
  int c__0 = 0;

  /* System generated locals */
  int i__1;
  double d__1, d__2, d__3;
  /* Local variables */
  double alfa, beta, wiii;
  int nfun;
  double prop;
  double c__;
  int i__;
  int k;
  int indic, iecri;
  double delta;
  int logic;
  double z__;
  double theta;
  int isign;
  double romin, romax;
  int i1, k1, k2, indic1, indic2;
  double fa, df, bi, di, gi;
  int nc;
  double sc;
  int nh, ir, np, nr;
  double wi, xi, ro;
  double rocand;
  int nr1;
  double dga;
  int nca;
  double roa, dnr;
  int ncs;
  double wii;
  int itr;
  double acc1;
  int nrp1;
  double dgaa;


  /* 
   *     initialisations 
   * 
   */
  /* Parameter adjustments */
  --ibloc;
  --indi;
  --bsup;
  --binf;
  --g;
  --w1;
  --w;
  --d__;
  --scale;
  --ga;
  --x;
  --h__;
  --iz;

  alfa = .7;
  beta = .1;
  prop = 1.;
  nfun = 1;
  iecri = 0;
  itr = 0;
  np = *n + 1;
  indic2 = 1;
  logic = 0;
  df = *df1;
  if (*mode <= 3)
    {
      goto L1;
    }
  nr = iz[1];
  goto L400;
  /*               calcul des bornes actives 
   */
 L1:
  nr = 0;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (scale[i__] > 0.)
	{
	  goto L2;
	}
      *mode = 2;
      Sciprintf( "n2qn1   appel incoherent\n");
      goto L15;
    L2:
      bi = bsup[i__];
      if (x[i__] < bi - scale[i__])
	{
	  goto L4;
	}
      if (x[i__] <= bi)
	{
	  goto L3;
	}
      *mode = 2;
      Sciprintf( "n2qn1   appel incoherent\n");
      goto L15;
    L3:
      if (ga[i__] >= 0. && *mode == 1)
	{
	  goto L13;
	}
      ibloc[i__] = 1;
      goto L14;
    L4:
      bi = binf[i__];
      if (x[i__] > bi + scale[i__])
	{
	  goto L13;
	}
      if (x[i__] >= bi)
	{
	  goto L5;
	}
      *mode = 2;
      Sciprintf( "n2qn1   appel incoherent\n");
      goto L15;
    L5:
      if (ga[i__] <= 0. && *mode == 1)
	{
	  goto L13;
	}
      ibloc[i__] = -1;
      goto L14;
    L13:
      ++nr;
      ibloc[i__] = 0;
    L14:
      ;
    }
  goto L16;
 L15:
  *niter = 1;
  *acc = 0.;
  *nsim = 1;
  goto L999;
  /* 
   *     le point de depart est-il optimal? 
   * 
   */
 L16:
  c__ = 0.;
  dnr = sqrt ((double) nr);
  acc1 = *acc * dnr;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (ibloc[i__] != 0)
	{
	  goto L100;
	}
      gi = ga[i__] * scale[i__];
      c__ += gi * gi;
    L100:
      ;
    }
  c__ = sqrt (c__);
  if (c__ > acc1)
    {
      goto L200;
    }
  optim_fcomp1 (&indic2, &ibloc[1], &indi[1], &h__[1], &ga[1], &d__[1], &w[1],
		&w1[1], n, &nr, &ncs, &dga, &delta, &prop, acc, &scale[1]);
  if (ncs != 0)
    {
      goto L102;
    }
  itr = 1;
  *mode = 1;
  goto L900;
 L102:
  ibloc[ncs] = 0;
  ++nr;
  goto L200;
  /* 
   *        initialisation du hessien, en fonction de scale et de df1 
   * 
   */
 L200:
  switch (*mode)
    {
    case 1:
      goto L300;
    case 2:
      goto L310;
    case 3:
      goto L320;
    }
 L300:
  if (*df1 > 0.)
    {
      goto L301;
    }
  *mode = 2;
  Sciprintf( "n2qn1   appel incoherent\n");
  goto L15;
 L301:
  c__ = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (ibloc[i__] != 0)
	{
	  goto L302;
	}
      gi = ga[i__];
      sc = scale[i__];
      c__ += gi * gi * sc * sc;
    L302:
      ;
    }
  c__ = c__ * .5 / *df1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      sc = scale[i__];
      /* L303: */
      w[i__] = c__ / (sc * sc);
    }
  nh = *n * (*n + 1) / 2;
  i__1 = nh;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L304: */
      h__[i__] = 0.;
    }
  /*     permutation de la matrice h 
   */
  nr1 = nr + 1;
  k1 = 1;
  k2 = nr + 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (ibloc[i__] != 0)
	{
	  goto L305;
	}
      indi[i__] = k1;
      ++k1;
      goto L306;
    L305:
      indi[i__] = k2;
      ++k2;
    L306:
      ;
    }
  *mode = 1;
  optim_fmani1 (mode, n, &w[1], &d__[1], &indi[1]);
  if (nr == 0)
    {
      goto L308;
    }
  k = 1;
  i__1 = nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      h__[k] = d__[i__];
      /* L307: */
      k = k + nr1 - i__;
    }
 L308:
  if (nr == *n)
    {
      goto L400;
    }
  k = np * nr - nr1 * nr / 2 + 1;
  i__1 = *n;
  for (i__ = nr1; i__ <= i__1; ++i__)
    {
      h__[k] = d__[i__];
      /* L309: */
      k = k + np - i__;
    }
  goto L400;
  /* 
   *     verification de la definie positivite de h 
   *     permutation et factorisation 
   * 
   */
 L310:
  optim_fmc11b (&h__[1], n, &k);
  if (k >= *n)
    {
      goto L312;
    }
 L311:
  if (*iprint != 0)
    {
      Sciprintf("n2qn1 remplace le hessien initial (qui n'est, pas defini positif)/ par une diagonale positive\n");
    }
  goto L300;
 L312:
  nr = *n;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L313: */
      indi[i__] = i__;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (ibloc[i__] == 0)
	{
	  goto L314;
	}
      nc = i__;
      optim_fajc1 (n, &nc, &nr, &h__[1], &w[1], &indi[1]);
    L314:
      ;
    }
  goto L400;
  /* 
   *     verification que la diagonale est positive 
   * 
   */
 L320:
  k = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (h__[k] <= 0.)
	{
	  goto L311;
	}
      /* L321: */
      k = k + np - i__;
    }
  goto L312;
  /*                 on est pret a y aller 
   */
 L400:
  indic2 = 0;
  if (*iprint < 2)
    {
      goto L410;
    }
  Sciprintf("n2qn1 bornes initiales, (0 inactive, -1 binf active, +1 bsup active)\n\t");
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      Sciprintf("\t%d %d\n", i__,ibloc[i__]);
    }
 L410:
  dnr = sqrt ((double) nr);
  acc1 = *acc * dnr;
  /* 
   *                  iteration 
   * 
   */
 L500:
  ++itr;
  if (itr != 1)
    {
      df = fa - *f;
    }
  fa = *f;
  indic1 = 0;
  /* L501: */
  if (itr <= *niter)
    {
      goto L502;
    }
  *mode = 4;
  goto L900;
 L502:
  if (*iprint <= 2)
    {
      goto L503;
    }
  Sciprintf( "n2qn1, %d iters,%d simuls, f=%15.7f\n",
	     itr,
	     nfun,
	     (*f));
 L503:
  ++iecri;
  if (iecri != -(*iprint))
    {
      goto L510;
    }
  iecri = 0;
  indic = 1;
  (*simul) (&indic, n, &x[1], f, &g[1], optim_data);
  /*              calcul de la direction de recherche 
   *              et du test d arret 
   */
 L510:
  if (nr != 0)
    {
      goto L511;
    }
  indic2 = 1;
  goto L540;
 L511:
  *mode = 1;
  optim_fmani1 (mode, n, &ga[1], &w[1], &indi[1]);
  wii = 0.;
  i__1 = nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      wi = w[i__];
      wiii = wi * scale[i__];
      wii += wiii * wiii;
      /* L512: */
      w[i__] = -wi;
    }
  wii = sqrt (wii);
  if (wii > acc1)
    {
      goto L513;
    }
  indic2 = 1;
  goto L540;
 L513:
  optim_fmc11e (&h__[1], &nr, &w[1], &w1[1], &nr);
  if (nr == *n)
    {
      goto L520;
    }
  nrp1 = nr + 1;
  i__1 = *n;
  for (i__ = nrp1; i__ <= i__1; ++i__)
    {
      /* L514: */
      w[i__] = 0.;
    }
  /*          calcul de la derivee directionnelle 
   */
 L520:
  *mode = -1;
  optim_fmani1 (mode, n, &w[1], &d__[1], &indi[1]);
  dga = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L521: */
      dga += ga[i__] * d__[i__];
    }
  if (dga < 0.)
    {
      goto L522;
    }
  indic2 = 1;
  goto L540;
 L522:
  if (indic1 == 1)
    {
      goto L550;
    }
  /*              contrainte sortante 
   */
 L540:
  optim_fcomp1 (&indic2, &ibloc[1], &indi[1], &h__[1], &ga[1], &w[1], &d__[1],
		&g[1], n, &nr, &ncs, &dga, &delta, &prop, acc, &scale[1]);
  if (ncs != 0)
    {
      goto L543;
    }
  if (indic2 != 1)
    {
      goto L541;
    }
  *mode = 1;
  goto L900;
 L541:
  *mode = -1;
  optim_fmani1 (mode, n, &w[1], &d__[1], &indi[1]);
  goto L550;
 L543:
  if (*iprint < 2)
    {
      goto L544;
    }
  Sciprintf("n2qn1, %d iters, %d simuls, f=%15.7f, borne %d desactivee\n",
	    itr,
	    nfun,
	    (*f),
	    ncs);
 L544:
  indic1 = 1;
  logic = 6;
  /*               mise a jour de ibloc et de h 
   */
  ibloc[ncs] = 0;
  optim_fretc1 (mode, n, &ncs, &nr, &h__[1], &w[1], &indi[1], &indic2);
  indic2 = 0;
  dnr = sqrt ((double) nr);
  acc1 = *acc * dnr;
  if (*mode == 0)
    {
      goto L511;
    }
  *mode = 7;
  if (*iprint != 0)
    {
      Sciprintf( "n2qn1 erreur dans la mise a jour de l\n");
    }
  goto L900;
  /*               calcul de romax 
   */
 L550:
  romax = 1e20;
  nca = 0;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      di = d__[i__];
      if (di == 0.)
	{
	  goto L555;
	}
      if (di > 0.)
	{
	  goto L552;
	}
      bi = binf[i__];
      xi = bi - x[i__];
      if (-1. >= di)
	{
	  goto L551;
	}
      if (xi <= di * 1e20)
	{
	  goto L555;
	}
    L551:
      rocand = xi / di;
      i1 = -1;
      goto L554;
    L552:
      bi = bsup[i__];
      xi = bi - x[i__];
      if (di >= 1.)
	{
	  goto L553;
	}
      if (xi > di * 1e20)
	{
	  goto L555;
	}
    L553:
      rocand = xi / di;
      i1 = 1;
    L554:
      if (rocand > romax)
	{
	  goto L555;
	}
      nca = i__;
      romax = rocand;
      isign = i1;
    L555:
      ;
    }
  /*               romax est-il nul? 
   */
  if (nca == 0)
    {
      goto L570;
    }
  if ((d__1 = romax * d__[nca], Abs (d__1)) <= scale[nca])
    {
      goto L560;
    }
  goto L570;
  /*              addition d'une contrainte 
   */
 L560:
  ibloc[nca] = isign;
  indic1 = 1;
  optim_fajc1 (n, &nca, &nr, &h__[1], &w[1], &indi[1]);
  if (*iprint >= 2 && isign < 0)
    {
      Sciprintf("n2qn1 %d iters, %d simuls, f=%15.7f, binf %d activee\n",
		itr,
		nfun,
		(*f),
		nca);
    }
  if (*iprint >= 2 && isign > 0)
    {
      Sciprintf("n2qn1, %d iters, %d simuls, f=%15.7f, bsup %d  activee\n",
		itr,
		nfun,
		(*f),
		nca);
    }
  dnr = sqrt ((double) nr);
  acc1 = *acc * dnr;
  goto L510;
  /*               recherche lineaire 
   */
 L570:
  if (itr <= *n && itr != 1 && *mode == 1)
    {
      goto L571;
    }
  ro = 1.;
  goto L573;
 L571:
  if (logic == 1)
    {
      goto L573;
    }
  if (logic != 6)
    {
      goto L572;
    }
  ro = 1.;
  goto L573;
 L572:
  ro = df * -2. / dga;
 L573:
  roa = ro;
  ro = Min (ro, romax);
  romin = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      z__ = d__[i__];
      /* L574: */
      /*Computing MAX 
       */
      d__2 = romin, d__3 = (d__1 = z__ / scale[i__], Abs (d__1));
      romin = Max (d__2, d__3);
    }
  romin = 1. / romin;
  optim_nlis0 (n, simul, optim_fuclid, &x[1], f, &dga, &ro, &romin,
	       &romax, &d__[1], &g[1], &alfa, &beta, iprint, lp, &logic,
	       &nfun, nsim, &w[1], optim_data);
  if (*iprint > 3)
    {
      Sciprintf("()\n");
    }
  if (logic <= 1)
    {
      goto L575;
    }
  if (logic == 6)
    {
      *mode = 6;
    }
  if (logic == 4)
    {
      *mode = 5;
    }
  if (logic == 5)
    {
      *mode = 0;
    }
  if (logic == 7)
    {
      *mode = indic;
    }
  goto L900;
  /*              formule de bfgs 
   */
 L575:
  theta = 1.;
  if (logic == 0)
    {
      goto L580;
    }
  dgaa = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L576: */
      dgaa += g[i__] * d__[i__];
    }
  if (dgaa < alfa * dga)
    {
      theta = alfa * dga / dgaa;
    }
 L580:
  *mode = 1;
  optim_fmani1 (mode, n, &d__[1], &w[1], &indi[1]);
  ir = -nr;
  optim_fmani1 (mode, n, &ga[1], &d__[1], &indi[1]);
  i__1 = nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L581: */
      d__[i__] = -d__[i__];
    }
  optim_fmlag1 (n, &nr, &h__[1], &w[1], &d__[1]);
  dga = 0.;
  i__1 = nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L582: */
      dga -= w[i__] * d__[i__];
    }
  d__1 = 1. / dga;
  optim_fmc11z (&h__[1], n, &nr, &d__[1], &d__1, &w1[1], &ir, &c__1, &c_b109);
  ir = -ir;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      gi = g[i__];
      g[i__] = theta * gi - ga[i__];
      /* L583: */
      ga[i__] = gi;
    }
  optim_fmani1 (mode, n, &g[1], &d__[1], &indi[1]);
  dga = 0.;
  i__1 = nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L584: */
      dga += w[i__] * d__[i__];
    }
  dga *= ro;
  ro = roa;
  d__1 = 1. / dga;
  optim_fmc11z (&h__[1], n, &nr, &d__[1], &d__1, &w1[1], &ir, &c__0, &c_b109);
  /*               test du rang de la nouvelle 
   *               sous-matrice active 
   */
  if (ir >= nr)
    {
      goto L500;
    }
  *mode = 3;
  if (*iprint == 0)
    {
      goto L900;
    }
  Sciprintf( "n2qn1 probleme dans bfgs\n");
  /*               ici,tout est termine 
   */
 L900:
  if (*mode != 5 && *mode != 3 && *mode >= 0)
    {
      goto L910;
    }
  indic = 4;
  (*simul) (&indic, n, &x[1], f, &ga[1], optim_data);
 L910:
  iz[1] = nr;
  /*          calcul de la precision obtenue 
   */
  *acc = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (ibloc[i__] != 0)
	{
	  goto L920;
	}
      gi = ga[i__];
      *acc += gi * gi;
    L920:
      ;
    }
  if (dnr == 0.)
    {
      goto L921;
    }
  *acc = sqrt (*acc) / dnr;
 L921:
  *niter = itr;
  *nsim = nfun;
 L999:
  return 0;
}				/* n2qn1a_ */




/*    Copyright INRIA 
 * 
 */

static int optim_fajc1 (int *n, int *nc, int *nr, double *h__, double *w, int *indi)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1;

  /* Local variables */
  int incm1;
  double a, b, c__;
  int i__, j, k;
  double u, v;
  int nkkmj;
  double h1, h2;
  int nsaut;
  double ai, di;
  int ii, ij, ik, nh, nj, nk, ko, nw;
  double di1;
  int nh1, nr1, nr2, inc, nkk, nrr, inc1, kom1;
  /* Parameter adjustments */
  --indi;
  --w;
  --h__;

  /* Function Body */
  inc = indi[*nc];
  nr1 = *nr + 1;
  nr2 = *nr - 1;
  nrr = *n - *nr;
  nkk = *nr - inc;
  /* 
   *         recherche des composantes de h 
   */
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
      kom1 = ko - 1;
      i__2 = kom1;
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
      nkkmj = nkk - j;
      i__2 = nkkmj;
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
  incm1 = inc - 1;
  i__1 = incm1;
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
  return 0;
}				/* fajc1_ */




static int optim_fcomp1 (int *indic2, int *ibloc, int *indi, double *h__, double *g,
			 double *d__, double *w, double *w1, int *n, int *nr, int *ncs,
			 double *dga, double *delta, double *prop, double *acc,
			 double *scale)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1, d__2;

  /* Local variables */
  double winc;
  int i__, j, k;
  double z__;
  double am, gi;
  int nh;
  double zm;
  int nh1, ibi, inc;
  double dmu;
  int inr, nrr, inc1;
  double dmu1;

  /*    Copyright INRIA 
   * 
   */
  /* Parameter adjustments */
  --h__;
  --scale;
  --w1;
  --w;
  --d__;
  --g;
  --indi;
  --ibloc;

  /* Function Body */
  *ncs = 0;
  if (*nr == *n)
    {
      return 0;
    }
  zm = 0.;
  if (*indic2 == 1)
    {
      goto L900;
    }
  *delta = 0.;
  nh = *nr * (*nr + 1) / 2;
  nrr = *n - *nr;
  optim_fmlag1 (n, nr, &h__[1], &d__[1], &w[1]);
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ibi = ibloc[i__];
      if (ibi == 0)
	{
	  goto L500;
	}
      gi = g[i__];
      inc = indi[i__];
      inc1 = inc - 1;
      inr = inc - *nr;
      winc = w[inc];
      dmu = winc + gi;
      /*Computing MIN 
       */
      d__1 = Abs (gi), d__2 = Abs (dmu);
      am = Min (d__1, d__2);
      if (Abs (winc) * 2. >= am)
	{
	  goto L500;
	}
      if (ibi == -1 && dmu >= 0.)
	{
	  goto L500;
	}
      if (ibi == 1 && dmu <= 0.)
	{
	  goto L500;
	}
      dmu = Abs (dmu);
      if (dmu * scale[i__] <= *acc)
	{
	  goto L500;
	}
      dmu1 = dmu * dmu;
      k = inr;
      nh1 = inc1 * (*n + 1) - inc1 * inc / 2 + 1;
      z__ = h__[nh1];
      if (*nr == 0)
	{
	  goto L350;
	}
      i__2 = *nr;
      for (j = 1; j <= i__2; ++j)
	{
	  w1[j] = h__[nh + k];
	  /* L200: */
	  k += nrr;
	}
      optim_fmc11e (&h__[1], nr, &w1[1], &w1[1], nr);
      k = inr;
      i__2 = *nr;
      for (j = 1; j <= i__2; ++j)
	{
	  z__ -= w1[j] * h__[nh + k];
	  /* L250: */
	  k += nrr;
	}
    L350:
      dmu1 /= z__;
      if (dmu1 <= *delta)
	{
	  goto L500;
	}
      *delta = dmu1;
      *ncs = i__;
      zm = dmu;
    L500:
      ;
    }
  if (*ncs == 0)
    {
      return 0;
    }
  if (*delta <= -(*prop) * *dga)
    {
      *ncs = 0;
    }
  return 0;
 L900:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ibi = ibloc[i__];
      if (ibi == 0)
	{
	  goto L910;
	}
      dmu = g[i__];
      if (ibi == -1 && dmu >= 0.)
	{
	  goto L910;
	}
      if (ibi == 1 && dmu <= 0.)
	{
	  goto L910;
	}
      dmu = Abs (dmu) * scale[i__];
      if (dmu <= zm)
	{
	  goto L910;
	}
      zm = dmu;
      *ncs = i__;
    L910:
      ;
    }
  if (zm <= *acc)
    {
      *ncs = 0;
    }
  return 0;
}				/* fcomp1_ */



static int optim_fmc11z (double *a, int *n, int *nr, double *z__, 
			 double *sig, double *w,
			 int *ir, int *mk, double *eps)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int i__, j, nh, nr1;

  /* 
   */
  /* Parameter adjustments */
  --a;
  --w;
  --z__;

  /* Function Body */
  if (*nr == *n)
    {
      goto L45;
    }
  nr1 = *nr + 1;
  nh = *nr * nr1 / 2 + 1;
  if (*nr == 0)
    {
      goto L25;
    }
  i__1 = *nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i__2 = *n;
      for (j = nr1; j <= i__2; ++j)
	{
	  a[nh] += *sig * z__[i__] * z__[j];
	  /* L10: */
	  ++nh;
	}
      /* L20: */
    }
 L25:
  i__1 = *n;
  for (j = nr1; j <= i__1; ++j)
    {
      i__2 = *n;
      for (i__ = j; i__ <= i__2; ++i__)
	{
	  a[nh] += *sig * z__[i__] * z__[j];
	  /* L30: */
	  ++nh;
	}
      /* L40: */
    }
  if (*nr == 0)
    {
      return 0;
    }
 L45:
  optim_fmc11a (&a[1], nr, &z__[1], sig, &w[1], ir, mk, eps);
  return 0;
}				/* fmc11z_ */




static int optim_fmc11a (double *a, int *n, double *z__, double *sig, double *w, int *ir,
			 int *mk, double *eps)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1;

  /* Local variables */
  double b;
  int i__, j;
  double r__, v, y, al;
  int ij;
  double gm;
  int ip, mm;
  double ti;
  int np;
  double tim;

  /*  update factors given in a by   sig*z*ztranspose 
   */
  /* Parameter adjustments */
  --a;
  --w;
  --z__;

  /* Function Body */
  if (*n > 1)
    {
      goto L1;
    }
  /*Computing 2nd power 
   */
  d__1 = z__[1];
  a[1] += *sig * (d__1 * d__1);
  *ir = 1;
  if (a[1] > 0.)
    {
      return 0;
    }
  a[1] = 0.;
  *ir = 0;
  return 0;
 L1:
  np = *n + 1;
  if (*sig > 0.)
    {
      goto L40;
    }
  if (*sig == 0. || *ir == 0)
    {
      return 0;
    }
  ti = 1. / *sig;
  ij = 1;
  if (*mk == 0)
    {
      goto L10;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (a[ij] != 0.)
	{
	  /*Computing 2nd power 
	   */
	  d__1 = w[i__];
	  ti += d__1 * d__1 / a[ij];
	}
      /* L7: */
      ij = ij + np - i__;
    }
  goto L20;
 L10:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L11: */
      w[i__] = z__[i__];
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ip = i__ + 1;
      v = w[i__];
      if (a[ij] > 0.)
	{
	  goto L12;
	}
      w[i__] = 0.;
      ij = ij + np - i__;
      goto L15;
    L12:
      /*Computing 2nd power 
       */
      d__1 = v;
      ti += d__1 * d__1 / a[ij];
      if (i__ == *n)
	{
	  goto L14;
	}
      i__2 = *n;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ij;
	  /* L13: */
	  w[j] -= v * a[ij];
	}
    L14:
      ++ij;
    L15:
      ;
    }
 L20:
  if (*ir <= 0)
    {
      goto L21;
    }
  if (ti > 0.)
    {
      goto L22;
    }
  if (*mk - 1 <= 0)
    {
      goto L40;
    }
  else
    {
      goto L23;
    }
 L21:
  ti = 0.;
  *ir = -(*ir) - 1;
  goto L23;
 L22:
  ti = *eps / *sig;
  if (*eps == 0.)
    {
      --(*ir);
    }
 L23:
  mm = 1;
  tim = ti;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      j = np - i__;
      ij -= i__;
      if (a[ij] != 0.)
	{
	  /*Computing 2nd power 
	   */
	  d__1 = w[j];
	  tim = ti - d__1 * d__1 / a[ij];
	}
      w[j] = ti;
      /* L30: */
      ti = tim;
    }
  goto L41;
 L40:
  mm = 0;
  tim = 1. / *sig;
 L41:
  ij = 1;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ip = i__ + 1;
      v = z__[i__];
      if (a[ij] > 0.)
	{
	  goto L53;
	}
      if (*ir > 0 || *sig < 0. || v == 0.)
	{
	  goto L52;
	}
      *ir = 1 - *ir;
      /*Computing 2nd power 
       */
      d__1 = v;
      a[ij] = d__1 * d__1 / tim;
      if (i__ == *n)
	{
	  return 0;
	}
      i__2 = *n;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ij;
	  /* L51: */
	  a[ij] = z__[j] / v;
	}
      return 0;
    L52:
      ti = tim;
      ij = ij + np - i__;
      goto L66;
    L53:
      al = v / a[ij];
      if (mm <= 0)
	{
	  goto L54;
	}
      else
	{
	  goto L55;
	}
    L54:
      ti = tim + v * al;
      goto L56;
    L55:
      ti = w[i__];
    L56:
      r__ = ti / tim;
      a[ij] *= r__;
      if (r__ == 0.)
	{
	  goto L70;
	}
      if (i__ == *n)
	{
	  goto L70;
	}
      b = al / ti;
      if (r__ > 4.)
	{
	  goto L62;
	}
      i__2 = *n;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ij;
	  z__[j] -= v * a[ij];
	  /* L61: */
	  a[ij] += b * z__[j];
	}
      goto L64;
    L62:
      gm = tim / ti;
      i__2 = *n;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ij;
	  y = a[ij];
	  a[ij] = b * z__[j] + y * gm;
	  /* L63: */
	  z__[j] -= v * y;
	}
    L64:
      tim = ti;
      ++ij;
    L66:
      ;
    }
 L70:
  if (*ir < 0)
    {
      *ir = -(*ir);
    }
  return 0;
}

static int optim_fmc11b (double *a, int *n, int *ir)
{
  /* System generated locals */
  int i__1, i__2, i__3;

  /* Local variables */
  int i__;
  double v, aa;
  int ii, ij, ik, jk, ni, ip, np;

  /*  factorize a matrix given in a 
   */
  /* Parameter adjustments */
  --a;

  /* Function Body */
  *ir = *n;
  if (*n > 1)
    {
      goto L100;
    }
  if (a[1] > 0.)
    {
      return 0;
    }
  a[1] = 0.;
  *ir = 0;
  return 0;
 L100:
  np = *n + 1;
  ii = 1;
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      aa = a[ii];
      ni = ii + np - i__;
      if (aa > 0.)
	{
	  goto L101;
	}
      a[ii] = 0.;
      --(*ir);
      ii = ni + 1;
      goto L104;
    L101:
      ip = ii + 1;
      ii = ni + 1;
      jk = ii;
      i__2 = ni;
      for (ij = ip; ij <= i__2; ++ij)
	{
	  v = a[ij] / aa;
	  i__3 = ni;
	  for (ik = ij; ik <= i__3; ++ik)
	    {
	      a[jk] -= a[ik] * v;
	      /* L102: */
	      ++jk;
	    }
	  /* L103: */
	  a[ij] = v;
	}
    L104:
      ;
    }
  if (a[ii] > 0.)
    {
      return 0;
    }
  a[ii] = 0.;
  --(*ir);
  return 0;
}				/* fmc11b_ */




static int optim_fmc11e (double *a, int *n, double *z__, double *w, int *ir)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int i__, j;
  double v;
  int i1, ii, ij, ip, np, nip;

  /*  multiply a vector z by the inverse of the factors given in a 
   */
  /* Parameter adjustments */
  --a;
  --w;
  --z__;

  /* Function Body */
  if (*ir < *n)
    {
      return 0;
    }
  w[1] = z__[1];
  if (*n > 1)
    {
      goto L400;
    }
  z__[1] /= a[1];
  return 0;
 L400:
  i__1 = *n;
  for (i__ = 2; i__ <= i__1; ++i__)
    {
      ij = i__;
      i1 = i__ - 1;
      v = z__[i__];
      i__2 = i1;
      for (j = 1; j <= i__2; ++j)
	{
	  v -= a[ij] * z__[j];
	  /* L401: */
	  ij = ij + *n - j;
	}
      w[i__] = v;
      /* L402: */
      z__[i__] = v;
    }
  z__[*n] /= a[ij];
  np = *n + 1;
  i__1 = *n;
  for (nip = 2; nip <= i__1; ++nip)
    {
      i__ = np - nip;
      ii = ij - nip;
      v = z__[i__] / a[ii];
      ip = i__ + 1;
      ij = ii;
      i__2 = *n;
      for (j = ip; j <= i__2; ++j)
	{
	  ++ii;
	  /* L410: */
	  v -= a[ii] * z__[j];
	}
      /* L411: */
      z__[i__] = v;
    }
  return 0;
}				/* fmc11e_ */



static int optim_fmani1 (int *mode, int *n, double *d__, double *w, int *indi)
{
  /* System generated locals */
  int i__1;

  /* Local variables */
  int i__;

  /*    Copyright INRIA 
   * 
   */
  /* Parameter adjustments */
  --indi;
  --w;
  --d__;

  /* Function Body */
  if (*mode == -1)
    {
      goto L20;
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L10: */
      w[indi[i__]] = d__[i__];
    }
  return 0;
 L20:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L30: */
      w[i__] = d__[indi[i__]];
    }
  return 0;
}				/* fmani1_ */



static int optim_fmlag1 (int *n, int *nr, double *a, double *z__, double *w)
{
  /* System generated locals */
  int i__1, i__2;

  /* Local variables */
  int i__, j;
  double u;
  int nh, nj, nh1, nr1, nrr;

  /* 
   */
  /* Parameter adjustments */
  --w;
  --z__;
  --a;

  /* Function Body */
  if (*nr == *n)
    {
      return 0;
    }
  nr1 = *nr + 1;
  if (*nr != 0)
    {
      goto L20;
    }
  i__1 = *n;
  for (i__ = nr1; i__ <= i__1; ++i__)
    {
      /* L10: */
      w[i__] = 0.;
    }
  return 0;
 L20:
  nrr = *n - *nr;
  nh1 = *nr * nr1 / 2;
  nh = nh1 + 1;
  i__1 = *n;
  for (j = nr1; j <= i__1; ++j)
    {
      u = 0.;
      nj = nh;
      i__2 = *nr;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  u += a[nj] * z__[i__];
	  /* L40: */
	  nj += nrr;
	}
      ++nh;
      w[j] = u;
      /* L30: */
    }
  return 0;
}				/* fmlag1_ */




static int optim_fretc1 (int *mode, int *n, int *nc, int *nr, double *h__, double *w,
			 int *indi, int *indic2)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1;

  /* Local variables */
  int nr1p1, i__, j;
  double v;
  int incmr, nsaut, i1, ii, ij, nh, nl;
  double wi;
  int nw, nh1, nr1, nr2, inc;
  double hij;
  int nii, nrr, nrm1, nmr1;

  /*    Copyright INRIA 
   * 
   */
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
  incmr = inc - nr1;
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
  i__1 = incmr;
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
  i__1 = incmr;
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
      nrm1 = *n - nr1;
      i__2 = nrm1;
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
  incmr = inc - nr1;
  i__1 = incmr;
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
  if (*indic2 != 1)
    {
      goto L190;
    }
  i__1 = *nr;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L185: */
      w[i__] = 0.;
    }
  if (*n == nr1)
    {
      goto L190;
    }
  nr1p1 = nr1 + 1;
  i__1 = *n;
  for (i__ = nr1p1; i__ <= i__1; ++i__)
    {
      /* L187: */
      w[i__] = 0.;
    }
  /*         stockage de w dans h 
   */
 L190:
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
  nmr1 = *n - nr1;
  i__1 = nmr1;
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
}				/* fretc1_ */
