#include "optim.h"
/*    Copyright INRIA 
 *---- 
 * 
 *    nlis0 + minuscules + commentaires 
 *    --------------------------------- 
 * 
 *       en sortie logic = 
 * 
 *       0          descente serieuse 
 *       1          descente bloquee 
 *       4          nap > napmax 
 *       5          retour a l'utilisateur 
 *       6          fonction et gradient pas d'accord 
 *       < 0        contrainte implicite active 
 * 
 *---- 
 * 
 *--- arguments 
 * 
 * 
 *--- variables locales 
 * 
 */

int optim_nlis0 (int *n, opt_simul simul, opt_prosca prosca, double *xn, double *fn,
		 double *fpn, double *t, double *tmin, double *tmax, double *d__,
		 double *g, double *amd, double *amf, int *imp, int *io,
		 int *logic, int *nap, int *napmax, double *x, opt_simul_data *optim_data)
{
  /* System generated locals */
  int i__1;
  double d__1, d__2, d__3, d__4;

  /* Local variables */
  double tesd, tesf, test, f;
  int i__, indic;
  double z__, d2, z1, fa, fd, fg, ta, fp;
  int indica;
  double td;
  int indicd;
  double tg, fpa, ffn, fpd, fpg;

  /* Parameter adjustments */
  --x;
  --g;
  --d__;
  --xn;

  if (*n > 0 && *fpn < 0. && *t > 0. && *tmax > 0. && *amf > 0. && *amd > *amf
      && *amd < 1.)
    {
      goto L5;
    }
  *logic = 6;
  goto L999;
 L5:
  tesf = *amf * *fpn;
  tesd = *amd * *fpn;
  td = 0.;
  tg = 0.;
  fg = *fn;
  fpg = *fpn;
  ta = 0.;
  fa = *fn;
  fpa = *fpn;
  (*prosca) (n, &d__[1], &d__[1], &d2, optim_data);
  /* 
   *              elimination d'un t initial ridiculement petit 
   * 
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
      Sciprintf("nlis0: tmin force a tmax\n");
    }
  *tmin = *tmax;
 L20:
  if (*fn + *t * *fpn < *fn + *t * .9 * *fpn)
    {
      goto L30;
    }
  *t *= 2.;
  goto L20;
 L30:
  indica = 1;
  *logic = 0;
  if (*t > *tmax)
    {
      *t = *tmax;
      *logic = 1;
    }
  if (*imp >= 4)
    { 
      Sciprintf("nlis0: fpn=%10.3g, d2=%9.2g, tmin=%9.2g, tmax= %9.2g\n",*fpn,d2,*tmin,*tmax);
    }
  /* 
   *    --- nouveau x 
   * 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      x[i__] = xn[i__] + *t * d__[i__];
      /* L50: */
    }
  /* 
   *--- boucle 
   * 
   */
 L100:
  ++(*nap);
  if (*nap > *napmax)
    {
      *logic = 4;
      *fn = fg;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  xn[i__] += tg * d__[i__];
	  /* L120: */
	}
      goto L999;
    }
  indic = 4;
  /* 
   *    --- appel simulateur 
   * 
   */
  (*simul) (&indic, n, &x[1], &f, &g[1], optim_data);
  if (indic == 0)
    {
      /* 
       *        --- arret demande par l'utilisateur 
       * 
       */
      *logic = 5;
      *fn = f;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  xn[i__] = x[i__];
	  /* L170: */
	}
      goto L999;
    }
  if (indic < 0)
    {
      /* 
       *        --- les calculs n'ont pas pu etre effectues par le simulateur 
       * 
       */
      td = *t;
      indicd = indic;
      *logic = 0;
      if (*imp >= 4)
	{
	  Sciprintf("nlis0: %10.3f, indic=%d\n",*t,indic);
	}
      *t = tg + (td - tg) * .1;
      goto L905;
    }
  /* 
   *    --- les tests elementaires sont faits, on y va 
   * 
   */
  (*prosca) (n, &d__[1], &g[1], &fp, optim_data);
  /* 
   *    --- premier test de Wolfe 
   * 
   */
  ffn = f - *fn;
  if (ffn > *t * tesf)
    {
      td = *t;
      fd = f;
      fpd = fp;
      indicd = indic;
      *logic = 0;
      if (*imp >= 4)
	{
	  Sciprintf("nlis0: %14.3f, %11.3f, %11.3f\n",*t,ffn,fp);
	}
      goto L500;
    }
  /* 
   *    --- test 1 ok, donc deuxieme test de Wolfe 
   * 
   */
  if (*imp >= 4)
    {
      Sciprintf("nlis0: %14.3f, %11.3f, %11.3f\n",*t,ffn,fp);
    }
  if (fp > tesd)
    {
      *logic = 0;
      goto L320;
    }
  if (*logic == 0)
    {
      goto L350;
    }
  /* 
   *    --- test 2 ok, donc pas serieux, on sort 
   * 
   */
 L320:
  *fn = f;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      xn[i__] = x[i__];
      /* L330: */
    }
  goto L999;
  /* 
   * 
   * 
   */
 L350:
  tg = *t;
  fg = f;
  fpg = fp;
  if (td != 0.)
    {
      goto L500;
    }
  /* 
   *             extrapolation 
   * 
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
   *             interpolation 
   * 
   */
 L500:
  if (indica <= 0)
    {
      ta = *t;
      *t = tg * .9 + td * .1;
      goto L900;
    }
  z__ = fp + fpa - (fa - f) * 3. / (ta - *t);
  z1 = z__ * z__ - fp * fpa;
  if (z1 < 0.)
    {
      ta = *t;
      *t = (td + tg) * .5;
      goto L900;
    }
  if (*t < ta)
    {
      z1 = z__ - sqrt (z1);
    }
  if (*t > ta)
    {
      z1 = z__ + sqrt (z1);
    }
  z__ = fp / (fp + z1);
  z__ = *t + z__ * (ta - *t);
  ta = *t;
  test = (td - tg) * .1;
  /*Computing MAX 
   */
  d__1 = z__, d__2 = tg + test;
  *t = Max (d__1, d__2);
  /*Computing MIN 
   */
  d__1 = *t, d__2 = td - test;
  *t = Min (d__1, d__2);
  /* 
   *--- fin de boucle 
   *    - t peut etre bloque sur tmax 
   *      (venant de l'extrapolation avec logic=1) 
   * 
   */
 L900:
  fa = f;
  fpa = fp;
 L905:
  indica = indic;
  /* 
   *--- faut-il continuer ? 
   * 
   */
  if (td == 0.)
    {
      goto L950;
    }
  if (td - tg < *tmin)
    {
      goto L920;
    }
  /* 
   *    --- limite de precision machine (arret de secours) ? 
   * 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      z__ = xn[i__] + *t * d__[i__];
      if (z__ != xn[i__] && z__ != x[i__])
	{
	  goto L950;
	}
      /* L910: */
    }
  /* 
   *--- arret sur dxmin ou de secours 
   * 
   */
 L920:
  *logic = 6;
  /* 
   *    si indicd<0, les derniers calculs n'ont pas pu etre fait par simul 
   * 
   */
  if (indicd < 0)
    {
      *logic = indicd;
    }
  /* 
   *    si tg=0, xn = xn_depart, 
   *    sinon on prend xn=x_gauche qui fait decroitre f 
   * 
   */
  if (tg == 0.)
    {
      goto L940;
    }
  *fn = fg;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L930: */
      xn[i__] += tg * d__[i__];
    }
 L940:
  if (*imp <= 0)
    {
      goto L999;
    }
  Sciprintf("nlis0: fin sur tmin,pas,fonctions,derivees\n");
  Sciprintf("nlis0: %18.8f, %18.8f,%11.3f\n",tg,fg,fpg);
  if (*logic == 6)
    {
      Sciprintf("nlis0: %18.8f, %18.8f,%11.3f\n",td,fd,fpd);
    }
  if (*logic == 7)
    {
      Sciprintf("nlis0: %18.8f, indic=%d\n",td,indicd);
    }
  goto L999;
  /* 
   *              recopiage de x et boucle 
   * 
   */
 L950:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      /* L960: */
      x[i__] = xn[i__] + *t * d__[i__];
    }
  goto L100;
 L999:
  return 0;
}

