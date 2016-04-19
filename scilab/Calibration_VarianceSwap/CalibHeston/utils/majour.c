#include "optim.h"

int
optim_majour (double *hm, double *hd, double *dd, int *n, double *hno,
	      int *ir, int *indic, double *eps)
{
  /* System generated locals */
  int i__1, i__2;
  double d__1;

  /* Local variables */
  double honm, b;
  int i__, j;
  double r__, y;
  int iplus;
  double gm;
  int ll, mm, np;
  double del, hml, hon;

  /*    Copyright INRIA 
   */
  /* Parameter adjustments */
  --hm;
  --dd;
  --hd;

  /* Function Body */
  if (*n == 1)
    {
      goto L100;
    }
  /* 
   */
  np = *n + 1;
  if (*hno > 0.)
    {
      goto L99;
    }
  /* 
   */
  if (*hno == 0.)
    {
      goto L999;
    }
  if (*ir == 0)
    {
      goto L999;
    }
  hon = 1. / *hno;
  ll = 1;
  if (*indic == 0)
    {
      goto L1;
    }
  /* 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      if (hm[ll] == 0.)
	{
	  goto L2;
	}
      /*Computing 2nd power 
       */
      d__1 = dd[i__];
      hon += d__1 * d__1 / hm[ll];
    L2:
      ll = ll + np - i__;
    }
  goto L3;
  /* 
   */
 L1:
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      dd[i__] = hd[i__];
      /* L4: */
    }
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      iplus = i__ + 1;
      del = dd[i__];
      if (hm[ll] > 0.)
	{
	  goto L6;
	}
      dd[i__] = 0.;
      ll = ll + np - i__;
      goto L5;
    L6:
      /*Computing 2nd power 
       */
      d__1 = del;
      hon += d__1 * d__1 / hm[ll];
      if (i__ == *n)
	{
	  goto L7;
	}
      i__2 = *n;
      for (j = iplus; j <= i__2; ++j)
	{
	  ++ll;
	  /* L8: */
	  dd[j] -= del * hm[ll];
	}
    L7:
      ++ll;
    L5:
      ;
    }
  /* 
   */
 L3:
  if (*ir <= 0)
    {
      goto L9;
    }
  if (hon > 0.)
    {
      goto L10;
    }
  if (*indic - 1 <= 0)
    {
      goto L99;
    }
  else
    {
      goto L11;
    }
 L9:
  hon = 0.;
  *ir = -(*ir) - 1;
  goto L11;
 L10:
  hon = *eps / *hno;
  if (*eps == 0.)
    {
      --(*ir);
    }
 L11:
  mm = 1;
  honm = hon;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      j = np - i__;
      ll -= i__;
      if (hm[ll] != 0.)
	{
	  /*Computing 2nd power 
	   */
	  d__1 = dd[j];
	  honm = hon - d__1 * d__1 / hm[ll];
	}
      dd[j] = hon;
      /* L12: */
      hon = honm;
    }
  goto L13;
  /* 
   */
 L99:
  mm = 0;
  honm = 1. / *hno;
 L13:
  ll = 1;
  /* 
   */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      iplus = i__ + 1;
      del = hd[i__];
      if (hm[ll] > 0.)
	{
	  goto L14;
	}
      if (*ir > 0)
	{
	  goto L15;
	}
      if (*hno < 0.)
	{
	  goto L15;
	}
      if (del == 0.)
	{
	  goto L15;
	}
      *ir = 1 - *ir;
      /*Computing 2nd power 
       */
      d__1 = del;
      hm[ll] = d__1 * d__1 / honm;
      if (i__ == *n)
	{
	  goto L999;
	}
      i__2 = *n;
      for (j = iplus; j <= i__2; ++j)
	{
	  ++ll;
	  /* L16: */
	  hm[ll] = hd[j] / del;
	}
      goto L999;
    L15:
      hon = honm;
      ll = ll + np - i__;
      goto L98;
    L14:
      hml = del / hm[ll];
      if (mm <= 0)
	{
	  goto L17;
	}
      else
	{
	  goto L18;
	}
    L17:
      hon = honm + del * hml;
      goto L19;
    L18:
      hon = dd[i__];
    L19:
      r__ = hon / honm;
      hm[ll] *= r__;
      if (r__ == 0.)
	{
	  goto L20;
	}
      if (i__ == *n)
	{
	  goto L20;
	}
      b = hml / hon;
      if (r__ > 4.)
	{
	  goto L21;
	}
      i__2 = *n;
      for (j = iplus; j <= i__2; ++j)
	{
	  ++ll;
	  hd[j] -= del * hm[ll];
	  /* L22: */
	  hm[ll] += b * hd[j];
	}
      goto L23;
    L21:
      gm = honm / hon;
      i__2 = *n;
      for (j = iplus; j <= i__2; ++j)
	{
	  ++ll;
	  y = hm[ll];
	  hm[ll] = b * hd[j] + y * gm;
	  /* L24: */
	  hd[j] -= del * y;
	}
    L23:
      honm = hon;
      ++ll;
    L98:
      ;
    }
  /* 
   */
 L20:
  if (*ir < 0)
    {
      *ir = -(*ir);
    }
  goto L999;
 L100:
  /*Computing 2nd power 
   */
  d__1 = hd[1];
  hm[1] += *hno * (d__1 * d__1);
  *ir = 1;
  if (hm[1] > 0.)
    {
      goto L999;
    }
  hm[1] = 0.;
  *ir = 0;
 L999:
  return 0;
}				/* majour_ */
