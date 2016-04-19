#include "optim.h"

/*    Copyright INRIA 
 * 
 *          Using f and fp at t and ta, computes new t by cubic formula 
 *          safeguarded inside [tlower,tupper]. 
 * 
 */

int optim_fcube (double *t, double *f, double *fp, double *ta, double *fa,
		 double *fpa, double *tlower, double *tupper)
{
  double d__1;
  double sign, anum, b, z1, discri, den;

  z1 = *fp + *fpa - (*fa - *f) * 3. / (*ta - *t);
  b = z1 + *fp;
  /* 
   *             first compute the discriminant (without overflow) 
   */
  if (Abs (z1) <= 1.)
    {
      discri = z1 * z1 - *fp * *fpa;
    }
  else
    {
      discri = *fp / z1;
      discri *= *fpa;
      discri = z1 - discri;
      if (z1 >= 0. && discri >= 0.)
	{
	  discri = sqrt (z1) * sqrt (discri);
	  goto L200;
	}
      if (z1 <= 0. && discri <= 0.)
	{
	  discri = sqrt (-z1) * sqrt (-discri);
	  goto L200;
	}
      discri = -1.;
    }
  if (discri < 0.)
    {
      if (*fp < 0.)
	{
	  *t = *tupper;
	}
      if (*fp >= 0.)
	{
	  *t = *tlower;
	}
      goto L990;
    }
  /* 
   *      discriminant nonnegative, stable solution formula 
   */
  discri = sqrt (discri);
 L200:
  if (*t - *ta < 0.)
    {
      discri = -discri;
    }
  sign = (*t - *ta) / (d__1 = *t - *ta, Abs (d__1));
  if (b * sign > 0.)
    {
      anum = (*ta - *t) * *fp;
      den = b + discri;
    }
  else
    {
      den = z1 + b + *fpa;
      anum = (*ta - *t) * (b - discri);
    }
  /* 
   *              now compute the ratio (without overflow) 
   * 
   */
  if (Abs (den) >= 1.)
    {
      *t += anum / den;
    }
  else
    {
      if (Abs (anum) < (*tupper - *tlower) * Abs (den))
	{
	  *t += anum / den;
	}
      else
	{
	  if (*fp < 0.)
	    {
	      *t = *tupper;
	    }
	  if (*fp >= 0.)
	    {
	      *t = *tlower;
	    }
	}
    }
  /* 
   *                      finally, safeguard 
   * 
   */
  *t = Max (*t, *tlower);
  *t = Min (*t, *tupper);
 L990:
  return 0;
}
