#include "optim.h"

#ifndef d_sign 
#define d_sign(a,b) ((*b) >= 0 ? (( (*a) >= 0 ? (*a) : - (*a))) : -(((*a) >= 0 ? (*a) : - (*a))))
#endif 

int
optim_fpq2 (int *inout, double *x, double *cx, double *fx, double *gx,
	    double *d__, double *sthalf, double *penlty, int *iyflag,
	    double *y, double *cy, double *fy, double *gy, double *z__,
	    double *cz, double *fz, double *gz, double *gg, double *hh,
	    double *s)
{
  double zero = 0.;
  double half = .5;
  double d__1, d__2;
  double absd, p, denom, absgx, smallh, dlower, dupper, gyplus, xminsy;

  absd = Abs (*d__);
  if (*inout == 0)
    {
      *iyflag = 0;
      *gg = zero;
      *hh = zero;
      *s = absd;
      if (*sthalf <= zero || *sthalf >= half)
	{
	  *sthalf = half * half;
	}
      if (*penlty <= zero)
	{
	  *penlty = half + half;
	}
      if (*gx != zero)
	{
	  *d__ = -d_sign (&absd, gx);
	}
      *inout = 1;
    }
  else
    {
      if (*cz > zero || *fz >= *fx)
	{
	  *inout = 3;
	  if (*iyflag == 0)
	    {
	      /*           if (cz.gt.zero) go to 100 
	       */
	      *gg = (*gz - *gx) / *d__;
	      *hh = *gg;
	      /* L100: */
	      *s = *sthalf / absd;
	      *iyflag = 1;
	    }
	  else
	    {
	      *hh = (*gz - *gy) / (*d__ - (*y - *x));
	    }
	  *y = *z__;
	  *cy = *cz;
	  *fy = *fz;
	  *gy = *gz;
	}
      else
	{
	  if (*gx * *gz < zero)
	    {
	      *inout = 2;
	      *hh = *gg;
	      if (*iyflag == 0)
		{
		  *gg = (*gz - *gx) / *d__;
		  *s = *sthalf / absd;
		  *iyflag = 1;
		}
	      else
		{
		  *gg = (*gz - *gy) / (*d__ - (*y - *x));
		}
	      *y = *x;
	      *cy = *cx;
	      *fy = *fx;
	      *gy = *gx;
	    }
	  else
	    {
	      *inout = 1;
	      *gg = (*gz - *gx) / *d__;
	    }
	  *x = *z__;
	  *cx = *cz;
	  *fx = *fz;
	  *gx = *gz;
	}
      if (*iyflag == 0)
	{
	  dlower = *s;
	  dupper = absd / *sthalf;
	  xminsy = -(*d__);
	}
      else
	{
	  xminsy = *x - *y;
	  smallh = Min (zero, *hh) * xminsy * half;
	  gyplus = *gy + smallh;
	  /*         if (cy.le.zero) then 
	   */
	  p = *fx - *fy - gyplus * xminsy;
	  denom = (d__1 = gyplus + smallh - *gx, Abs (d__1));
	  /*         else 
	   *           p=-cy-gyplus*xminsy 
	   *           denom=abs(gyplus+smallh-gx*p*penlty) 
	   *         end if 
	   */
	  if (p >= zero)
	    {
	      goto L500;
	    }
	  p = zero;
	  *s = *sthalf / Abs (xminsy);
	L500:
	  dlower = *s * xminsy * xminsy;
	  dupper = Abs (xminsy) - dlower;
	  if (Abs (p) < denom * dupper)
	    {
	      /*Computing MAX 
	       */
	      d__1 = dlower, d__2 = p / denom;
	      dupper = Max (d__1, d__2);
	    }
	}
      absgx = Abs (*gx);
      absd = dupper;
      if (absgx < *gg * dupper)
	{
	  /*Computing MAX 
	   */
	  d__1 = dlower, d__2 = absgx / *gg;
	  absd = Max (d__1, d__2);
	}
      *d__ = -d_sign (&absd, &xminsy);
    }
  *z__ = *x + *d__;
  return 0;
}
