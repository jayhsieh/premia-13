#include "optim.h"


/* 
 *    calcul de l'etat initial dans ICSE : cas standard : 
 *    controle par l'etat initial 
 *    Copyright INRIA 
 *! 
 * 
 */

int
optim_icsei (int *indi, int *nui, double *u, double *y0, double *y0u,
	     int *itu, double *dtu, double *t0, double *tf, double *dti,
	     double *dtf, double *ermx, int *iu, int *nuc, int *nuv,
	     int *ilin, int *nti, int *ntf, int *ny, int *nea, int *itmx,
	     int *nex, int *nob, int *ntob, int *ntobi, int *nitu, int *ndtu)
{
  double c_b3 = 0.;
  int c__1 = 1;
  /* System generated locals */
  int y0u_dim1, y0u_offset, i__1;

  /* Local variables */
  int i__;

  /* Parameter adjustments */
  --u;
  --iu;
  y0u_dim1 = *ny;
  y0u_offset = y0u_dim1 + 1;
  y0u -= y0u_offset;
  --y0;
  --itu;
  --dtu;

  /* Function Body */
  if (*indi == 1)
    {
      i__1 = *ny;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  y0[i__] = u[i__];
	}
    }
  /* 
   */
  if (*indi == 2)
    {
      /*      cas ou y0u est l identite 
       */
      i__1 = *ny * *nui;
      nsp_dset (&i__1, &c_b3, &y0u[y0u_offset], &c__1);
      i__1 = *ny;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  y0u[i__ + i__ * y0u_dim1] = 1.;
	}
    }
  return 0;
}
