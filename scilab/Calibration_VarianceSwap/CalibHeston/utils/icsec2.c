#include "optim.h"




/* 
 * 
 *    critere standard des moindres carres 
 * 
 * 
 *     Cout ponctuel : 
 *      Parametres d'entree : 
 *       indc     : 1 si on desire calculer c2,2 si on desire 
 *                  calculer c2y,c2u 
 *       tob      : instants de mesure 
 *       obs      : matrice d'observation 
 *       cof      : coefficients de ponderation du cout 
 *       ytob     : valeur de l'etat aux instants d'observation 
 *       ob       : mesures 
 *       u(nu)    : controle.Le controle variable est stocke a la 
 *                  suite du controle suite du constant. 
 *      Parametres de sortie : 
 *       indc     : comme pour icsec1 
 *       c2       : cout 
 *       c2y(ny,ntob) : derivee de c2 par rapport a y 
 *       g(nu)  : derivee de c2 par rapport a u 
 *    Copyright INRIA 
 *! 
 * 
 *    critere standard des moindres carres 
 * 
 */


int optim_icsec2 (int *indc, int *nu, double *tob, double *obs, double *cof,
		  double *ytob, double *ob, double *u, double *c__, double *cy,
		  double *g, double *yob, double *d__, int *itu, double *dtu,
		  double *t0, double *tf, double *dti, double *dtf, double *ermx,
		  int *iu, int *nuc, int *nuv, int *ilin, int *nti, int *ntf,
		  int *ny, int *nea, int *itmx, int *nex, int *nob, int *ntob,
		  int *ntobi, int *nitu, int *ndtu)
{
  int c__1 = 1;
  /* System generated locals */
  int obs_dim1, obs_offset, cof_dim1, cof_offset, ytob_dim1, ytob_offset,
    ob_dim1, ob_dim2, ob_offset, cy_dim1, cy_offset, yob_dim1, yob_offset,
    i__1, i__2, i__3;
  double d__1;

  /* Local variables */
  int i__, j, k;

  /* Parameter adjustments */
  --g;
  --u;
  --iu;
  --d__;
  obs_dim1 = *nob;
  obs_offset = obs_dim1 + 1;
  obs -= obs_offset;
  yob_dim1 = *nob;
  yob_offset = yob_dim1 + 1;
  yob -= yob_offset;
  cy_dim1 = *ny;
  cy_offset = cy_dim1 + 1;
  cy -= cy_offset;
  ob_dim1 = *nex;
  ob_dim2 = *ntob;
  ob_offset = ob_dim1 * (ob_dim2 + 1) + 1;
  ob -= ob_offset;
  ytob_dim1 = *ny;
  ytob_offset = ytob_dim1 + 1;
  ytob -= ytob_offset;
  cof_dim1 = *nob;
  cof_offset = cof_dim1 + 1;
  cof -= cof_offset;
  --tob;
  --itu;
  --dtu;

  /* Function Body */
  dmmul_scicos (&obs[obs_offset], nob, &ytob[ytob_offset], ny,
		&yob[yob_offset], nob, nob, ny, ntob);
  if (*indc == 1)
    {
      *c__ = 0.;
      i__1 = *nob;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  i__2 = *ntob;
	  for (j = 1; j <= i__2; ++j)
	    {
	      i__3 = *nex;
	      for (k = 1; k <= i__3; ++k)
		{
		  /*Computing 2nd power 
		   */
		  d__1 =
		    yob[i__ + j * yob_dim1] - ob[k +
						 (j +
						  i__ * ob_dim2) * ob_dim1];
		  *c__ += cof[i__ + j * cof_dim1] * .5 * (d__1 * d__1);
		  /* L10: */
		}
	      /* L11: */
	    }
	  /* L12: */
	}
    }
  else
    {
      i__1 = *ntob;
      for (j = 1; j <= i__1; ++j)
	{
	  i__2 = *nob;
	  for (i__ = 1; i__ <= i__2; ++i__)
	    {
	      d__[i__] = 0.;
	      i__3 = *nex;
	      for (k = 1; k <= i__3; ++k)
		{
		  d__[i__] +=
		    cof[i__ + j * cof_dim1] * (yob[i__ + j * yob_dim1] -
					       ob[k +
						  (j +
						   i__ * ob_dim2) * ob_dim1]);
		  /* L24: */
		}
	      /* L25: */
	    }
	  dmmul_scicos (&d__[1], &c__1, &obs[obs_offset], nob,
			&cy[j * cy_dim1 + 1], &c__1, &c__1, nob, ny);
	  /* L20: */
	}
    }
  return 0;
}				/* icsec2_ */
