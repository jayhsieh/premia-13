#include "optim.h"

/*
 *    Copyright INRIA 
 */

int
optim_icscof (int *ico, int *ntob, int *nex, int *nob, double *yob,
	      double *ob, double *cof)
{
  /* System generated locals */
  int yob_dim1, yob_offset, ob_dim1, ob_dim2, ob_offset, cof_dim1, cof_offset,
    i__1, i__2, i__3;
  double d__1;

  /* Local variables */
  int i__, j, k;

  /*    ce programme est appele par les macros icsua (ico=1) et icsuq 
   *    (ico=2) de icse.bas pour le calcul initial des coefficients 
   *    de ponderation du cout 
   * 
   *    en entree:(pour ico=2) 
   * 
   *    yob      double precision (nob,ntob) 
   *             yob=obs*ytob,avec obs(nob,ny) matrice d'observation et 
   *             ytob(ny,ntob) valeurs calculees de l'etat aux instants 
   *             de mesure 
   * 
   *    ob       double precision (nex,ntob,nob) 
   *             mesures 
   * 
   *    en sortie: 
   * 
   *    cof      double precision (nob,ntob) 
   *             coefficients de ponderation du cout 
   * 
   */
  /* Parameter adjustments */
  cof_dim1 = *nob;
  cof_offset = cof_dim1 + 1;
  cof -= cof_offset;
  ob_dim1 = *nex;
  ob_dim2 = *ntob;
  ob_offset = ob_dim1 * (ob_dim2 + 1) + 1;
  ob -= ob_offset;
  yob_dim1 = *nob;
  yob_offset = yob_dim1 + 1;
  yob -= yob_offset;

  /* Function Body */
  i__1 = *nob;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      i__2 = *ntob;
      for (j = 1; j <= i__2; ++j)
	{
	  /* L5: */
	  cof[i__ + j * cof_dim1] = 0.;
	}
    }
  /*    si ico=1 (macro icsua:ponderation "arithmetique" du cout) 
   *    les coefficients de ponderation du cout cof(nob,ntob) 
   *    sont:cof(i,j)=nex/(|ob(1,j,i)|+..+|ob(nex,j,i)|) 
   */
  if (*ico == 1)
    {
      i__2 = *nob;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  i__1 = *ntob;
	  for (j = 1; j <= i__1; ++j)
	    {
	      i__3 = *nex;
	      for (k = 1; k <= i__3; ++k)
		{
		  /* L10: */
		  cof[i__ + j * cof_dim1] += (d__1 =
					      ob[k +
						 (j +
						  i__ * ob_dim2) * ob_dim1],
					      Abs (d__1));
		}
	    }
	}
      i__3 = *nob;
      for (i__ = 1; i__ <= i__3; ++i__)
	{
	  i__1 = *ntob;
	  for (j = 1; j <= i__1; ++j)
	    {
	      /* L15: */
	      cof[i__ + j * cof_dim1] =
		(double) (*nex) / cof[i__ + j * cof_dim1];
	    }
	}
      /*    si ico=2 (macro icsuq:ponderation "quadratique" du cout) 
       *    les coefficients de ponderation du cout cof(nob,ntob) sont: 
       *cof(i,j)=1/2*[(yob(i,j)-ob(1,j,i))**2+..+(yob(i,j)-ob(nex,j,i))**2] 
       */
    }
  else
    {
      i__1 = *nob;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  i__3 = *ntob;
	  for (j = 1; j <= i__3; ++j)
	    {
	      i__2 = *nex;
	      for (k = 1; k <= i__2; ++k)
		{
		  /* L20: */
		  /*Computing 2nd power 
		   */
		  d__1 =
		    yob[i__ + j * yob_dim1] - ob[k +
						 (j +
						  i__ * ob_dim2) * ob_dim1];
		  cof[i__ + j * cof_dim1] += d__1 * d__1;
		}
	    }
	}
      i__2 = *nob;
      for (i__ = 1; i__ <= i__2; ++i__)
	{
	  i__3 = *ntob;
	  for (j = 1; j <= i__3; ++j)
	    {
	      /* L25: */
	      cof[i__ + j * cof_dim1] = .5 / cof[i__ + j * cof_dim1];
	    }
	}
    }
  return 0;
}				/* icscof_ */
