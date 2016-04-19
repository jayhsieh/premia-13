#include "optim.h"

/*    Copyright INRIA 
 * 
 *    mise a jour des vecteurs ({y}(i),{s}(i),ys(i),i=1,np) 
 * 
 *    -----mise a jour de y(lb, ) , s(lb, ) , ys(lb) 
 */


int
optim_majysa (int *n, int *nt, int *np, double *y, double *s, double *ys,
	      int *lb, double *g, double *x, double *g1, double *x1,
	      int *index, int *ialg, int *nb)
{
  /* System generated locals */
  int y_dim1, y_offset, s_dim1, s_offset, i__1;

  /* Local variables */
  int i__, ij;

  /* Parameter adjustments */
  --index;
  --x1;
  --g1;
  --x;
  --g;
  --ys;
  s_dim1 = *nt;
  s_offset = s_dim1 + 1;
  s -= s_offset;
  y_dim1 = *nt;
  y_offset = y_dim1 + 1;
  y -= y_offset;
  --ialg;

  /* Function Body */
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      y[*lb + i__ * y_dim1] = g[i__] - g1[i__];
      s[*lb + i__ * s_dim1] = x[i__] - x1[i__];
      /* L100: */
    }
  ys[*lb] = 0.;
  i__1 = *n;
  for (i__ = 1; i__ <= i__1; ++i__)
    {
      ys[*lb] += y[*lb + i__ * y_dim1] * s[*lb + i__ * s_dim1];
      /* L200: */
    }
  /* 
   *    accumulation eventuelle 
   */
  if (ialg[8] == 5 && *np > 0)
    {
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  y[i__ * y_dim1 + 1] += y[*lb + i__ * y_dim1];
	  s[i__ * s_dim1 + 1] += s[*lb + i__ * s_dim1];
	  /* L20: */
	}
      ys[1] = 0.;
      i__1 = *n;
      for (i__ = 1; i__ <= i__1; ++i__)
	{
	  /* L30: */
	  ys[1] += y[i__ * y_dim1 + 1] * s[i__ * s_dim1 + 1];
	}
    }
  /* 
   * 
   *    -----mise a jour de np et index 
   */
  if (*np < *nt)
    {
      ++(*np);
      index[*lb] = *np;
    }
  else
    {
      ij = *lb;
      i__1 = *nt;
      for (i__ = *nb; i__ <= i__1; ++i__)
	{
	  ++ij;
	  if (ij > *nt)
	    {
	      ij = *nb;
	    }
	  index[i__] = ij;
	  /* L300: */
	}
    }
  /* 
   *    ------chercher la prochaine place libre 
   */
  if (*lb == *nt)
    {
      *lb = *nb;
    }
  else
    {
      ++(*lb);
    }
  /* 
   *    -------------- 
   */
  return 0;
}				/* majysa_ */
