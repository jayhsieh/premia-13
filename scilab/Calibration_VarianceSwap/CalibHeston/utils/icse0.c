#include "optim.h"

/* Common Block Declarations */

struct
{
  double t00, tf0, dti0, dtf0, ermx0;
  int iu0[5], nuc0, nuv0, ilin0, nti0, ntf0, ny0, nea0, itmx0, nex0, nob0,
    ntob0, ntobi0, nitu0, ndtu0;
} icsez_;

#define icsez_1 icsez_

struct
{
  int nitv0, nrtv0, ndtv0;
} nird_;

#define nird_1 nird_

int
optim_icse0 (int *nu, double *t0, double *tf, double *dti, double *dtf,
	     double *ermx, int *iu, int *nuc, int *nuv, int *ilin, int *nti,
	     int *ntf, int *ny, int *nea, int *itmx, int *nex, int *nob,
	     int *ntob, int *ntobi, int *nitu, int *ndtu, int *nitv,
	     int *nrtv, int *ndtv)
{
  int i__;
  double zz;
  int ind, izz;

  /* 
   *    programme d'initialisation appele par icse.bas 
   *    initialisation des commons icsez icsez0 et nird 
   *    Copyright INRIA 
   * 
   * 
   */
  /* Parameter adjustments */
  --iu;

  /* Function Body */
  icsez_1.t00 = *t0;
  icsez_1.tf0 = *tf;
  icsez_1.dti0 = *dti;
  icsez_1.dtf0 = *dtf;
  icsez_1.ermx0 = *ermx;
  for (i__ = 1; i__ <= 5; ++i__)
    {
      /* L10: */
      icsez_1.iu0[i__ - 1] = iu[i__];
    }
  icsez_1.nuc0 = *nuc;
  icsez_1.nuv0 = *nuv;
  icsez_1.ilin0 = *ilin;
  icsez_1.nti0 = *nti;
  icsez_1.ntf0 = *ntf;
  icsez_1.ny0 = *ny;
  icsez_1.nea0 = *nea;
  icsez_1.itmx0 = *itmx;
  icsez_1.nex0 = *nex;
  icsez_1.nob0 = *nob;
  icsez_1.ntob0 = *ntob;
  icsez_1.ntobi0 = *ntobi;
  icsez_1.nitu0 = *nitu;
  icsez_1.ndtu0 = *ndtu;
  nird_1.nitv0 = 0;
  nird_1.nrtv0 = 0;
  nird_1.ndtv0 = 0;
  ind = 0;
  optim_icse (&ind, nu, &zz, &zz, &zz, &izz, &zz, &zz, NULL, NULL, NULL);
  *nitv = Max (1, nird_1.nitv0);
  *nrtv = Max (1, nird_1.nrtv0);
  *ndtv = Max (1, nird_1.ndtv0);
  return 0;
}				/* icse0_ */
