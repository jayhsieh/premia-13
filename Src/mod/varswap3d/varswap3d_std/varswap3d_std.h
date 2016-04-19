#ifndef _VARSWAP3D_STD_H
#define _VARSWAP3D_STD_H

#include "varswap3d/varswap3d.h"
#include "std/std.h"
#include "numfunc.h"
#include "transopt.h"
#include "pnl/pnl_random.h"

typedef struct VARSWAP3D_MOD{
  int Nb_factor;
  int is_call;
  int Strike;
  double T;
  double S0;
  double Divid;
  double R;
  double F0;
  double Bond;
  double V0;
  double V0_time;
  double V0_sqr;
  PnlVect * Beta;
  PnlVect * MeanReversion;
  PnlVect * SqrtMeanReversion;
  double Rho;
  double Sum_Beta; 
}VARSWAP3D_MOD;


VARSWAP3D_MOD * svs_model_create_from_Model(VARSWAP3D * Model);
void svs_model_initialise_from_Option(VARSWAP3D_MOD *  M,TYPEOPT *ptOpt);
void svs_model_initialise(VARSWAP3D_MOD * M);
void svs_sigma_time(VARSWAP3D_MOD * M, double T);
void svs_model_free(VARSWAP3D_MOD ** M);

#endif
