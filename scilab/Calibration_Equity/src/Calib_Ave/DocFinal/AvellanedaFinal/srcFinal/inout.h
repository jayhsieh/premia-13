#ifndef INOUT_H
#define INOUT_H

#include "tree.h"
/*
  MATHFI Project, Inria Rocquencourt.
*/


void loadParamCalib(char *name_in_calib, double *pt_S_0, double *pt_r, double *pt_q, int *pt_N ,\
 double *pt_sigma_0, double *pt_sigma_min, double *pt_sigma_max, double *pt_sigma_bar, \
 double *pt_gradtol,double *pt_steptol,int *pt_verbosity,int *pt_saveSuccessiveXinFile,\
 int *pt_maxCounter,double *pt_lambda,double *pt_alpha, char **pt_name_in_data,char **pt_name_out_Vol_Loc);

void loadParamPricer(char *name_in_pricer,double *pt_K, double *pt_T, char **pt_option_type,char ** pt_name_in_grid_price_in,char
**pt_name_in_grid_price_out);


void stringToString(char *strcopy, char *string);

int stringToInt(char *string);

double stringToDouble(char *string);

#endif
