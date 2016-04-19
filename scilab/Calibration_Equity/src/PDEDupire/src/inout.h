#include "data.h"

/*
  MATHFI Project, Inria Rocquencourt.
  Sophie Volle and Jean-Marc Cognet, November 2002.
*/

void loadSigmaddl(char *name_in_sigma, int *pt_n, int *pt_m, double **pt_y_coarseGrid, double **pt_T_coarseGrid, double **pt_sigma);

void loadSigmaGrid(char *name_in_sigma, int *pt_n, int *pt_m, double **pt_y_coarseGrid, double **pt_T_coarseGrid);

void sigmaddlToSigmaFineGrid(double **sigmaFineGrid, int n, int m, double *y_coarseGrid, double *T_coarseGrid, double *sigma_param, double *y_fineGrid, double *T_fineGrid, int N, int M);

void loadDataPrices(char *name_in_data, struct marketData **data);

//void loadSigmaddl(char *name_in_sigma_ddl, int *pt_n, int *pt_m, double **y_coarseGrid, double **T_coarseGrid, double **sigma_ddl);
//void loadSigmaddl(char *name_in_sigma_ddl, int *pt_n, int *pt_m, double **y_coarseGrid);

void loadParamSimul(char *name_in_simul, double *pt_S_0, double *pt_r, double *pt_q, int *pt_optionType, int *pt_optionSimul, double *pt_t_0, double *pt_T_max, double *pt_y_min, double *pt_y_max, int *pt_N, int *pt_M, int *pt_gridType, double *pt_theta, char **pt_name_in_sigma_ddl, double *pt_sigmaCte, char **pt_name_in_visu, char **pt_name_out_visu, char **pt_name_in_data, char **pt_name_out_data);

void loadParamRafsig(char *name_in_rafsig, char **pt_name_in_sigma1_ddl, char **pt_nbsplit_y_char, char **pt_nbsplit_T_char, char **pt_name_out_sigma2_ddl);

void loadParamVisusig(char *name_in_visusig, char **pt_name_in_sigma_ddl, char **pt_name_in_sigma_visu, char **pt_name_out_sigma_visu);

void loadParamImpsig(char *name_in_impsig, double *pt_S_0, double *pt_r, double *pt_q, int *pt_optionType, double *pt_t_0, char **pt_name_in_data, char **pt_name_out_sigma);

void loadParamOptim(char *name_in_optim, double *pt_gradtol, double *pt_steptol, int *pt_verbosity, int *pt_saveSuccessiveXinFile, int *pt_maxCounter, double *pt_lambda);

void loadParamOptim2(char *name_in_optim, double *pt_pgtol, double *pt_factr, int *pt_iprint, int *pt_maxCounter, double *pt_sigma_min, double *pt_sigma_max, double *pt_lambda);

void loadParamCalib(char *name_in_calib, double *pt_S_0, double *pt_r, double *pt_q, int *pt_optionType, double *pt_t_0, double *pt_T_max, double *pt_y_min, double *pt_y_max, int *pt_N, int *pt_M, int *pt_gridType, double *pt_theta, int *pt_choice_optim, char **pt_name_in_optim, char **pt_name_in_data, char **pt_name_in_sigmainit_ddl, char **pt_name_out_sigmaest_ddl, char **pt_name_in_sigma_visu, char **pt_name_out_sigmainit_visu, char **pt_name_out_sigmaest_visu);

void stringToString(char *strcopy, char *string);

int stringToInt(char *string);

double stringToDouble(char *string);

void savePriceOrSigma(char *name_out_visu, double **priceOrSigma_visuGrid, int Nprice_visu, int Mprice_visu, double *Kprice_visu, double *Tprice_visu);

void saveSigmaddl(char *name_out_sigma_ddl, double *sigma_ddl, int n, int m, double *y_coarseGrid, double *T_coarseGrid);
