#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "pnl/pnl_mathtools.h"
#include "pnl/pnl_optim.h"

#include "levy_process.h"
#include "levy_diffusion.h"
#include "finance_tool_box.h"
#include "levy_calibration.h"
#include "levy_diffusion_calibration.h"
#include "carr.h"

// static void pnl_lbfgssolver_gradient_by_fd(double *f,
//                                            PnlVect *  grad,
//                                            const PnlVect * x,
//                                            PnlRnFuncGradR  *func)
// {
//   pnl_lbfgssolver_gradient_f(f,grad,x,func);
// };


int main()
{
  
  int iter_max=80;
  int print_algorithm_steps=1;
  int am;
  PnlVect *LevyParamsOutput;
  PnlVect *LevyParamsInit;
  Calibration_Data *data;
  List_Option_Eqd *op;
  PnlRnFuncR FuncToMinimize;
  PnlRnFuncRm Constraints;
  int type_diffusion;
  int type_model;
  double ri;
  double dividi;
  double tolerance = 1e-8;
  double S0;
  int size_model;
  

  am=0;
  

  printf(">> Initial Spot - (default value put 0) \n");
  scanf("%lf",&S0);
  if (S0==0)
    S0=2461.44;

  printf(">> Interest rate - (default value put -1) \n");
  scanf("%lf",&ri);
  if (ri==-1.)
    ri=0.03;

  printf(">> Dividend rate  - (default value put -1) \n");
  scanf("%lf",&dividi);
  if (dividi==-1.)
    dividi=0.0;
  
 
   

  //>>>>>  Initialisation :
  op=list_option_eqd_create(am,S0);
  list_option_eqd_set_rate(op,ri,dividi);
  list_option_eqd_readmarketdata(op,"Data/MarketData.dat");
  list_option_eqd_savemarketdata(op,"Data/MarketData_withvol.dat");

  //list_option_eqd_readmarketdata(op,"Data/MarketData1Yc.dat");
  // -- Minimisation Problem
  

  printf(">> Type of diffusion Levy 0, Stochastic Volatility 1 \n");
  scanf("%i",&type_diffusion);

  if(type_diffusion==0)
    {
      printf(">> Type of model : \n");
      printf(" CGMY            3 \n");
      printf(" Temperedstable  4 \n");
      printf(" VG              5 \n");
      printf(" Meixner         7 \n");
      scanf("%i",&type_model);
      switch (type_model)
        {
        case 3:
          {
            double default_parameters[4]={1.0,0.5,1.2,1.6};
            size_model=4;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 4:
          {
            double default_parameters[6]={0.5,0.5,6.,4.,1.4,1.4};
            size_model=6;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 5:
          {
            double default_parameters[3]={0.15,-0.5,0.2};
            size_model=3;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 7:
          {
            double default_parameters[3]={0.5,0.5,1.};
            size_model=3;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        default:
          size_model=0;
          break;
        }
      data=calibration_data_create(op,ri,dividi,type_model);
      FuncToMinimize.function  = &QuadraticError_ForLevyProcess;
      FuncToMinimize.params    = data;
      
      Constraints.function =&Constraints_ForLevyProcess;
      Constraints.params   = data;
      
    }
  
  if(type_diffusion==1)
    {
      printf(">> Type of model: \n");
      printf(" Heston  1 \n");
      printf(" Bates   2 \n");
      printf(" BNS     3 \n");
      printf(" DPS     4 \n");
      printf(" CIRNIG  5 \n");
      printf(" GOUNIG  6 \n");
      scanf("%i",&type_model);
      switch (type_model)
        {
        case 1:
          {
            double default_parameters[5]={0.1,0.5,0.1,0.1,0.2};
            size_model=5;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 2:
          {
            double default_parameters[8]={0.1,0.5,0.1,0.1,0.2,0.0,0.0,0.0};
            size_model=8;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 3:
          {
            double default_parameters[5]={0.1,-4.5,20.,1.0,0.2};
            size_model=5;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 4:
          {
            double default_parameters[15]={0.1,0.5,0.1,0.1,0.2,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01,0.01};
            size_model=15;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 5:
          {
            double default_parameters[7]={1.2,0.6,1.0,1.0,0.15,-0.5,0.2};
            size_model=7;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
        case 6:
          {
            double default_parameters[7]={2.0,0.5,0.8,1.0,0.15,-0.5,0.2};
            //Gamma CIR double default_parameters[7]={2.0,0.5,0.8,1.0,16.0,-4.0,1.2};
            size_model=7;
            LevyParamsInit = pnl_vect_create_from_ptr(size_model,default_parameters);
            break;
          }
      
        default:
          size_model=0;
        }
      
      data=calibration_data_create(op,ri,dividi,type_model);
      FuncToMinimize.function  = &QuadraticError_ForLevyDiffusion;
      FuncToMinimize.params    = data;
 
      Constraints.function =&Constraints_ForLevyDiffusion;
      Constraints.params   = data;
 
    }


  LevyParamsOutput = pnl_vect_create_from_double(size_model,0.);


  
  // Solve Quatratic minimisation problem
  pnl_optim_intpoints_bfgs_solve(&FuncToMinimize, NULL,
                                 &Constraints,
                                 NULL, NULL,
                                 LevyParamsInit,
                                 tolerance,
                                 iter_max,
                                 print_algorithm_steps,
                                 LevyParamsOutput);
  
  printf("Sortie : \n");
  printf("Erreur Quadratique = sqrt(f) = %f \n", sqrt(PNL_EVAL_RNFUNCR(&FuncToMinimize,LevyParamsOutput)));

  
  //GeneratePrices_ForLevyProcess(data,Levy);
  pnl_vect_print(LevyParamsOutput);

  //Compute implied volatility 
  list_option_eqd_compute_implied_vol(data->list_model);
  //store result in file ModelData
  list_option_eqd_savemarketdata(data->list_model,"Data/ModelData.dat");

  // Free Memory
  pnl_vect_free(&LevyParamsOutput);
  pnl_vect_free(&LevyParamsInit);  
  calibration_data_free(&data);
  list_option_eqd_free(&op);
  return 0;
};

