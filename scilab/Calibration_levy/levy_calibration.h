#ifndef _LEVY_CALIBRATION_H_
#define _LEVY_CALIBRATION_H_
#include <stdio.h>
#include <stdlib.h>

#include "pnl/pnl_matrix.h"
#include "pnl/pnl_vector.h"
#include "lbfgssolver.h"

#include "finance_tool_box.h"


extern double TT_interest_rate(double t);//0.03
extern double TT_volatility(double t); //0.15*0.15;}

typedef struct _Calibration_Data Calibration_Data;
struct _Calibration_Data
{
  List_Option_Eqd * list_input;
  List_Option_Eqd * list_model;
  int type_of_process;
};
  
extern Calibration_Data * calibration_data_create(List_Option_Eqd * list_input_,double r,double divid,int type_of_model_);
extern void calibration_data_free(Calibration_Data ** data);
extern double calibration_data_QuadraticError(const Calibration_Data * Data);

extern double QuadraticError_ForLevyProcess(const PnlVect *GenerationParams,void *data);
extern void Constraints_ForLevyProcess(const PnlVect *x, PnlVect *res, void *data);

extern double QuadraticError_ForLevyDiffusion(const PnlVect *GenerationParams,void *data);
extern void Constraints_ForLevyDiffusion(const PnlVect *x, PnlVect *res, void *data);


#endif // LEVY_CALIBRATION_H

