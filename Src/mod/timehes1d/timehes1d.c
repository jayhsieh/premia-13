#include "timehes1d.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"
#include "premia_obj.h"

static int MOD(Init)(Model *model)       
{
  TYPEMOD* pt=(TYPEMOD*)(model->TypeModel);
  
  if (model->init == 0 )
    {
      model->init = 1;
      model->nvar=0;
      pt->T.Vname = "Current Date";
      pt->T.Vtype=DATE;
      pt->T.Val.V_DATE=0.;
      pt->T.Viter=ALLOW;
      model->nvar++;

      pt->S0.Vname = "Spot";
      pt->S0.Vtype=PDOUBLE;
      pt->S0.Val.V_PDOUBLE=100.;
      pt->S0.Viter=ALLOW;
      model->nvar++;

      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=DOUBLE;
      pt->Divid.Val.V_DOUBLE=0.;
      pt->Divid.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=5.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->Sigma0.Vname = "Current Variance";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=0.04;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;

      pt->MeanReversion.Vname = "Mean Reversion";
      pt->MeanReversion.Vtype=DOUBLE;
      pt->MeanReversion.Val.V_DOUBLE=3.;
      pt->MeanReversion.Viter=ALLOW;
      model->nvar++;


      pt->TimeDepParameters.Vname = "Piecewise Constant Parameters";
      pt->TimeDepParameters.Vtype=FILENAME;
      pt->TimeDepParameters.Val.V_FILENAME=NULL;
      pt->TimeDepParameters.Viter=FORBID;
      pt->TimeDepParameters.Vsetable = SETABLE;
      model->nvar++;

      if (( pt->TimeDepParameters.Val.V_FILENAME=malloc(sizeof(char)*MAX_PATH_LEN))==NULL)
        return MEMORY_ALLOCATION_FAILURE;
      sprintf( pt->TimeDepParameters.Val.V_FILENAME, "%s%sHeston_TimeDepParameters.dat", premia_data_dir,path_sep);
      
      pt->TimeStep.Vname = "Interval of constance";
      pt->TimeStep.Vtype = PDOUBLE;
      pt->TimeStep.Val.V_PDOUBLE = 0.25;
      pt->TimeStep.Viter = FORBID;
      pt->TimeStep.Vsetable = SETABLE;
      model->nvar++;
    }

  return OK;
}

TYPEMOD TimeHeston1dim;
MAKEMOD(TimeHeston1dim);

