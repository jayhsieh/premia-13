#include "bscir2d.h"
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

      pt->Sigma.Vname = "Volatility";
      pt->Sigma.Vtype=DOUBLE;
      pt->Sigma.Val.V_DOUBLE=0.1358;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->r0.Vname = "Current Rate";
      pt->r0.Vtype=PDOUBLE;
      pt->r0.Val.V_PDOUBLE=0.04;
      pt->r0.Viter=ALLOW;
      model->nvar++;

      pt->k.Vname = "Speed of Mean Reversion";
      pt->k.Vtype=PDOUBLE;
      pt->k.Val.V_PDOUBLE=1;
      pt->k.Viter=ALLOW;
      model->nvar++;

      pt->theta.Vname = "Long Term Mean";
      pt->theta.Vtype=PDOUBLE;
      pt->theta.Val.V_PDOUBLE=0.04;
      pt->theta.Viter=ALLOW;
      model->nvar++;
      
      pt->SigmaR.Vname = "Volatility of Interest Rate";
      pt->SigmaR.Vtype=PDOUBLE;
      pt->SigmaR.Val.V_PDOUBLE=0.2;
      pt->SigmaR.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Correlation";
      pt->Rho.Vtype=PDOUBLE;
      pt->Rho.Val.V_PDOUBLE=0;
      pt->Rho.Viter=ALLOW;
      model->nvar++;
      
      pt->Mortality.Vname = "Mortality Table";
      pt->Mortality.Vtype=FILENAME;
      pt->Mortality.Val.V_FILENAME=NULL;
      pt->Mortality.Viter=FORBID;
      pt->Mortality.Vsetable =SETABLE;

      model->nvar++;
      
      if ((pt->Mortality.Val.V_FILENAME=malloc(sizeof(char)*MAX_PATH_LEN))==NULL)
        return MEMORY_ALLOCATION_FAILURE;
      sprintf( pt->Mortality.Val.V_FILENAME, "%s%sMortality.dat", premia_data_dir,path_sep);
      
    }

  return OK;
}

TYPEMOD bscir2d;
MAKEMOD(bscir2d);
