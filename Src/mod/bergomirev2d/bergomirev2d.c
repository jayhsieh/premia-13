#include "bergomirev2d.h"
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
      pt->R.Val.V_DOUBLE=4.;
      pt->R.Viter=ALLOW;
      model->nvar++;
      
      pt->theta.Vname = "Theta";
      pt->theta.Vtype=PDOUBLE;
      pt->theta.Val.V_PDOUBLE=0.3;
      pt->theta.Viter=ALLOW;
      model->nvar++;

      pt->k1.Vname = "Mean Reversion Speed 1";
      pt->k1.Vtype=PDOUBLE;
      pt->k1.Val.V_PDOUBLE=0.3;
      pt->k1.Viter=ALLOW;
      model->nvar++;

      pt->k2.Vname = "Mean Reversion Speed 2";
      pt->k2.Vtype=PDOUBLE;
      pt->k2.Val.V_PDOUBLE=8.0;
      pt->k2.Viter=ALLOW;
      model->nvar++;

      /* pt->rhoxy.Vname = "Correlation between factors";
       * pt->rhoxy.Vtype=RGDOUBLEM11;
       * pt->rhoxy.Val.V_RGDOUBLEM11=0.;
       * pt->rhoxy.Viter=ALLOW;
       * model->nvar++; */

      pt->rhoSx.Vname = "Correlation between stock and factor 1";
      pt->rhoSx.Vtype=RGDOUBLEM11;
      pt->rhoSx.Val.V_RGDOUBLEM11=-0.6;
      pt->rhoSx.Viter=ALLOW;
      model->nvar++;
      
      pt->rhoSy.Vname = "Correlation between stock and factor 2";
      pt->rhoSy.Vtype=RGDOUBLEM11;
      pt->rhoSy.Val.V_RGDOUBLEM11=-0.7;
      pt->rhoSy.Viter=ALLOW;
      model->nvar++;

      pt->ForwardVarianceData.Vname = "ForwardVarianceData";
      pt->ForwardVarianceData.Vtype=FILENAME;
      pt->ForwardVarianceData.Val.V_INT=0;
      pt->ForwardVarianceData.Viter=FORBID;
      pt->ForwardVarianceData.Vsetable = SETABLE;

      model->nvar++;
      if (( pt->ForwardVarianceData.Val.V_FILENAME=malloc(sizeof(char)*MAX_PATH_LEN))==NULL)
        return MEMORY_ALLOCATION_FAILURE;
      sprintf( pt->ForwardVarianceData.Val.V_FILENAME, "%s%sForwardVarianceData.dat", premia_data_dir,path_sep);

    }

  return OK;
}

TYPEMOD bergomirev2d;
MAKEMOD(bergomirev2d);
