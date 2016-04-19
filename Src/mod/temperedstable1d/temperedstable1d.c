#include "temperedstable1d.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"

extern char* path_sep;



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

      pt->Mu.Vname = "Trend";
      pt->Mu.Vtype=DOUBLE;
      pt->Mu.Val.V_DOUBLE=0.;
      pt->Mu.Viter=ALLOW;
      model->nvar++;

      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=DOUBLE;
      pt->Divid.Val.V_DOUBLE=0.;
      pt->Divid.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=10.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->AlphaPlus.Vname = "AlphaPlus";
      pt->AlphaPlus.Vtype=RGDOUBLE02;
      pt->AlphaPlus.Val.V_RGDOUBLE02=0.5;
      pt->AlphaPlus.Viter=ALLOW;
      model->nvar++;

      pt->AlphaMinus.Vname = "AlphaMinus";
      pt->AlphaMinus.Vtype=RGDOUBLE02;
      pt->AlphaMinus.Val.V_RGDOUBLE02=0.5;
      pt->AlphaMinus.Viter=ALLOW;
      model->nvar++;

      pt->LambdaPlus.Vname = "LambdaPlus";
      pt->LambdaPlus.Vtype=RGDOUBLE1;
      pt->LambdaPlus.Val.V_RGDOUBLE1=6.;
      pt->LambdaPlus.Viter=ALLOW;
      model->nvar++;

      pt->LambdaMinus.Vname = "LambdaMinus";
      pt->LambdaMinus.Vtype=SPDOUBLE;
      pt->LambdaMinus.Val.V_SPDOUBLE=4.;
      pt->LambdaMinus.Viter=ALLOW;
      model->nvar++;

      pt->CPlus.Vname = "CPlus";
      pt->CPlus.Vtype=PDOUBLE;
      pt->CPlus.Val.V_PDOUBLE=1.;
      pt->CPlus.Viter=ALLOW;
      model->nvar++;

      pt->CMinus.Vname = "Cminus";
      pt->CMinus.Vtype=PDOUBLE;
      pt->CMinus.Val.V_PDOUBLE=1.;
      pt->CMinus.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}


TYPEMOD TemperedStable1dim;

MAKEMOD(TemperedStable1dim);


