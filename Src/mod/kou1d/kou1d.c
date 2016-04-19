#include "kou1d.h"
#include "chk.h"
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
      //pt->R.Val.V_DOUBLE=10.;
      pt->R.Val.V_DOUBLE=5.12711;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Sigma";
      pt->Sigma.Vtype=DOUBLE;
      //pt->Sigma.Val.V_DOUBLE=0.2;
      pt->Sigma.Val.V_DOUBLE=0.3;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->Lambda.Vname = "Intensity of Jump Lambda";
      pt->Lambda.Vtype=SPDOUBLE;
      //pt->Lambda.Val.V_SPDOUBLE=1.;
      //pt->Lambda.Val.V_SPDOUBLE=0.33;
      pt->Lambda.Val.V_SPDOUBLE=7;
      pt->Lambda.Viter=ALLOW;
      model->nvar++;

      pt->LambdaPlus.Vname = "LambdaPlus";
      pt->LambdaPlus.Vtype=RGDOUBLE1;
      //pt->LambdaPlus.Val.V_RGDOUBLE1=6.;
      //pt->LambdaPlus.Val.V_RGDOUBLE1=9.6;
      pt->LambdaPlus.Val.V_RGDOUBLE1=50.;
      pt->LambdaPlus.Viter=ALLOW;
      model->nvar++;

      pt->LambdaMinus.Vname = "LambdaMinus";
      pt->LambdaMinus.Vtype=SPDOUBLE;
      //pt->LambdaMinus.Val.V_SPDOUBLE=4.;
      //pt->LambdaMinus.Val.V_SPDOUBLE=3.1;
      pt->LambdaMinus.Val.V_SPDOUBLE=25;
      pt->LambdaMinus.Viter=ALLOW;
      model->nvar++;

      pt->P.Vname = "Probability of Positive Jump";
      pt->P.Vtype=RGDOUBLE;
      //pt->P.Val.V_RGDOUBLE=0.5;
      //pt->P.Val.V_RGDOUBLE=0.2;
      pt->P.Val.V_RGDOUBLE=0.6;
      pt->P.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}


TYPEMOD Kou1dim;

MAKEMOD(Kou1dim);


