#include "garch1d.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"
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


      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=10.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->alpha0.Vname = "alpha0";
      pt->alpha0.Vtype=DOUBLE;
      pt->alpha0.Val.V_DOUBLE=0.00001;
      pt->alpha0.Viter=ALLOW;
      model->nvar++;

      pt->alpha1.Vname = "alpha1";
      pt->alpha1.Vtype=DOUBLE;
      pt->alpha1.Val.V_DOUBLE=0.2;
      pt->alpha1.Viter=ALLOW;
      model->nvar++;

      pt->lambda.Vname = "lambda";
      pt->lambda.Vtype=DOUBLE;
      pt->lambda.Val.V_DOUBLE=0.01;
      pt->lambda.Viter=ALLOW;
      model->nvar++;

      pt->beta1.Vname = "beta1";
      pt->beta1.Vtype=DOUBLE;
      pt->beta1.Val.V_DOUBLE=0.7;
      pt->beta1.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}

TYPEMOD garch1d;
MAKEMOD(garch1d);
