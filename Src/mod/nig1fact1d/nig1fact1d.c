#include "nig1fact1d.h"
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

      pt->S0.Vname = "Initial Forward Price";
      pt->S0.Vtype=PDOUBLE;
      pt->S0.Val.V_PDOUBLE=100.;
      pt->S0.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=5.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->alpha.Vname = "alpha";
      pt->alpha.Vtype=SPDOUBLE;
      pt->alpha.Val.V_SPDOUBLE=15.81;
      pt->alpha.Viter=ALLOW;
      model->nvar++;

      pt->beta.Vname = "beta";
      pt->beta.Vtype=DOUBLE;
      pt->beta.Val.V_DOUBLE=-1.581;
      pt->beta.Viter=ALLOW;
      model->nvar++;

      pt->delta.Vname = "delta";
      pt->delta.Vtype=SPDOUBLE;
      pt->delta.Val.V_SPDOUBLE=15.57;
      pt->delta.Viter=ALLOW;
      model->nvar++;

      pt->mu.Vname = "mu";
      pt->mu.Vtype=SPDOUBLE;
      pt->mu.Val.V_SPDOUBLE=1.56;
      pt->mu.Viter=ALLOW;
      model->nvar++;

      
      pt->Sigma.Vname = "Sigma";
      pt->Sigma.Vtype=SPDOUBLE;
      pt->Sigma.Val.V_SPDOUBLE=0.5747;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->lambda.Vname = "lambda";
      pt->lambda.Vtype=SPDOUBLE;
      pt->lambda.Val.V_SPDOUBLE=3;
      pt->lambda.Viter=ALLOW;
      model->nvar++;
    }

  return OK;
}

TYPEMOD nig1fact1d;
MAKEMOD(nig1fact1d);
