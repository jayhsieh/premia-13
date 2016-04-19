#include "Hawkes_Intensity.h"
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

      pt->lambda0.Vname = "lambda_0";
      pt->lambda0.Vtype=DOUBLE;
      pt->lambda0.Val.V_DOUBLE=1.;
      pt->lambda0.Viter=ALLOW;
      model->nvar++;

      pt->kappa.Vname = "Long-Run Variance";
      pt->kappa.Vtype=DOUBLE;
      pt->kappa.Val.V_DOUBLE=2.;
      pt->kappa.Viter=ALLOW;
      model->nvar++;

      pt->c.Vname = "Mean reversion";
      pt->c.Vtype=DOUBLE;
      pt->c.Val.V_DOUBLE=1.;
      pt->c.Viter=ALLOW;
      model->nvar++;


      pt->delta.Vname = "Loss impact";
      pt->delta.Vtype=DOUBLE;
      pt->delta.Val.V_DOUBLE=1.;
      pt->delta.Viter=ALLOW;
      model->nvar++;
      
      pt->r.Vname = "Interest Rate";
      pt->r.Vtype=DOUBLE;
      pt->r.Val.V_DOUBLE=0.05;
      pt->r.Viter=ALLOW;
      model->nvar++;

      pt->Ncomp.Vname = "Number of Companies";
      pt->Ncomp.Vtype=PINT;
      pt->Ncomp.Val.V_PINT=125;
      pt->Ncomp.Viter=ALLOW;
      model->nvar++;
      
    }

  return OK;
}

TYPEMOD Hawkes_Intensity;
MAKEMOD(Hawkes_Intensity);
