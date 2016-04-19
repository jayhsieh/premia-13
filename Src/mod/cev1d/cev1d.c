#include "cev1d.h"
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

      pt->Divid.Vname = "Annual Dividend Rate";
      pt->Divid.Vtype=DOUBLE;
      pt->Divid.Val.V_DOUBLE=0.;
      pt->Divid.Viter=ALLOW;
      model->nvar++;

      pt->v.Vname = "v";
      pt->v.Vtype=PDOUBLE;
      pt->v.Val.V_PDOUBLE=0.2;
      pt->v.Viter=ALLOW;
      model->nvar++;

      pt->beta.Vname = "beta";
      pt->beta.Vtype=PDOUBLE;
      pt->beta.Val.V_PDOUBLE=0.8;
      pt->beta.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}

TYPEMOD cev1d;
MAKEMOD(cev1d);
