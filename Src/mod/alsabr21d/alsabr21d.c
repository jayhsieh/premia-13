#include "alsabr21d.h"
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
      model->nvar++;;

      pt->z0.Vname = "z0";
      pt->z0.Vtype=SPDOUBLE;
      pt->z0.Val.V_SPDOUBLE=0.05;
      pt->z0.Viter=ALLOW;
      model->nvar++;;

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

      pt->mu.Vname = "mu";
      pt->mu.Vtype=SNDOUBLE;
      pt->mu.Val.V_SNDOUBLE=-0.00005;
      pt->mu.Viter=ALLOW;
      model->nvar++;

      pt->eta.Vname = "eta";
      pt->eta.Vtype=SPDOUBLE;
      pt->eta.Val.V_SPDOUBLE=0.1;
      pt->eta.Viter=ALLOW;
      model->nvar++;

      pt->a1.Vname = "a1";
      pt->a1.Vtype=SDOUBLE2;
      pt->a1.Val.V_SDOUBLE2=3.0;
      pt->a1.Viter=ALLOW;
      model->nvar++;

      pt->a2.Vname = "a2";
      pt->a2.Vtype=SPDOUBLE;
      pt->a2.Val.V_SPDOUBLE=0.5;
      pt->a2.Viter=ALLOW;
      model->nvar++;

      pt->c1.Vname = "c1";
      pt->c1.Vtype=PDOUBLE;
      pt->c1.Val.V_PDOUBLE=1.5;
      pt->c1.Viter=ALLOW;
      model->nvar++;

      pt->c2.Vname = "c2";
      pt->c2.Vtype=PDOUBLE;
      pt->c2.Val.V_PDOUBLE=1.5;
      pt->c2.Viter=ALLOW;
      model->nvar++;
    }

  return OK;
}

TYPEMOD alsabr21d;
MAKEMOD(alsabr21d);
