#include "variancegamma2d.h"
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

      pt->S01.Vname = "Spot 1";
      pt->S01.Vtype=PDOUBLE;
      pt->S01.Val.V_PDOUBLE=100.;
      pt->S01.Viter=ALLOW;
      model->nvar++;

       pt->S02.Vname = "Spot 2";
      pt->S02.Vtype=PDOUBLE;
      pt->S02.Val.V_PDOUBLE=100.;
      pt->S02.Viter=ALLOW;
      model->nvar++;

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=10.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->ap.Vname = "a+";
      pt->ap.Vtype=PDOUBLE;
      pt->ap.Val.V_PDOUBLE=20.4499;
      pt->ap.Viter=ALLOW;
      model->nvar++;

      pt->am.Vname = "a-";
      pt->am.Vtype=PDOUBLE;
      pt->am.Val.V_PDOUBLE=24.4499;
      pt->am.Viter=ALLOW;
      model->nvar++;

      pt->lambda.Vname = "lambda";
      pt->lambda.Vtype=PDOUBLE;
      pt->lambda.Val.V_PDOUBLE=10.;
      pt->lambda.Viter=ALLOW;
      model->nvar++;

      pt->alpha.Vname = "alpha";
      pt->alpha.Vtype=RGDOUBLE;
      pt->alpha.Val.V_RGDOUBLE=0.4;
      pt->alpha.Viter=ALLOW;
      model->nvar++;
    
    }

  return OK;
}

TYPEMOD variancegamma2d;
MAKEMOD(variancegamma2d);
