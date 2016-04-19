#include "ou1d.h"
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

      pt->Speed.Vname = "Speed";
      pt->Speed.Vtype=SPDOUBLE;
      pt->Speed.Val.V_SPDOUBLE=0.301;
      pt->Speed.Viter=ALLOW;
      model->nvar++;
     
      pt->Sigma.Vname = "Sigma";
      pt->Sigma.Vtype=SPDOUBLE;
      pt->Sigma.Val.V_SPDOUBLE=0.334;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->a1.Vname = "a1";
      pt->a1.Vtype=SPDOUBLE;
      pt->a1.Val.V_SPDOUBLE=100;
      pt->a1.Viter=ALLOW;
      model->nvar++;

      pt->a2.Vname = "a2";
      pt->a2.Vtype=SPDOUBLE;
      pt->a2.Val.V_SPDOUBLE=1;
      pt->a2.Viter=ALLOW;
      model->nvar++;

      pt->a3.Vname = "a3";
      pt->a3.Vtype=SPDOUBLE;
      pt->a3.Val.V_SPDOUBLE=1;
      pt->a3.Viter=ALLOW;
      model->nvar++;
      
      model->HelpFilenameHint = "OU1D";

    }

  return OK;
}

TYPEMOD ou1d;
MAKEMOD(ou1d);
