#include "cgmy1d.h"
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

      pt->C.Vname = "C";
      pt->C.Vtype=PDOUBLE;
      pt->C.Val.V_PDOUBLE=0.5;
      pt->C.Val.V_PDOUBLE=0.024400;
      pt->C.Viter=ALLOW;
      model->nvar++;

      pt->G.Vname = "G";
      pt->G.Vtype=PDOUBLE;
      pt->G.Val.V_PDOUBLE=2;
      pt->G.Val.V_PDOUBLE=0.076500;
      pt->G.Viter=ALLOW;
      model->nvar++;

      pt->M.Vname = "M";
      pt->M.Vtype=PDOUBLE;
      pt->M.Val.V_PDOUBLE=3.5;
      pt->M.Val.V_PDOUBLE=7.551500;
      pt->M.Viter=ALLOW;
      model->nvar++;

      pt->Y.Vname = "Y";
      pt->Y.Vtype=PDOUBLE;
      pt->Y.Val.V_PDOUBLE=0.5;
      pt->Y.Val.V_PDOUBLE=1.294500;
      pt->Y.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}

TYPEMOD CGMY1dim;

MAKEMOD(CGMY1dim);


