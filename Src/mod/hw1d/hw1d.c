#include "hw1d.h"
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

      pt->Sigma0.Vname = "Current Variance";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=0.01;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;

      pt->Mean.Vname = "Mean";
      pt->Mean.Vtype=DOUBLE;
      pt->Mean.Val.V_DOUBLE=0.;
      pt->Mean.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Volatility of Volatility";
      pt->Sigma.Vtype=DOUBLE;
      pt->Sigma.Val.V_DOUBLE=0.2;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Rho";
      pt->Rho.Vtype=DOUBLE;
      pt->Rho.Val.V_DOUBLE=0.5;
      pt->Rho.Viter=ALLOW;
      model->nvar++;

      model->HelpFilenameHint = "HW1D";

    }

  return OK;
}


TYPEMOD HullWhite1dim;

MAKEMOD(HullWhite1dim);


