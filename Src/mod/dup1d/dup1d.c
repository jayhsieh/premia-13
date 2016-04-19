#include "dup1d.h"
#include "chk.h"
#include "error_msg.h"
#include "model.h"
#include "enums.h"

extern char* path_sep;

static PremiaEnumMember volatility_members[] = 
  {
    { "15/x", 0 },
    { "0.01+0.1*exp(-x/100)+0.01*t", 1 },
    { NULL, NULLINT }
};

static DEFINE_ENUM(volatility_type, volatility_members);

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

      pt->Sigma.Vname = "Volatility type";
      pt->Sigma.Vtype=ENUM;
      pt->Sigma.Val.V_ENUM.value=0;
      pt->Sigma.Val.V_ENUM.members=&volatility_type;
      pt->Sigma.Viter=ALLOW;
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

      model->HelpFilenameHint = "DUP1D";

    }

  return OK;
}


TYPEMOD Dupire1dim;

MAKEMOD(Dupire1dim);


