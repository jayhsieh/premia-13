#include "purejump1d.h"
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

      pt->Sigma.Vname = "Volatility";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.2;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      /*pt->Divid.Vname = "Annual Dividend Rate";
	pt->Divid.Vtype=DOUBLE;
	pt->Divid.Val.V_DOUBLE=0.;
	pt->Divid.Viter=ALLOW;
	model->nvar++;
      */

      pt->R.Vname = "Annual Interest Rate";
      pt->R.Vtype=DOUBLE;
      pt->R.Val.V_DOUBLE=10.;
      pt->R.Viter=ALLOW;
      model->nvar++;

      pt->Beta.Vname = "Beta";
      pt->Beta.Vtype=DOUBLE;
      pt->Beta.Val.V_DOUBLE=1.;
      pt->Beta.Viter=ALLOW;
      model->nvar++;

      pt->Nu.Vname = "Nu";
      pt->Nu.Vtype=DOUBLE;
      pt->Nu.Val.V_DOUBLE=500.;
      pt->Nu.Viter=ALLOW;
      model->nvar++;

    }

  return OK;
}


TYPEMOD PureJump1dim;

MAKEMOD(PureJump1dim);


