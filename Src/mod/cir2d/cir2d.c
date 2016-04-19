#include "cir2d.h"
#include "chk.h"
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
      pt->T.Val.V_DATE=0.0;
      pt->T.Viter=ALLOW;
      model->nvar++;

      pt->x01.Vname = "Current X1";
      pt->x01.Vtype=PDOUBLE;
      pt->x01.Val.V_PDOUBLE=0.04;
      pt->x01.Viter=ALLOW;
      model->nvar++;

      pt->x02.Vname = "Current X2";
      pt->x02.Vtype=PDOUBLE;
      pt->x02.Val.V_PDOUBLE=0.02;
      pt->x02.Viter=ALLOW;
      model->nvar++;

      pt->k1.Vname = "Speed of Mean Reversion 1";
      pt->k1.Vtype=PDOUBLE;
      pt->k1.Val.V_PDOUBLE=0.02;
      pt->k1.Viter=ALLOW;
      model->nvar++;

      pt->k2.Vname = "Speed of Mean Reversion 2";
      pt->k2.Vtype=PDOUBLE;
      pt->k2.Val.V_PDOUBLE=0.02;
      pt->k2.Viter=ALLOW;
      model->nvar++;

      pt->theta1.Vname = "Long Term Mean 1";
      pt->theta1.Vtype=PDOUBLE;
      pt->theta1.Val.V_PDOUBLE=0.03;
      pt->theta1.Viter=ALLOW;
      model->nvar++;

      pt->theta2.Vname = "Long Term Mean 2";
      pt->theta2.Vtype=PDOUBLE;
      pt->theta2.Val.V_PDOUBLE=0.01;
      pt->theta2.Viter=ALLOW;
      model->nvar++;

      pt->Sigma1.Vname = "Volatility 1";
      pt->Sigma1.Vtype=PDOUBLE;
      pt->Sigma1.Val.V_PDOUBLE=0.04;
      pt->Sigma1.Viter=ALLOW;
      model->nvar++;

      pt->Sigma2.Vname = "Volatility 2";
      pt->Sigma2.Vtype=PDOUBLE;
      pt->Sigma2.Val.V_PDOUBLE=0.02;
      pt->Sigma2.Viter=ALLOW;
      model->nvar++;

      pt->shift.Vname = "Initial Shift";
      pt->shift.Vtype=PDOUBLE;
      pt->shift.Val.V_PDOUBLE=0.02;
      pt->shift.Viter=ALLOW;
      model->nvar++;

    }
  return OK;
}
TYPEMOD Cir2d;
MAKEMOD(Cir2d);


