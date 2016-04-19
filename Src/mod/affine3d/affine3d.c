#include "affine3d.h"
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
      pt->x01.Vtype=DOUBLE;
      pt->x01.Val.V_DOUBLE=0.01;
      pt->x01.Viter=ALLOW;
      model->nvar++;

      pt->x02.Vname = "Current X2";
      pt->x02.Vtype=DOUBLE;
      pt->x02.Val.V_DOUBLE=0.005;
      pt->x02.Viter=ALLOW;
      model->nvar++;

      pt->x03.Vname = "Current X3";
      pt->x03.Vtype=DOUBLE;
      pt->x03.Val.V_DOUBLE=-0.02;
      pt->x03.Viter=ALLOW;
      model->nvar++;

      pt->k1.Vname = "Speed of Mean Reversion 1";
      pt->k1.Vtype=PDOUBLE;
      pt->k1.Val.V_PDOUBLE=1;
      pt->k1.Viter=ALLOW;
      model->nvar++;

      pt->k2.Vname = "Speed of Mean Reversion 2";
      pt->k2.Vtype=PDOUBLE;
      pt->k2.Val.V_PDOUBLE=0.2;
      pt->k2.Viter=ALLOW;
      model->nvar++;
      
      pt->k3.Vname = "Speed of Mean Reversion 3";
      pt->k3.Vtype=PDOUBLE;
      pt->k3.Val.V_PDOUBLE=0.5;
      pt->k3.Viter=ALLOW;
      model->nvar++;
     
      pt->Sigma1.Vname = "Volatility 1";
      pt->Sigma1.Vtype=PDOUBLE;
      pt->Sigma1.Val.V_PDOUBLE=0.01;
      pt->Sigma1.Viter=ALLOW;
      model->nvar++;

      pt->Sigma2.Vname = "Volatility 2";
      pt->Sigma2.Vtype=PDOUBLE;
      pt->Sigma2.Val.V_PDOUBLE=0.005;
      pt->Sigma2.Viter=ALLOW;
      model->nvar++;

      pt->Sigma3.Vname = "Volatility 3";
      pt->Sigma3.Vtype=PDOUBLE;
      pt->Sigma3.Val.V_PDOUBLE=0.002;
      pt->Sigma3.Viter=ALLOW;
      model->nvar++;

      pt->shift.Vname = "Initial Shift";
      pt->shift.Vtype=PDOUBLE;
      pt->shift.Val.V_PDOUBLE=0.06;
      pt->shift.Viter=ALLOW;
      model->nvar++;

      pt->Rho12.Vname = "Rho 12";
      pt->Rho12.Vtype=RGDOUBLEM11;
      pt->Rho12.Val.V_RGDOUBLEM11=-0.2;
      pt->Rho12.Viter=ALLOW;
      model->nvar++;

      pt->Rho13.Vname = "Rho 13";
      pt->Rho13.Vtype=RGDOUBLEM11;
      pt->Rho13.Val.V_RGDOUBLEM11=-0.1;
      pt->Rho13.Viter=ALLOW;
      model->nvar++;

      pt->Rho23.Vname = "Rho 23";
      pt->Rho23.Vtype=RGDOUBLEM11;
      pt->Rho23.Val.V_RGDOUBLEM11=0.3;
      pt->Rho23.Viter=ALLOW;
      model->nvar++;
    }
  return OK;
}
TYPEMOD Affine3d;
MAKEMOD(Affine3d);


