#include "cirpp2d.h"
#include "chk.h"
#include "model.h"
#include "enums.h"

extern char* path_sep;

static PremiaEnumMember flatint_members[] =
{
    { "to Initial Yields and Intensity in data/ directory", 0 },
    { "to Initial Yields and Spread in data/ directory", 1 },
    { NULL, NULLINT }
};

DEFINE_ENUM(flatint, flatint_members);

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

      pt->flat_flag.Vname = "Calibration";
      pt->flat_flag.Vtype=ENUM;
      pt->flat_flag.Val.V_ENUM.value=0;
      pt->flat_flag.Val.V_ENUM.members=&flatint;
      pt->flat_flag.Viter=ALLOW;
      model->nvar++;

      pt->InitialYieldsR.Vname = "Initial R0";
      pt->InitialYieldsR.Vtype=PDOUBLE;
      pt->InitialYieldsR.Val.V_PDOUBLE=0.05;
      pt->InitialYieldsR.Viter=ALLOW;
      model->nvar++;
      
      pt->aR.Vname = "Speed of Mean Reversion Interest Rate";
      pt->aR.Vtype=PDOUBLE;
      pt->aR.Val.V_PDOUBLE=0.15;
      pt->aR.Viter=ALLOW;
      model->nvar++;
      
      pt->bR.Vname = "Long Term Mean Interest Rate";
      pt->bR.Vtype=PDOUBLE;
      pt->bR.Val.V_PDOUBLE=0.05;
      pt->bR.Viter=ALLOW;
      model->nvar++;

      pt->SigmaR.Vname = "Volatility Interest Rate";
      pt->SigmaR.Vtype=PDOUBLE;
      pt->SigmaR.Val.V_PDOUBLE=0.1;
      pt->SigmaR.Viter=ALLOW;
      model->nvar++;

      pt->InitialYieldsI.Vname = "Initial I0";
      pt->InitialYieldsI.Vtype=PDOUBLE;
      pt->InitialYieldsI.Val.V_PDOUBLE=0.05;
      pt->InitialYieldsI.Viter=ALLOW;
      model->nvar++;

      pt->aI.Vname = "Speed of Mean Reversion Intensity";
      pt->aI.Vtype=PDOUBLE;
      pt->aI.Val.V_PDOUBLE=0.15;
      pt->aI.Viter=ALLOW;
      model->nvar++;
      
      pt->bI.Vname = "Long Term Mean Interest Intensity";
      pt->bI.Vtype=PDOUBLE;
      pt->bI.Val.V_PDOUBLE=0.05;
      pt->bI.Viter=ALLOW;
      model->nvar++;

      pt->SigmaI.Vname = "Volatility Intensity";
      pt->SigmaI.Vtype=PDOUBLE;
      pt->SigmaI.Val.V_PDOUBLE=0.1;
      pt->SigmaI.Viter=ALLOW;
      model->nvar++;

      pt->Rho.Vname = "Rho";
      pt->Rho.Vtype=DOUBLE;
      pt->Rho.Val.V_DOUBLE=0.5;
      pt->Rho.Viter=ALLOW;
      model->nvar++;
      
    }
  return OK;
}
TYPEMOD CirPlus2d;
MAKEMOD(CirPlus2d);


