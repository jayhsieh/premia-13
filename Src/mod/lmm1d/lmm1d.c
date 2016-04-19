#include "lmm1d.h"
#include "chk.h"
#include "model.h"

extern char* path_sep;

static PremiaEnumMember nbfact_members[] =
{
     {"1: Flat Volatility",1},
     {"2: Second Volatility factor",2},
    { NULL, NULLINT}
};

static DEFINE_ENUM(nbfact,nbfact_members);

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

      pt->NbFactors.Vname = "Number of Factors";
      pt->NbFactors.Vtype=ENUM;
      pt->NbFactors.Val.V_ENUM.value=1;
      pt->NbFactors.Val.V_ENUM.members=&nbfact;
      pt->NbFactors.Viter=ALLOW;
      model->nvar++;

      pt->l0.Vname = "Flat Initial Libor Rates";
      pt->l0.Vtype=PDOUBLE;
      pt->l0.Val.V_PDOUBLE=0.05;
      pt->l0.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Flat Volatility Libor Rates";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.2;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

    }
  return OK;
}
TYPEMOD LMM1d;
MAKEMOD(LMM1d);


