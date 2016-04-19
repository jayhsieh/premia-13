#include "lmm_heston1d.h"
#include "chk.h"
#include "model.h"

extern char* path_sep;

static PremiaEnumMember nbfacthes_members[] =
{
     {"1:Flat Volatility",1},
     {"2:Second Volatility factor: 1./sqrt(0.04+0.00075*t) * (0.01 - 0.05*exp(-0.1*(T-t)))",2},
    { NULL, NULLINT}
};

static DEFINE_ENUM(nbfacthes,nbfacthes_members);

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
      pt->NbFactors.Val.V_ENUM.members=&nbfacthes;
      pt->NbFactors.Viter=ALLOW;
      model->nvar++;

      pt->l0.Vname = "Flat Initial Libor Rates";
      pt->l0.Vtype=PDOUBLE;
      pt->l0.Val.V_PDOUBLE=0.05;
      pt->l0.Viter=ALLOW;
      model->nvar++;

      pt->Sigma.Vname = "Flat Volatility Libor Rates ";
      pt->Sigma.Vtype=PDOUBLE;
      pt->Sigma.Val.V_PDOUBLE=0.2;
      pt->Sigma.Viter=ALLOW;
      model->nvar++;

      pt->Sigma0.Vname = "Current Variance";
      pt->Sigma0.Vtype=DOUBLE;
      pt->Sigma0.Val.V_DOUBLE=1.0;
      pt->Sigma0.Viter=ALLOW;
      model->nvar++;

      pt->MeanReversion.Vname = "Mean Reversion";
      pt->MeanReversion.Vtype=DOUBLE;
      pt->MeanReversion.Val.V_DOUBLE=1.;
      pt->MeanReversion.Viter=ALLOW;
      model->nvar++;

      pt->LongRunVariance.Vname = "Long-Run Variance";
      pt->LongRunVariance.Vtype=DOUBLE;
      pt->LongRunVariance.Val.V_DOUBLE=1.;
      pt->LongRunVariance.Viter=ALLOW;
      model->nvar++;

      pt->Sigma2.Vname = "Volatility of Volatility";
      pt->Sigma2.Vtype=DOUBLE;
      pt->Sigma2.Val.V_DOUBLE=0.6;
      pt->Sigma2.Viter=ALLOW;
      model->nvar++;

      pt->Rho1.Vname = "Rho 1";
      pt->Rho1.Vtype=DOUBLE;
      pt->Rho1.Val.V_DOUBLE=0.5;
      pt->Rho1.Viter=ALLOW;
      model->nvar++;

      pt->Rho2.Vname = "Rho 2: Only in the Second Factor Case";
      pt->Rho2.Vtype=DOUBLE;
      pt->Rho2.Val.V_DOUBLE=0.2;
      pt->Rho2.Viter=ALLOW;
      model->nvar++;


    }
  return OK;
}
TYPEMOD LMM_HESTON1d;
MAKEMOD(LMM_HESTON1d);


