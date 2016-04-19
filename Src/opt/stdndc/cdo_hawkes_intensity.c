#include "stdndc.h"
#include "error_msg.h"
#include "premia_obj.h"

static TYPEOPT CDO_HAWKES_INTENSITY =
{
	{"Number of Companies", PINT, {0}, FORBID, UNSETABLE},
	{"Maturity",DATE,{0},ALLOW,SETABLE},
	{"Homogeneous Nominals", ENUM, {0}, FORBID, SETABLE},
	{"Tranches", PNLVECT, {0}, FORBID, SETABLE},
	{"Type of Recovery",ENUM,{0},FORBID,SETABLE},   
	{"Recovery Parameters",PNLVECT,{0},FORBID,SETABLE},   
	{"Number of coupon payments per year",INT,{0},ALLOW,SETABLE},
	{"Current date", DATE, {0}, IRRELEVANT, UNSETABLE},
	{"Number of defaults at current date", INT, {0}, IRRELEVANT, UNSETABLE}

};

static int OPT(Init)(Option *opt,Model *mod)
{
	TYPEOPT* pt=(TYPEOPT*)(opt->TypeOpt);
	VAR* ptMod=(VAR*)(mod->TypeModel);

	/* get the size from the model */
	mod->Init(mod);
	pt->Ncomp.Val.V_PINT = ptMod[0].Val.V_PINT;

	if (opt->init == 0 )
	{
		opt->init = 1;
		opt->nvar = 9;

		pt->maturity.Val.V_DATE=5.0;
		pt->date.Val.V_DATE = 0.; /* useless but needs to be initialiased */
		pt->n_defaults.Val.V_INT = 0; /* useless but needs to be initialiased */
		pt->NbPayment.Val.V_INT=4;

		opt->nvar_setable = 2;

		pt->t_nominal.Viter = IRRELEVANT;
		pt->t_nominal.Vsetable = UNSETABLE;
		pt->t_nominal.Val.V_ENUM.value=1;
		pt->t_nominal.Val.V_ENUM.members=NULL;

		pt->tranch.Val.V_PNLVECT=NULL;

		pt->t_recovery.Viter = IRRELEVANT;
		pt->t_recovery.Vsetable = UNSETABLE;
		pt->t_recovery.Val.V_ENUM.members = NULL;

		pt->p_recovery.Viter=IRRELEVANT;
		pt->p_recovery.Vsetable = UNSETABLE;
		pt->p_recovery.Vtype = DOUBLE;
		pt->p_recovery.Val.V_DOUBLE=0.4;
	}

	/* tranches */
	if ((pt->tranch).Val.V_PNLVECT == NULL)
	{
		double tranches[5] = {0, 0.03, 0.06, 0.1, 1};
		if ((pt->tranch.Val.V_PNLVECT =
			pnl_vect_create_from_ptr (5, tranches))==NULL)
			return WRONG;
	}

	return OK;
}

MAKEOPT(CDO_HAWKES_INTENSITY);
