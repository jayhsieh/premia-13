#include  "exoi.h"
static NumFunc_1 call=
{
    Call,
    {{"Strike",PDOUBLE,{100},FORBID,UNSETABLE},{" ",PREMIA_NULLTYPE,{0},FORBID,SETABLE}},
    CHK_call
};

static TYPEOPT CallableRangeAccrual=
{
    {"Payoff",NUMFUNC_1,{0},FORBID,SETABLE},     /* PayOff; */
    {"First Exercise Date",DATE,{0},ALLOW,SETABLE}, /* FirstExerciseDate;*/
    {"Last Payment Date",DATE,{0},ALLOW,SETABLE},/* LastPaymentDate;*/
    {"Reset Period",PDOUBLE,{0},ALLOW,SETABLE},  /* Reset Period;*/
    {"Nominal Value",PDOUBLE,{0},ALLOW,SETABLE}, /* Nominal;*/
    {"Spread Rate",PDOUBLE,{0},ALLOW,UNSETABLE},     /* Spread Rate;*/
    {"Cap Rate",PDOUBLE,{0},ALLOW,UNSETABLE},      /* Cap Rate;*/
    {"Strike",PDOUBLE,{0},ALLOW,UNSETABLE},      /* Strike;*/
    {"Gearing",PDOUBLE,{0},ALLOW,UNSETABLE},     /* Gearing;*/
    {"Floor Rate",PDOUBLE,{0},ALLOW,UNSETABLE},      /* Floor;*/
    {"Fixed Rate",PDOUBLE,{0},ALLOW,SETABLE},      /* FixedRate;*/
    {"Lower Range Bound",PDOUBLE,{0},ALLOW,SETABLE},      /* LowerRangeBound;*/
    {"Upper Range Bound",PDOUBLE,{0},ALLOW,SETABLE},      /* UpperRangeBound;*/
    {"CMS1 Maturity",PDOUBLE,{0},ALLOW,UNSETABLE},      /* CMSMat1;*/
    {"CMS2 Maturity",PDOUBLE,{0},ALLOW,UNSETABLE},      /* CMSMat2;*/
};

static int OPT(Init)(Option *opt,Model *mod)
{
    TYPEOPT* pt=(  TYPEOPT*)(opt->TypeOpt);

    if ( opt->init == 0)
    {
        opt->init = 1;
        opt->nvar = 15;
        opt->nvar_setable=8;

        pt->PayOff.Val.V_NUMFUNC_1=&call;

        (pt->FirstExerciseDate).Val.V_DATE=1.0;
        (pt->LastPaymentDate).Val.V_DATE=10.0;
        (pt->ResetPeriod).Val.V_PDOUBLE=1.0;
        (pt->Nominal).Val.V_PDOUBLE=1.0;
        (pt->Spread).Val.V_PDOUBLE=0.05;
        (pt->Cap).Val.V_PDOUBLE=10;
        (pt->Strike).Val.V_PDOUBLE=0.12;
        (pt->Gearing).Val.V_PDOUBLE=2.;
        (pt->Floor).Val.V_PDOUBLE=0.;
        (pt->FixedRate).Val.V_PDOUBLE=0.07;
        (pt->LowerRangeBound).Val.V_PDOUBLE=0.;
        (pt->UpperRangeBound).Val.V_PDOUBLE=0.045;
        (pt->CMSMat1).Val.V_PDOUBLE=10;
        (pt->CMSMat2).Val.V_PDOUBLE=2;
    }

    return OK;
}

MAKEOPT(CallableRangeAccrual);
