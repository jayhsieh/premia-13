#include "exoi.h"
#include "error_msg.h"

extern Option OPT(CallableCappedFloater);
extern Option OPT(CallableInverseFloater);
extern Option OPT(CallableRangeAccrual);
extern Option OPT(CallableCMSSpread);

Option* OPT(family)[]=
{
    &OPT(CallableCappedFloater),
    &OPT(CallableInverseFloater),
    &OPT(CallableRangeAccrual),
    &OPT(CallableCMSSpread),
    NULL
};
