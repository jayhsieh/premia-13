#include "vol.h"

extern Option OPT(CallVixEuro);
extern Option OPT(PutVixEuro);
extern Option OPT(CallRealVarEuro);
extern Option OPT(PutRealVarEuro);
extern Option OPT(VarianceSwap);
extern Option OPT(CorrelationSwap);
extern Option OPT(CallVolatilityEuro);
extern Option OPT(PutVolatilityEuro);
extern Option OPT(VolatilitySwap);
extern Option OPT(Timer);

Option* OPT(family)[]=
{
  &OPT(CallRealVarEuro),
  &OPT(PutRealVarEuro),
  &OPT(CallVixEuro),
  &OPT(PutVixEuro),
  &OPT(VarianceSwap),
  &OPT(CorrelationSwap),
  &OPT(CallVolatilityEuro),
  &OPT(PutVolatilityEuro),
  &OPT(VolatilitySwap),
  &OPT(Timer),
  NULL
};
