#include "std2d.h"
#include "error_msg.h"

extern Option OPT(BestOfEuro);
extern Option OPT(CallMaximumAmer);
extern Option OPT(CallMaximumEuro);
extern Option OPT(ExchangeAmer);
extern Option OPT(ExchangeEuro);
extern Option OPT(PutMinimumAmer);
extern Option OPT(PutMinimumEuro);
extern Option OPT(BestOfAmer);

Option* OPT(family)[]={
  &OPT(BestOfEuro),
  &OPT(CallMaximumEuro),
  &OPT(PutMinimumEuro),
  &OPT(ExchangeEuro),
  &OPT(BestOfAmer),
  &OPT(CallMaximumAmer),
  &OPT(PutMinimumAmer),
  &OPT(ExchangeAmer),
  NULL
};
