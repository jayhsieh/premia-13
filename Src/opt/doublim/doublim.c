#include "doublim.h"
#include "error_msg.h"

extern Option OPT(DoubleCallOutEuro);
extern Option OPT(DoublePutInEuro);
extern Option OPT(DoublePutOutEuro);
extern Option OPT(DoubleCallInEuro);
extern Option OPT(DoubleCallOutAmer);
extern Option OPT(DoublePutInAmer);
extern Option OPT(DoublePutOutAmer);
extern Option OPT(DoubleCallInAmer);
extern Option OPT(ParisianDoubleCallOutEuro);
extern Option OPT(ParisianDoubleCallInEuro);

Option* OPT(family)[]={
  &OPT(DoubleCallOutEuro),
  &OPT(DoublePutOutEuro),
  &OPT(DoubleCallInEuro),
  &OPT(DoublePutInEuro),
  &OPT(DoubleCallOutAmer),
  &OPT(DoublePutOutAmer),
  &OPT(DoubleCallInAmer),
  &OPT(DoublePutInAmer),
  &OPT(ParisianDoubleCallOutEuro),
  &OPT(ParisianDoubleCallInEuro),
  NULL,
};
