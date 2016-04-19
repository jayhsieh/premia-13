#include "lim.h"
#include "error_msg.h"

extern Option OPT(PutUpOutEuro);
extern Option OPT(CallDownOutAmer);
extern Option OPT(CallDownOutEuro);
extern Option OPT(CallUpInEuro);
extern Option OPT(CallUpInAmer);
extern Option OPT(CallUpOutAmer);
extern Option OPT(CallUpOutEuro);
extern Option OPT(PutDownInEuro);
extern Option OPT(PutDownInAmer);
extern Option OPT(PutDownOutAmer);
extern Option OPT(PutDownOutEuro);
extern Option OPT(PutUpInEuro);
extern Option OPT(PutUpInAmer);
extern Option OPT(PutUpOutAmer);
extern Option OPT(CallDownInEuro);
extern Option OPT(CallDownInAmer);
extern Option OPT(ParisianCallDownOutEuro);
extern Option OPT(ParisianCallDownInEuro);
extern Option OPT(ParisianPutDownOutEuro);
extern Option OPT(ParisianPutDownInEuro);
extern Option OPT(ParisianCallUpOutEuro);
extern Option OPT(ParisianCallUpInEuro);
extern Option OPT(ParisianPutUpOutEuro);
extern Option OPT(ParisianPutUpInEuro);

Option* OPT(family)[]={
  &OPT(CallDownOutEuro),
  &OPT(PutDownOutEuro),
  &OPT(CallUpOutEuro),
  &OPT(PutUpOutEuro),
  &OPT(CallDownInEuro),
  &OPT(CallUpInEuro),
  &OPT(PutDownInEuro),
  &OPT(PutUpInEuro),
  &OPT(CallDownOutAmer),
  &OPT(PutDownOutAmer),
  &OPT(CallUpOutAmer),
  &OPT(PutUpOutAmer),
  &OPT(CallDownInAmer),
  &OPT(PutDownInAmer),
  &OPT(CallUpInAmer),
  &OPT(PutUpInAmer),
  &OPT(ParisianCallDownOutEuro),
  &OPT(ParisianCallDownInEuro),
  &OPT(ParisianPutDownOutEuro),
  &OPT(ParisianPutDownInEuro),
  &OPT(ParisianCallUpOutEuro),
  &OPT(ParisianCallUpInEuro),
  &OPT(ParisianPutUpOutEuro),
  &OPT(ParisianPutUpInEuro),
  NULL
};
