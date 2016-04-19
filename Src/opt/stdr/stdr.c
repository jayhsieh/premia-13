#include "stdr.h"

extern Option OPT(VaRisk);
extern Option OPT(CreditVaRisk);

Option* OPT(family)[]=
{
  &OPT(VaRisk),
  &OPT(CreditVaRisk),
  NULL
};
