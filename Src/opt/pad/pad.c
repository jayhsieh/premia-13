#include "pad.h"
#include "error_msg.h"

extern Option OPT(AsianCallFloatingEuro);
extern Option OPT(AsianPutFloatingEuro);
extern Option OPT(AsianPutFixedEuro);
extern Option OPT(AsianCallFixedEuro);
extern Option OPT(LookBackCallFixedEuro);
extern Option OPT(LookBackCallFloatingEuro);
extern Option OPT(LookBackPutFixedEuro);
extern Option OPT(LookBackPutFloatingEuro);
extern Option OPT(AsianCallFloatingAmer);
extern Option OPT(AsianPutFloatingAmer);
extern Option OPT(AsianPutFixedAmer);
extern Option OPT(AsianCallFixedAmer);
extern Option OPT(LookBackCallFixedAmer);
extern Option OPT(LookBackCallFloatingAmer);
extern Option OPT(LookBackPutFixedAmer);
extern Option OPT(LookBackPutFloatingAmer);
extern Option OPT(Cliquet);
extern Option OPT(MovingAverageCallFixedAmer);
extern Option OPT(MovingAveragePutFixedAmer);
extern Option OPT(MovingAverageCallFloatingAmer);
extern Option OPT(MovingAveragePutFloatingAmer);

Option* OPT(family)[]={
  &OPT(LookBackCallFloatingEuro),
  &OPT(LookBackPutFloatingEuro),
  &OPT(LookBackPutFixedEuro),
  &OPT(LookBackCallFixedEuro),
  &OPT(AsianCallFixedEuro),
  &OPT(AsianPutFixedEuro),
  &OPT(AsianCallFloatingEuro),
  &OPT(AsianPutFloatingEuro),
  &OPT(LookBackCallFloatingAmer),
  &OPT(LookBackPutFloatingAmer),
  &OPT(LookBackPutFixedAmer),
  &OPT(LookBackCallFixedAmer),
  &OPT(AsianCallFixedAmer),
  &OPT(AsianPutFixedAmer),
  &OPT(AsianCallFloatingAmer),
  &OPT(AsianPutFloatingAmer),
  &OPT(Cliquet),
  &OPT(MovingAverageCallFixedAmer),
  &OPT(MovingAveragePutFixedAmer),
  &OPT(MovingAverageCallFloatingAmer),
  &OPT(MovingAveragePutFloatingAmer),
  NULL
};
