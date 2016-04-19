#include "stdi.h"

extern Option OPT(ZeroCouponBond);
extern Option OPT(ZeroCouponCallBondEuro);
extern Option OPT(ZeroCouponPutBondEuro);
extern Option OPT(ZeroCouponCallBondAmer);
extern Option OPT(ZeroCouponPutBondAmer);
extern Option OPT(CouponBearingCallEuro);
extern Option OPT(Cap);
extern Option OPT(Floor);
extern Option OPT(Caplet);
extern Option OPT(PayerSwaption);
extern Option OPT(ReceiverSwaption);
extern Option OPT(PayerBermudanSwaption);
extern Option OPT(ReceiverBermudanSwaption);


Option* OPT(family)[]=
{
  &OPT(ZeroCouponBond),
  &OPT(ZeroCouponCallBondEuro),
  &OPT(ZeroCouponPutBondEuro),
  &OPT(ZeroCouponCallBondAmer),
  &OPT(ZeroCouponPutBondAmer),
  &OPT(CouponBearingCallEuro),
  &OPT(Cap),
  &OPT(Floor),
  &OPT(Caplet),
  &OPT(PayerSwaption),
  &OPT(ReceiverSwaption),
  &OPT(PayerBermudanSwaption),
  &OPT(ReceiverBermudanSwaption),
  NULL
};
