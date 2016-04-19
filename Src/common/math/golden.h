#ifndef GOLDEN_H_INCLUDED
#define GOLDEN_H_INCLUDED

#include "pnl/pnl_mathtools.h"

double golden(PnlFunc * F, double ax, double bx, double cx, double tol, double *xmin);

#endif // GOLDEN_H_INCLUDED
