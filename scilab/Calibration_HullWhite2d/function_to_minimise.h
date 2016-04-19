#ifndef FUNCTION_TO_MINIMISE_H_INCLUDED
#define FUNCTION_TO_MINIMISE_H_INCLUDED

typedef struct {
  PnlVect *x_lower_bound;
  PnlVect *x_upper_bound;
} Boundaries;

void BoxConstraints_HW2d(const PnlVect *x, PnlVect *res, void *params);

double fitting_error_cap(const PnlVect *x, void *MktATMCapVol);

void  print_model_cap_vol(const PnlVect *x, MktCapZCData *MktATMCapVol);

double fitting_error_swaption(const PnlVect *x, void *MktATMSwaptionVol);

void  print_model_swaption_vol(const PnlVect *x, MktSwpZCData *MktATMSwaptionVol);

void ChooseInputParameters(PnlVect* x_input, PnlVect*lower_bounds, PnlVect* upper_bounds, int type_generator);

void ChooseRandomInputParameters(PnlVect* x_input, PnlVect*lower_bounds, PnlVect* upper_bounds, int type_generator);


#endif // FUNCTION_TO_MINIMISE_H_INCLUDED
