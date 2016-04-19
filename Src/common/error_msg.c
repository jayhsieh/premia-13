#include <stdlib.h>
#include "error_msg.h"


char *error_msg[MAX_ERROR_MSG];

int InitErrorMsg()
{
  error_msg[MEMORY_ALLOCATION_FAILURE]="MEMORY_ALLOCATION_FAILURE";
  error_msg[time_bigger_than_the_last_time_value_entered_in_initialyield]="time_bigger_than_the_last_time_value_entered_in_initialyield";
  error_msg[UNABLE_TO_OPEN_FILE]="UNABLE_TO_OPEN_FILE";
  error_msg[OPTION_IRRELEVANT_TO_THIS_METHOD]="OPTION_IRRELEVANT_TO_THIS_METHOD";
  error_msg[NEGATIVE_PROBABILITY]="NEGATIVE_PROBABILITY";
  error_msg[UNABLE_TO_FIND_ACROBAT]="UNABLE_TO_FIND_ACROBAT";
  error_msg[DIMENSION_QMC_SEQ_EXCEDEED]="DIMENSION_QMC_SEQ_EXCEDEED";
  error_msg[STEP_NUMBER_TOO_SMALL]="STEP_NUMBER_TOO_SMALL";
  error_msg[PATH_TOO_LONG]="PATH_TOO_LONG";
  error_msg[BAD_ALPHA_TEMPSTABLE]="ALPHA+,ALPHA-,Y SHOULD BE DIFFERENT FROM 1";
  error_msg[NON_DEFINITE_MATRIX]="NON_DEFINITE_MATRIX";
  error_msg[BASIS_MAX_SIZE_EXCEEDED]="BASIS_MAX_SIZE_EXCEEDED";
  error_msg[BAD_TESSELATION_FORMAT]="BAD TESSELATION FORMAT";
  error_msg[UNTREATED_CASE]="UNTREATED_CASE";
  error_msg[AVAILABLE_IN_FULL_PREMIA]="This function is currently only available in the full \n version of Premia. For further information please go on Premia's Website";
  error_msg[PREMIA_UNTREATED_COPULA]="THIS TYPE OF COPULA IS NOT HANDLED";
  error_msg[PREMIA_UNTREATED_TAU_BHAR_CHIARELLA]="TAU VOLATILITY PARAMETER GREATER THAN BOND MATURITY";
  error_msg[ONLY_HOMOGENEOUS_CDO]="Only homogeneous CDOs can be priced";
  error_msg[HEGDING_MATURITY_GREATER_THAN_MATURITY]="HEGDING_MATURITY_GREATER_THAN_MATURITY";
  error_msg[MODEL_MAX_SIZE_EXCEEDED]="This method cannot be run with the current model size.\nTry with a smaller model size.";
  error_msg[ GAP_SHOULD_BE_NEGATIVE]=" GAP SHOULD BE NEGATIVE";
  return OK;
}

