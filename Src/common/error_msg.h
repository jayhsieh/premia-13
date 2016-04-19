#ifndef  _ERROR_MSG_H
#define _ERROR_MSG_H

#include "pnl/pnl_mathtools.h"

/*DONT USE 1 RESERVED TO OK*/
#define MEMORY_ALLOCATION_FAILURE 2
#define UNABLE_TO_OPEN_FILE 3
#define time_bigger_than_the_last_time_value_entered_in_initialyield 4
#define STEP_NUMBER_TOO_SMALL 5
#define PATH_TOO_LONG 6
#define DIMENSION_QMC_SEQ_EXCEDEED 7
#define BAD_ALPHA_TEMPSTABLE 8
#define NON_DEFINITE_MATRIX 9
#define NEGATIVE_PROBABILITY 10
#define BASIS_MAX_SIZE_EXCEEDED 11
#define BAD_TESSELATION_FORMAT 12
#define UNTREATED_CASE 13
#define AVAILABLE_IN_FULL_PREMIA 14
#define PREMIA_UNTREATED_COPULA 15
#define PREMIA_UNTREATED_TAU_BHAR_CHIARELLA 16
#define ONLY_HOMOGENEOUS_CDO 17
#define HEGDING_MATURITY_GREATER_THAN_MATURITY 18
#define OPTION_IRRELEVANT_TO_THIS_METHOD 20
#define MODEL_MAX_SIZE_EXCEEDED 21
#define GAP_SHOULD_BE_NEGATIVE 22

#define UNABLE_TO_FIND_ACROBAT 30

/* This last flag should be greater than the maximum of the previous ones */
#define MAX_ERROR_MSG 31

extern char *error_msg[MAX_ERROR_MSG];

int InitErrorMsg();

#endif
