#ifndef  _ERROR_MSG_H
#define _ERROR_MSG_H


#define MEMORY_ALLOCATION_FAILURE 2
#define UNABLE_TO_OPEN_FILE 3

#define NEGATIVE_PROBABILITY 10

#define OPTION_IRRELEVANT_TO_THIS_METHOD 20

#define UNABLE_TO_FIND_ACROBAT 30

#define DIMENSION_QMC_SEQ_EXCEDEED 7

#define STEP_NUMBER_TOO_SMALL 5

#define PATH_TOO_LONG 6


/* This last flag should be greater than the maximum of the previous ones */
#define MAX_ERROR_MSG 31


int InitErrorMsg();
void FreeErrorMsg();


#endif
