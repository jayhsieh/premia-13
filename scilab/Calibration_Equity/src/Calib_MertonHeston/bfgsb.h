#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

#include "costFunction.h"
#include "gradFunction.h"
#include "paramsbfgs.h"
#include "./utils/minix.h"


void initbfgsb(int _minix);
int bfgsb(int,double*,int*,double*,double*);

