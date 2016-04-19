#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <float.h>
#include "mt19937.h"
#include <math.h>
#include "../src/cdo.h"

int			main(void)
{
    copula		*cop;
    double		y;

    cop = init_nig_copula(0.30, 1.2558, -0.2231);
    for (y = - 6.; y < 6.; y+=0.1) {
//	printf("%g\t%g\n", y, cop->density(cop, y));
    }
    return (0);
}

