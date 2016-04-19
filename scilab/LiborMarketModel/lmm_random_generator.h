
#ifndef LMM_RANDOM_GENERATOR_H 
#define LMM_RANDOM_GENERATOR_H

#include"lmm_header.h"

int mallocRandomGenerator(RandomGenerator **ptRand ,int numOfFac );
int freeRandomGenerator(RandomGenerator **ptRand);
int randomVector(RandomGenerator *ptRand);
int printRandomGenerator(RandomGenerator *ptRand);
double getRandom(RandomGenerator *ptRand, int factorNumber);

#endif  




