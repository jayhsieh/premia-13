CC0=gcc

CFLAGS=-g

LINK=-lm

MAINDIR=./

SOURCE=${MAINDIR}lmm_cev_pricer.c ${MAINDIR}lmm_cev_sci.c
TEST=${SOURCE} ${MAINDIR}lmm_cev_test.c

OBJ=${SOURCE:.c=.o}
OBJTEST=${TEST:.c=.o}


test:${OBJTEST}
	${CC0} -o test_cev_pricer ${OBJTEST} ${LINK}


test_obj:${TEST} 
	${CC0} ${CFLAGS} -c ${TEST}


scilab:${SOURCE}
	${CC0} ${CFLAGS} -c ${SOURCE}



