CC0=gcc
CFLAGS=-g 

LINK=-lm 

MAINDIR=./
SUBDIR = ./premia_files
p6=$(PREMIADIR)

LMM_TEST=lmm_stochastic_volatility_test lmm_martingaleXtest \
lmm_martingaleVtest lmm_bermuda_swaption_test 
LMM_EXAM=lmm_stochastic_volatility_example lmm_martingaleX_example \
lmm_martingaleV_example


arbitrageFreeV=Vtest
arbitrageFreeX=Xtest

###### BGM Stochastic volatility
BGMSTO_C =  ${MAINDIR}lmm_stochastic_volatility.c ${MAINDIR}lmm_libor.c ${MAINDIR}lmm_random_generator.c ${MAINDIR}lmm_products.c  ${MAINDIR}lmm_volatility.c  ${MAINDIR}lmm_mathtools.c
BGMSTO_O = ${BGMSTO_C:.c=.o}
BGMSTO_H = ${BGMSTO_C:.c=.h}


###### ARBITRAGE FREE 
ARBITRAGEFREE_V_SRC = lmm_martingaleV.c  lmm_libor.c lmm_random_generator.c lmm_products.c  lmm_volatility.c  lmm_mathtools.c lmm_numerical.c cumulfunc.c  lmm_zero_bond.c ipmpar.c dcdflib.c 
ARBITRAGEFREE_V_O = ${ARBITRAGEFREE_V_SRC:.c=.o}
ARBITRAGEFREE_V_H = ${ARBITRAGEFREE_V_SRC:.c=.h}

ARBITRAGEFREE_X_SRC = lmm_martingaleX.c lmm_libor.c lmm_random_generator.c lmm_products.c  lmm_volatility.c  lmm_mathtools.c lmm_numerical.c cumulfunc.c lmm_zero_bond.c ipmpar.c dcdflib.c  
ARBITRAGEFREE_X_O = ${ARBITRAGEFREE_X_SRC:.c=.o}
ARBITRAGEFREE_X_H = ${ARBITRAGEFREE_X_SRC:.c=.h}



##### BERMUDA 
COMMON_SRC=${MAINDIR}lmm_mathtools.c ${MAINDIR}lmm_volatility.c ${MAINDIR}lmm_libor.c ${MAINDIR}lmm_random_generator.c ${MAINDIR}lmm_numerical.c ${MAINDIR}lmm_products.c ${MAINDIR}lmm_zero_bond.c  
COMMON_O=${COMMON_SRC:.c=.o}
COMMON_H=${COMMON_SRC:.c=.h}

BERMUDA_SRC=${MAINDIR}lmm_basis.c ${MAINDIR}lmm_bermudaprice.c
BERMUDA_O=${BERMUDA_SRC:.c=.o}
BERMUDA_H=${BERMUDA_SRC:.c=.h}




####### MATHTOOLS
MATHTOOLS_SRC = ${SUBDIR}/complex.c  
MATHTOOLS_H=${MATHTOOLS_SRC:.c=.h}
MATHTOOLS_O=${MATHTOOLS_SRC:.c=.o}



LIBS=-L/usr/lib/ -lm -lc

INC1=-I./
INCS=-I$(SUBDIR)/

CC=$(CC0)  $(INC1)  $(INCS)


%.o : %.c
	$(CC) -c $(CFLAGS) $< -o ${<:%.c=%.o}

all: $(LMM_TEST) $(LMM_EXAM) lmm_swpt_stovol_sci.o lmm_martingaleX_sci.o lmm_martingaleV_sci.o lmm_bermuda_LS_sci.o



#######bermudan pedersen

lmm_bermuda_swaption_test: ${COMMON_O} ${BERMUDA_O} $(MATHTOOLS_O) lmm_bermudatest.o 
	$(CC) $^ -o  lmm_bermuda_swaption_test $(LINK)

lmm_bermuda_example: ${COMMON_O} ${BERMUDA_O} $(MATHTOOLS_O)  lmm_bermuda_example.o 
	$(CC) $^ -o  lmm_bermuda_example   ${LINK}

lmm_bermudatest.o:${COMMON_H} ${BERMUDA_H}  lmm_bermudatest.c
lmm_bermuda_example.o:${COMMON_H} ${BERMUDA_H}  lmm_bermuda_example.c
lmm_bermudaprice.o:${COMMON_H} lmm_basis.h lmm_bermudaprice.c
lmm_basis.o: lmm_basis.c 
lmm_bermuda_LS_sci.o: ${COMMON_H} ${BERMUDA_H}  lmm_bermuda_LS_sci.c

##### stochastic volatility

lmm_stochastic_volatility_test:  $(BGMSTO_O) $(MATHTOOLS_O)  lmm_stochastic_volatility_test.o
	$(CC) $(LIBS) $^  -o lmm_stochastic_volatility_test 


lmm_stochastic_volatility_example:  $(BGMSTO_O) $(MATHTOOLS_O) lmm_stochastic_volatility_example.o
	$(CC) $(LIBS) $^  -o lmm_stochastic_volatility_example


lmm_stochastic_volatility_test.o: $(BGMSTO_H) $(MATHTOOLS_H)  lmm_stochastic_volatility_test.c
lmm_stochastic_volatility_example.o: $(BGMSTO_H) $(MATHTOOLS_H) lmm_stochastic_volatility_example.c
lmm_swpt_stovol_sci.o: $(BGMSTO_H) $(MATHTOOLS_H)  lmm_swpt_stovol_sci.c


###### Arbitrage free 


lmm_martingaleX_example: $(ARBITRAGEFREE_X_O)      lmm_martingaleX_example.o
	$(CC) $(LIBS) $^ -o lmm_martingaleX_example 

lmm_martingaleXtest: $(ARBITRAGEFREE_X_O)   lmm_martingaleXtest.o
	$(CC) $(LIBS) $^ -o lmm_martingaleXtest 

lmm_martingaleV_example: $(ARBITRAGEFREE_V_O) lmm_martingaleV_example.o
	$(CC) $(LIBS) $^ -o lmm_martingaleV_example 


lmm_martingaleVtest: $(ARBITRAGEFREE_V_O)    lmm_martingaleVtest.o
	$(CC) $(LIBS) $^  -o lmm_martingaleVtest

lmm_martingaleX_example.o:  $(ARBITRAGEFREE_X_H)   lmm_martingaleX_example.c
lmm_martingaleXtest.o:  $(ARBITRAGEFREE_X_H)    lmm_martingaleXtest.c
lmm_martingaleX_sci.o: $(ARBITRAGEFREE_X_H)     lmm_martingaleX_sci.c
lmm_martingaleV_example.o: $(ARBITRAGEFREE_V_H)   lmm_martingaleV_example.c
lmm_martingaleVtest.o:  $(ARBITRAGEFREE_V_H)   lmm_martingaleVtest.c
lmm_martingaleV_sci.o: $(ARBITRAGEFREE_V_H)     lmm_martingaleV_sci.c

clean:
	rm -f   $(OBJS0) $(OBJS1) $(OBJS2)  *~

cleanall:
	rm -f *.o *~ $(LMM_TEST) $(LMM_EXAM)


