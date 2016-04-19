CC0=gcc
CFLAGS=-Wall -g 
#CFLAGS=-O3

p6=./from_premia
f77=f95
incf=$(p6)

svj=calib
fichier=$(svj) 

SVJ =  $(svj).c svj_std/ft_call_heston.c svj_std/ft_put_heston.c svj_std/ft_call_merton.c svj_std/ft_put_merton.c svj_std/ft_call_hestmert.c svj_std/ft_put_hestmert.c svj_std/svj.c svj_std/grad_svj.c svj_std/grad_ft_call_merton.c svj_std/grad_ft_call_heston.c svj_std/grad_ft_put_merton.c svj_std/grad_ft_put_heston.c  svj_std/grad_ft_call_hestmert.c svj_std/grad_ft_put_hestmert.c 

SRC0 = $(p6)/mathtools.c $(p6)/complex.c  ./utils/my_integral.c
SRC1 =  costFunction.c gradFunction.c   ./utils/datamarket.c  bfgsb.c ./utils/minix.c ./from_premia/routines.c
SRC2 =  parameter.c ./utils/plot.c ./utils/bs.c ./utils/utils1.c

LIBS1=-L/usr/lib/ -lm  -L/home/usr/lib/ -lg2c  

INC1=-I./
INC2=-I$(p6)

OBJSVJ=$(patsubst %.c,%.o,$(SVJ))
OBJS0=$(patsubst %.c,%.o,$(SRC0))
OBJS1=$(patsubst %.c,%.o,$(SRC1))
OBJS2=$(patsubst %.c,%.o,$(SRC2))

CC=$(CC0)  $(INC1) $(INC2)



all: $(fichier) 

$(svj): $(OBJSVJ) $(OBJS0) $(OBJS1) $(OBJS2) *.h svj_std/*.h 
	$(CC)  $(LIBS1) $(OBJSVJ) $(OBJS0) $(OBJS1) $(OBJS2) -o $(svj) 

# routines.o: $(incf)/routines.f
# 	$(f77) -c $(incf)/routines.f -o routines.o 



clean:
	echo "============================================"
	echo  
	echo "============================================"
	rm -f  $(OBJSVJ) $(OBJS0) $(OBJS1) $(OBJS2)  routines.o *~

veryclean:
	rm -f *.o *~ $(fichier)


