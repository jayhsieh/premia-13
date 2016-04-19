CC0=gcc
CFLAGS=-Wall 

p6=../from_premia


SVJ =  ft_call_heston.c ft_put_heston.c ft_call_merton.c ft_put_merton.c ft_call_hestmert.c ft_put_hestmert.c svj.c grad_svj.c grad_ft_call_merton.c grad_ft_call_heston.c grad_ft_put_merton.c grad_ft_put_heston.c grad_ft_call_hestmert.c grad_ft_put_hestmert.c

LIBS=-L/usr/lib/ -lm 

INC1=-I./
INC2=-I$(p6)



OBJS=$(patsubst %.c,%.o,$(SVJ))

CC=$(CC0)  $(INC1) $(INC2)



all: $(OBJS) 

doc:	svj_std.pdf

svj_std.pdf: svj_std.dvi
	dvipdf svj_std.dvi 
svj_std.dvi: svj_std.tex
	latex svj_std.tex

clean:
	rm -f *.o *~ 

veryclean:
	rm -f *.o *~ $(fichier)


