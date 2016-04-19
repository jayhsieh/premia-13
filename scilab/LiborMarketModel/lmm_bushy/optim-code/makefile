LDFLAGS =       -lm
CC =    gcc

BIN =   optim

main=	main.c

#SRCS =  QuasiNewton.c BFGSupdate.c lineSearch.c cholesky.c  stopping.c testing.c
SRCS :=	$(wildcard *.c)

#HEADERS = $(SRCS:%.c=%.h)
HEADERS:= $(wildcard *.h) 

OBJECTS  = $(SRCS:%.c=%.o)
#	   $(main:%.c=%.o)

%.o :   %.c $(HEADERS)
	$(CC) -c  $<

all: $(OBJECTS)
	$(CC)  -o $(BIN) $(OBJECTS)  $(LDFLAGS) 
	$(clean)

clean:
	rm *.o
