SRC_C = main.c maths.c  tryschon.c spread_dynamique.c initial_calibration.c spread_calib.c
OBJ_C = ${SRC_C:%.c=%.o}

.SUFFIXES: .c.o

.c.o: gcc -c $< -o $@


all:: $(OBJ_C) $(OBJ_CPP)
	gcc -lm  $(OBJ_C) -o calibration

clean::
	rm *.o 
#	find -name .\ex -print -exec rm {} \
 
