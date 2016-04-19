# -------------------------------------------------------------------------
# Les compilateurs a utiliser:

CC      = gcc -ggdb

# -------------------------------------------------------------------------
# le nom du programme a compiler, votre fichier .c doit avoir la meme
# premiere partie de nom, par exemple si vous voulez que votre programme
# s'appelle "bonjour", il faut que le code source s'appelle "bonjour.c"
# Rem: ceci est du a la configuration de ce Makefile, vous pouvez le
# modifier si vous voulez.

PROGRAM1 = pricer_Avellaneda
PROGRAM2 = calib_Avellaneda

# -------------------------------------------------------------------------
# repertoire des fichiers d'optimisation (QUASINEWTON)

SUBDIRS = ./optim-code

# -------------------------------------------------------------------------
# une liste d'options de compilations

CFLAGS = -g
FFLAGS = -O

# -------------------------------------------------------------------------
# une liste de librairies a inclure pour l'edition des liens (m est la
# librairie mathematique et c la librairie standard de C)

LIBS   = -lm -lc

# -------------------------------------------------------------------------
# la liste d'objets a compiler, ie. au moins le programme principal:

OBJECTS1  =   pricer_Avellaneda.o inout.o tree.o
OBJECTS2  =   calib_Avellaneda.o tree.o inout.o rollback.o \
$(SUBDIRS)/testing.o $(SUBDIRS)/stopping.o $(SUBDIRS)/QuasiNewton.o \
$(SUBDIRS)/lineSearch.o $(SUBDIRS)/cholesky.o $(SUBDIRS)/BFGSupdate.o

# -------------------------------------------------------------------------
# la maniere dont les objets (*.o) doivent se compiler (edition des liens)

all: $(PROGRAM1) $(PROGRAM2)

$(PROGRAM1): $(OBJECTS1)
	@echo " "
	@echo "creation de l'executable" $(PROGRAM1)
	$(CC) $(CFLAGS) -o $(PROGRAM1) $(OBJECTS1) $(LIBS)
	@echo " "
$(PROGRAM2): $(OBJECTS2)
	@echo " "
	@echo "creation de l'executable" $(PROGRAM2)
	$(CC) $(CFLAGS) -o $(PROGRAM2) $(OBJECTS2) $(LIBS)
	@echo " "
# -------------------------------------------------------------------------
# la maniere dont les sources (*.c) doivent se compiler

tree.o: tree.c tree.h
	$(CC) $(CFLAGS) -c tree.c

rollback.o: rollback.c rollback.h tree.h
	$(CC) $(CFLAGS) -c rollback.c

inout.o: inout.c inout.h tree.h
	$(CC) $(CFLAGS) -c inout.c

calib_Avellaneda.o: calib_Avellaneda.c tree.h rollback.h inout.h ./optim-code/QuasiNewton.h
	$(CC) $(CFLAGS) -c calib_Avellaneda.c

pricer_Avellaneda.o: pricer_Avellaneda.c pricer_Avellaneda.h inout.h tree.h
	$(CC) $(CFLAGS) -c pricer_Avellaneda.c

testing.o: testing.c testing.h
	$(CC) $(CFLAGS) -c testing.c

stopping.o: stopping.c stopping.h
	$(CC) $(CFLAGS) -c stopping.c

QuasiNewton.o: QuasiNewton.c QuasiNewton.h lineSearch.h BFGSupdate.h cholesky.h stopping.h testing.h
	$(CC) $(CFLAGS) -c QuasiNewton.c

lineSearch.o: lineSearch.c lineSearch.h
	$(CC) $(CFLAGS) -c lineSearch.c

cholesky.o: cholesky.c cholesky.h
	$(CC) $(CFLAGS) -c cholesky.c

BFGSupdate.o: BFGSupdate.c BFGSupdate.h cholesky.h
	$(CC) $(CFLAGS) -c BFGSupdate.c

# -------------------------------------------------------------------------
# une cible pour nettoyer le repertoire de compilation du programme
# avec "make clean" tous les fichiers temporaires sont effaces, sauf
# le code source que vous avez ecrit, bien sur.

clean:
	rm -f *.o ./optim-code/*.o *~ $(PROGRAM1) $(PROGRAM2)
