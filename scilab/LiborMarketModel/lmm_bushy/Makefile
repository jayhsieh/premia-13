# -------------------------------------------------------------------------
# Les compilateurs a utiliser:

CC      = gcc -ggdb

# -------------------------------------------------------------------------
# le nom du programme a compiler, votre fichier .c doit avoir la meme
# premiere partie de nom, par exemple si vous voulez que votre programme
# s'appelle "bonjour", il faut que le code source s'appelle "bonjour.c"
# Rem: ceci est du a la configuration de ce Makefile, vous pouvez le
# modifier si vous voulez.

PROGRAM1 = lmm_bushy_pricer

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

OBJECTS1  =   ManageMemory.o MMsimulBGM.o spline.o $(SUBDIRS)/testing.o $(SUBDIRS)/stopping.o $(SUBDIRS)/QuasiNewton.o \
$(SUBDIRS)/lineSearch.o $(SUBDIRS)/cholesky.o $(SUBDIRS)/BFGSupdate.o

# -------------------------------------------------------------------------
# la maniere dont les objets (*.o) doivent se compiler (edition des liens)

all: $(PROGRAM1)

$(PROGRAM1): $(OBJECTS1)
	@echo " "
	@echo "creation de l'executable" $(PROGRAM1)
	$(CC) $(CFLAGS) -o $(PROGRAM1) $(OBJECTS1) $(LIBS)
	@echo " "
# -------------------------------------------------------------------------
# la maniere dont les sources (*.c) doivent se compiler

ManageMemory.o: ManageMemory.c MMsimulBGM.h
	$(CC) $(CFLAGS) -c ManageMemory.c

MMsimulBGM.o: MMsimulBGM.c MMsimulBGM.h spline.h ./optim-code/QuasiNewton.h
	$(CC) $(CFLAGS) -c MMsimulBGM.c

spline.o: spline.c spline.h
	$(CC) $(CFLAGS) -c spline.c

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
