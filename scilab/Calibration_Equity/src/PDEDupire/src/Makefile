# -------------------------------------------------------------------------
# Les compilateurs a utiliser:

CC      = gcc

# -------------------------------------------------------------------------
# le nom du programme a compiler, votre fichier .c doit avoir la meme
# premiere partie de nom, par exemple si vous voulez que votre programme 
# s'appelle "bonjour", il faut que le code source s'appelle "bonjour.c"
# Rem: ceci est du a la configuration de ce Makefile, vous pouvez le
# modifier si vous voulez.

PROGRAM1 = simul
PROGRAM2 = calib1
PROGRAM3 = rafsig
PROGRAM4 = visusig
PROGRAM5 = impsig

# -------------------------------------------------------------------------
# une liste d'options de compilations

CFLAGS = -g
FFLAGS = -O

# -------------------------------------------------------------------------
# une liste de librairies a inclure pour l'edition des liens (m est la
# librairie mathematique et c la librairie standard de C)

LIBS   = -lm -lc -lg2c

# -------------------------------------------------------------------------
# la liste d'objets a compiler, ie. au moins le programme principal:

OBJECTS1  =   simul.o DupirePDE.o solveSystem.o spline.o sparse.o inout.o
OBJECTS2  =   calib1.o DupirePDE.o solveSystem.o spline.o gradFiniteDiff.o sparse.o objFunction.o testing.o stopping.o QuasiNewton.o lineSearch.o cholesky.o BFGSupdate.o inout.o routines.o
OBJECTS3  =   rafsig.o DupirePDE.o solveSystem.o spline.o sparse.o inout.o
OBJECTS4  =   visusig.o solveSystem.o spline.o inout.o
OBJECTS5  =   impsig.o DupirePDE.o inout.o solveSystem.o spline.o

# -------------------------------------------------------------------------
# la maniere dont les objets (*.o) doivent se compiler (edition des liens)

all: $(PROGRAM1) $(PROGRAM2) $(PROGRAM3) $(PROGRAM4) $(PROGRAM5)

$(PROGRAM1): $(OBJECTS1)
	$(CC) $(CFLAGS) -o $(PROGRAM1) $(OBJECTS1) $(LIBS)
$(PROGRAM2): $(OBJECTS2)
	$(CC) $(CFLAGS) -o $(PROGRAM2) $(OBJECTS2) $(LIBS)
$(PROGRAM3): $(OBJECTS3)
	$(CC) $(CFLAGS) -o $(PROGRAM3) $(OBJECTS3) $(LIBS)
$(PROGRAM4): $(OBJECTS4)
	$(CC) $(CFLAGS) -o $(PROGRAM4) $(OBJECTS4) $(LIBS)
$(PROGRAM5): $(OBJECTS5)
	$(CC) $(CFLAGS) -o $(PROGRAM5) $(OBJECTS5) $(LIBS)

# -------------------------------------------------------------------------
# la maniere dont les sources (*.c) doivent se compiler

simul.o: simul.c DupirePDE.h spline.h
	$(CC) $(CFLAGS) -c simul.c

calib.o: calib1.c objFunction.h gradFiniteDiff.h spline.h inout.h
	$(CC) $(CFLAGS) -c calib1.c 

rafsig.o: rafsig.c DupirePDE.h solveSystem.h spline.h sparse.h inout.h
	$(CC) $(CFLAGS) -c rafsig.c

visusig.o: visusig.c solveSystem.h spline.h inout.h
	$(CC) $(CFLAGS) -c visusig.c

impsig.o: impsig.c DupirePDE.h inout.h
	$(CC) $(CFLAGS) -c impsig.c

DupirePDE.o: DupirePDE.c DupirePDE.h solveSystem.h
	$(CC) $(CFLAGS) -c DupirePDE.c

solveSytem.o: solveSystem.c solveSystem.h
	$(CC) $(CFLAGS) -c solveSystem.c

spline.o: spline.c spline.h solveSystem.h
	$(CC) $(CFLAGS) -c spline.c

gradFiniteDiff.o: gradFiniteDiff.c gradFiniteDiff.h data.h spline.h DupirePDE.h objFunction.h
	$(CC) $(CFLAGS) -c gradFiniteDiff.c

grad.o: grad.c grad.h data.h spline.h DupirePDE.h sparse.h utilsGrad.h solveSystem.h 
	$(CC) $(CFLAGS) -c grad.c

utilsGrad.o: utilsGrad.c utilsGrad.h solveSystem.h data.h sparse.h DupirePDE.h spline.h
	$(CC) $(CFLAGS) -c utilsGrad.c

sparse.o: sparse.c sparse.h solveSystem.h
	$(CC) $(CFLAGS) -c sparse.c

objFunction.o: objFunction.c objFunction.h data.h spline.h DupirePDE.h
	$(CC) $(CFLAGS) -c objFunction.c

inout.o: inout.c inout.h data.h spline.h
	$(CC) $(CFLAGS) -c inout.c

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

routines.o: routines.c
	$(CC) $(FFLAGS) -c routines.c

# -------------------------------------------------------------------------
# comment faire une edition des lien incluant l'outil purify (l'executable
# contenant purify sera nomme avec pure_<mon_programme>).

pure1: $(OBJECTS1) 
	/usr/Rational/releases/purify-5.3-solaris2/purify -cache-dir=/afs/inria.fr/rocq/home/mathfi/cognet/purify/cache $(CC) $(CFLAGS) -o pure_$(PROGRAM1) $(OBJECTS1) $(LIBS)

pure2: $(OBJECTS2) 
	/usr/Rational/releases/purify-5.3-solaris2/purify -cache-dir=/afs/inria.fr/rocq/home/mathfi/cognet/purify/cache $(CC) $(CFLAGS) -o pure_$(PROGRAM2) $(OBJECTS2) $(LIBS)

pure3: $(OBJECTS3) 
	/usr/Rational/releases/purify-5.3-solaris2/purify -cache-dir=/afs/inria.fr/rocq/home/mathfi/cognet/purify/cache $(CC) $(CFLAGS) -o pure_$(PROGRAM3) $(OBJECTS3) $(LIBS)

pure4: $(OBJECTS4) 
	/usr/Rational/releases/purify-5.3-solaris2/purify -cache-dir=/afs/inria.fr/rocq/home/mathfi/cognet/purify/cache $(CC) $(CFLAGS) -o pure_$(PROGRAM4) $(OBJECTS4) $(LIBS)

pure5: $(OBJECTS5) 
	/usr/Rational/releases/purify-5.3-solaris2/purify -cache-dir=/afs/inria.fr/rocq/home/mathfi/cognet/purify/cache $(CC) $(CFLAGS) -o pure_$(PROGRAM5) $(OBJECTS5) $(LIBS)

# -------------------------------------------------------------------------
# une cible pour nettoyer le repertoire de compilation du programme
# avec "make clean" tous les fichiers temporaires sont effaces, sauf
# le code source que vous avez ecrit, bien sur.

clean:
	rm -f *.o *~ $(PROGRAM1) $(PROGRAM2) $(PROGRAM3) $(PROGRAM4) $(PROGRAM5) pure_$(PROGRAM1) pure_$(PROGRAM2) pure_$(PROGRAM3) pure_$(PROGRAM4) pure_$(PROGRAM5) core *.lo libcalib.a libcalib.la libcalib.so 
