INCLUDES = @INCLUDES@ -I$(top_srcdir)/src 
DEFS = @DEFS@ 
noinst_PROGRAMS = test test_bessel

test_SOURCES = test.c
test_LDADD = ../src/.libs/libcdo.a -lm -lfftw3

test_bessel_SOURCES = test_bessel.c
test_bessel_LDADD = ../src/.libs/libcdo.a -lm -lfftw3
