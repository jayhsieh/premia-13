INCLUDES = @INCLUDES@ -I$(top_srcdir)/src 
DEFS = @DEFS@

lib_LTLIBRARIES = libcdo.la
libcdo_la_SOURCES = mt19937.h \
		   mt19937.c \
		   maths.h \
		   maths.c \
		   structs.h \
		   structs.c \
		   saddlepoint.c\
		   Implied_Copula.c \
		   Base_Correlation.c \
		   copulas.h \
		   gaussian_copula.c \
		   clayton_copula.c \
		   student.c \
		   double_t_copula.c \
		   nig_copula.c \
		   company.h \
		   company.c \
		   cdo.h \
		   cdo.c \
		   hull_white.c \
		   laurent_gregory.c \
		   monte_carlo.c 

libcdo_la_CFLAGS="-std=c99"
