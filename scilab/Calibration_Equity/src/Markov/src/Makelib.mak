# generated by builder.sce: Please do not edit this file
# ------------------------------------------------------
SHELL = /bin/sh
SCIDIR =C:/Scilab-2.7.2
SCIDIR1 =C:\Scilab-2.7.2
# name of the dll to be built
LIBRARY = libmanu
# list of objects file
OBJS = implied.obj implied_volatility.obj libmanu.obj
# added libraries 
OTHERLIBS = 
!include $(SCIDIR1)\Makefile.incl.mak
CFLAGS = $(CC_OPTIONS) -DFORDLL -I"$(SCIDIR)/routines" -Dmexfunction_=mex$*_  -DmexFunction=mex_$*  
FFLAGS = $(FC_OPTIONS) -DFORDLL -I"$(SCIDIR)/routines" -Dmexfunction=mex$* 
EXTRA_LDFLAGS = 
!include $(SCIDIR1)\config\Makedll.incl 
