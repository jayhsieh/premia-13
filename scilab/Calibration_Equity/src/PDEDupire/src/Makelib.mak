# generated by builder.sce: Please do not edit this file
# ------------------------------------------------------
SHELL = /bin/sh
SCIDIR =C:/Scilab-2.7.2
SCIDIR1 =C:\Scilab-2.7.2
# name of the dll to be built
LIBRARY = libcalib
# list of objects file
OBJS = dcalib.obj calib.obj objFunction.obj gradFiniteDiff.obj spline.obj inout.obj BFGSupdate.obj cholesky.obj sparse.obj DupirePDE.obj routines.obj QuasiNewton.obj lineSearch.obj stopping.obj solveSystem.obj libcalib.obj
# added libraries 
OTHERLIBS = 
!include $(SCIDIR1)\Makefile.incl.mak
CFLAGS = $(CC_OPTIONS) -DFORDLL -I"$(SCIDIR)/routines" -Dmexfunction_=mex$*_  -DmexFunction=mex_$*  
FFLAGS = $(FC_OPTIONS) -DFORDLL -I"$(SCIDIR)/routines" -Dmexfunction=mex$* 
EXTRA_LDFLAGS = 
!include $(SCIDIR1)\config\Makedll.incl 
