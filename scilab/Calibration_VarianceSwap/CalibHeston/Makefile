
BINS=CalibHeston
CalibHeston_OBJS = calib_heston.o calibration1.o utils/intg.o utils/ctcab.o utils/ctonb.o \
	costFunction.o gradFunction.o utils/dcube.o utils/fcube.o \
	utils/my_integral.o utils1.o varswaps.o bsvanillas.o utils/fpq2.o  \
	utils/fuclid.o utils/gcbd.o  utils/majour.o utils/majysa.o utils/n1fc1.o \
	utils/n1gc2.o utils/n1qn1.o \
	utils/n1qn2.o utils/n1qn3.o utils/n2qn1.o utils/nlis0.o  \
	utils/proj.o utils/qnbd.o  utils/rlbd.o  


include ../../../Src/libpnl/Makeuser.incl


clean::
	$(RM) utils/*.o utils/*.lo *~
