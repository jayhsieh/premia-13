link('/usr/lib/libfftw3.so');
ilib_name='libnsp_cdo';
files=['nsp_cdo.o'];
libs=['../src/.libs/libcdo','../demo/.libs/libdemo'];
table=['price_cdo','nsp_premia_cdo'];
ilib_build(ilib_name,table,files,libs);

