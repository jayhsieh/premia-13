//link('/opt/local/lib/libfftw3.dylib');
ilib_name='libsci_cdo';
files=['sci_cdo.o'];
libs=['../src/.libs/libcdo','../demo/.libs/libdemo'];
table=['price_cdo','sci_cdo'];
ilib_build(ilib_name,table,files,libs);

