ilib_name="libcalib"

files=["dcalib.o";"calib_ave_scilab.o";"tree.o";"inout.o";"rollback.o";"./optim-code/testing.o";"./optim-code/stopping.o";"./optim-code/QuasiNewton.o";"./optim-code/lineSearch.o";"./optim-code/cholesky.o";"./optim-code/BFGSupdate.o"];

libs=[];

table=["calib_ave_scilab","dcalib"];




// extra parameters can be transmited to the linker 
// and to the C and Fortran compilers with 
// ldflags,cflags,fflags 
// for example to link a set of routines using the 
// ImageMagick library 
//  ldflags = "`Magick-config --ldflags --libs`"; 
//  cflags  = "`Magick-config --cppflags`"; 
//  fflags   = ""; 

ldflags = "";
cflags ="";
fflags ="";

// do not modify below 
// ----------------------------------------------
ilib_build(ilib_name,table,files,libs,'Makelib',ldflags,cflags,fflags)
