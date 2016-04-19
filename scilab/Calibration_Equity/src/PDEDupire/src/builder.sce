ilib_name="libcalib"

files=["dcalib.o";"calib.o";'objFunction.o';'gradFiniteDiff.o';'spline.o';'inout.o';'BFGSupdate.o';'cholesky.o';'sparse.o';'DupirePDE.o';'routines.o';'QuasiNewton.o';'lineSearch.o';'stopping.o';'solveSystem.o'];

libs=[];

table=["calib","dcalib"];




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










