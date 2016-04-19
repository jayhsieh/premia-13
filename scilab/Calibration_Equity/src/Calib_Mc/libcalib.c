#include <mex.h> 
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc dcalib;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,dcalib,"calib"},
};
 
int C2F(libcalib)()
{
  Rhs = Max(0, Rhs);
  (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  return 0;
}
