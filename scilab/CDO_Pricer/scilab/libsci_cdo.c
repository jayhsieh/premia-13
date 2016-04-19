#include <mex.h> 
static int direct_gateway(char *fname,void F(void)) { F();return 0;};
extern Gatefunc sci_cdo;
static GenericTable Tab[]={
  {(Myinterfun)sci_gateway,sci_cdo,"price_cdo"},
};
 
int C2F(libsci_cdo)()
{
  Rhs = Max(0, Rhs);
  (*(Tab[Fin-1].f))(Tab[Fin-1].name,Tab[Fin-1].F);
  return 0;
}
