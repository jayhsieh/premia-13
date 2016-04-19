#ifndef _basis_H_
#define _basis_H_

extern void Name_To_Basis(char*, char*, double (**)(double*, int),int Basis_Dimension);
extern double Tensor_Delta(int j,int k, int BSdim,int *power);

#endif
