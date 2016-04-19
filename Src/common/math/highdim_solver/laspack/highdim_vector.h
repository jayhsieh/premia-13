/****************************************************************************/
/*                                 highdim_vector.h                                 */
/****************************************************************************/
/*                                                                          */
/* type VECTOR                                                              */
/*                                                                          */
/* Copyright (C) 1992-1995 Tomas Skalicky. All rights reserved.             */
/*                                                                          */
/****************************************************************************/
/*                                                                          */
/*        ANY USE OF THIS CODE CONSTITUTES ACCEPTANCE OF THE TERMS          */
/*              OF THE COPYRIGHT NOTICE (SEE FILE COPYRGHT.H)               */
/*                                                                          */
/****************************************************************************/

#ifndef HIGHDIM_VECTOR_H
#define HIGHDIM_VECTOR_H

#include <stdlib.h>

#include "lastypes.h"
#include "elcmp.h"
#include "copyrght.h"

typedef struct {
    char *Name;
    size_t Dim;
    InstanceType Instance;
    int LockLevel;
    double Multipl;
    Boolean OwnData;
    double *Cmp;
} Vector;

void V_Constr(Vector *V, char *Name, size_t Dim, InstanceType Instance,
	      Boolean OwnData);
void V_Destr(Vector *V);
void V_SetName(Vector *V, char *Name);
char *V_GetName(Vector *V);
size_t V_GetDim(Vector *V);
void V_SetCmp(Vector *V, size_t Ind, double Val);
void V_SetAllCmp(Vector *V, double Val);
void V_SetRndCmp(Vector *V);
double V_GetCmp(Vector *V, size_t Ind);
void V_AddCmp(Vector *V, size_t Ind, double Val);

/* macros for fast access */
#define     V__SetCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] = (Val)
#define     V__GetCmp(PtrV, Ind)            (PtrV)->Cmp[Ind]
#define     V__AddCmp(PtrV, Ind, Val)       (PtrV)->Cmp[Ind] += (Val)

void V_Lock(Vector *V);
void V_Unlock(Vector *V);

#endif /* HIGHDIM_VECTOR_H */

