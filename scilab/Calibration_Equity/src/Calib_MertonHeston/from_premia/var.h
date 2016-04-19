#ifndef  _VAR_H
#define _VAR_H
#include "ftools.h"
int Fprintf(int user,const char *s,...);
int FprintfVar(int user,const char s[],VAR *x);
int FScanVar(char **InputFile,Planning *pt_plan,int user,VAR *x);

int Valid (int user,int status,char *helpfile);
int FGetParVar(char **InputFile,Planning *pt_plan,int user,VAR *x);

int InitVar(void);
void ExitVar(void);
int PrintVar(Planning *pt_plan,int user,VAR*);
int ScanVar(Planning *pt_plan,int user,VAR*);
int ChkVar(Planning *pt_plan,VAR *x);
int GetParVar(Planning *pt_plan,int user,VAR *x);
int ShowParVar(Planning *pt_plan,int user,VAR *x);
int ChkParVar(Planning *pt_plan,VAR *x);
int LowerVar(int user,VAR *x, VAR*y);
void CopyVar(VAR *srce,VAR *dest);
int CheckIterationValue(char **InputFile,Planning *pt_plan,int user,VAR *x,int nbline,int nbchar);

void ResetPlanning(Planning *pt_plan);
void ShowPlanning(int user,Planning *pt_plan);
void ShrinkPlanning(int index,Planning *pt_plan);
int ChkStepNumber(int user,Iterator *pt_iterator,int step);
void NextValue(int count,Iterator* pt_iterator);


int ShowParVarTestRes(Planning *pt_plan, int user, VAR *x);

#ifndef _WIN32
int Spawnlp( int mode, const char *cmdname, const char *arg0, const char *arg1, const char *arg2 );
#endif



#endif

