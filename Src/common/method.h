#ifndef  _METHOD_H
#define _METHOD_H

int GetMethod(int user,Planning *pt_plan,Pricing *Pr,PricingMethod *met,Option *opt);
int FGetMethod(char **InputFile,int user,Planning *pt_plan,Pricing *Pr,PricingMethod *met,Option *opt);
int ShowMethod(int user,Planning *pt_plan,Pricing *Pr,PricingMethod *met,Option *opt);
int ShowResultMethod(int user,Planning *pt_plan,int error,PricingMethod *met);
int FixMethod(Planning *pt_plan, PricingMethod *pt_method);
int CompareParameterNames(const char *s1, const char *s2);
#endif

