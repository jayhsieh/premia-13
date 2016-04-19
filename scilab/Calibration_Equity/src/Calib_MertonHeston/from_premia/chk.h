#ifndef _CHK_H
#define _CHK_H

int CHK_ok(int user, Planning *pt_plan,void*);

int CHK_call(int user, Planning *pt_plan,void*);
int CHK_put(int user, Planning *pt_plan,void*);
int CHK_callspread(int user, Planning *pt_plan,void*);
int CHK_digit(int user, Planning *pt_plan,void*);

int CHK_tree(int user, Planning *pt_plan,void*);
int CHK_mc(int user, Planning *pt_plan,void*);
int CHK_mcBaldi(int user, Planning *pt_plan,void*);
int CHK_fdiff(int user, Planning *pt_plan,void*);
int CHK_split(int user, Planning *pt_plan,void*);
int CHK_psor(int user, Planning *pt_plan,void*);
int CHK_mc_generator(int user, Planning *pt_plan,void*);

#endif




