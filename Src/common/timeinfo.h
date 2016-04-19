#ifndef  _TIMEINFO_H
#define _TIMEINFO_H

int GetTimeInfo(int user,Planning *pt_plan,TimeInfo *met);
int ShowTimeInfo(int user,Planning *pt_plan,TimeInfo *met);
int ShowResultTimeInfo(int user,Planning *pt_plan,int error,TimeInfo *met);

#endif

