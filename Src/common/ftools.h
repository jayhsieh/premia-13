#ifndef  _FTOOLS_H
#define _FTOOLS_H
#define MAX_LINE 100
#define MAX_CHAR_LINE 200

int   FGetTimeInfo(char InputFile[MAX_LINE][MAX_CHAR_LINE],int user,Planning *pt_plan,TimeInfo *Met);
int StrCasecmp(const char *chaine1,const char *chaine2);

void ReadInputFile(char *InputFile, char RedFile[MAX_LINE][MAX_CHAR_LINE]);
char FChooseAction(char InputFile[MAX_LINE][MAX_CHAR_LINE]);
char FChooseProduct(char InputFile[MAX_LINE][MAX_CHAR_LINE]);
int FMoreAction(char InputFile[MAX_LINE][MAX_CHAR_LINE],int *count);

int FSelectModel(char InputFile[MAX_LINE][MAX_CHAR_LINE],int user,Planning *pt_plan,
                 Model * *listmod,Model **mod);
int FSelectOption(char InputFile[MAX_LINE][MAX_CHAR_LINE],int user,Planning *pt_plan,
                  Family **listopt,Model* pt_model,Pricing **pricing,Option **opt);

int FSelectPricing(char InputFile[MAX_LINE][MAX_CHAR_LINE],int user,Model *pt_model,
                   Option *pt_option,Pricing **family,Pricing **result);
int FSelectMethod(char InputFile[MAX_LINE][MAX_CHAR_LINE],int user,Planning *pt_plan,
                  Pricing *pt_pricing, Option *opt,Model *mod,PricingMethod **met);
int FSelectTest(char InputFile[MAX_LINE][MAX_CHAR_LINE],int user, Planning *pt_plan,
                Pricing *pt_pricing, Option *opt,Model *mod,PricingMethod *met,DynamicTest **test);

#endif

