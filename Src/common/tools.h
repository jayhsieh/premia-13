#ifndef  _TOOLS_H
#define _TOOLS_H

void premia_spawnlp(const char *fhelp_name);
int Iterate(Planning *pt_plan,Iterator* pt_iterator,int count,char action,
            Model *pt_model,Option *pt_option,Pricing *pt_pricing,PricingMethod *pt_method, 
            DynamicTest *pt_test,int user,TimeInfo *pt_time_info);
void Action(Model *pt_model,Option *pt_option,Pricing *pt_pricing,PricingMethod *pt_method, 
            DynamicTest *pt_test,int user,Planning *pt_plan,TimeInfo *pt_time_info);

int OutputFile(FILE **pt_file);
void InputMode(int *pt_user);

char ChooseAction(char);
char ChooseProduct(void);
int MoreAction(int *count);

int SelectModel(int user,Planning *pt_plan,Model * *listmod,Family **families, Pricing **pricings,Model **mod);
int Premia_model_has_products (Model *pt_model, Family **families, Pricing **pricings);
int SelectOption(int user,Planning *pt_plan,Family **listopt,Model* pt_model,Pricing **pricing,Option **opt);
int MatchingPricing(Model *pt_model,Option *pt_option,Pricing **pricing);
int SelectPricing(int user,Model *pt_model,Option *pt_option,Pricing **family,Pricing **result);
int SelectMethod(int user,Planning *pt_plan,Pricing *pt_pricing, Option *opt,Model *mod,PricingMethod **met);
int SelectTest(int user,Planning *pt_plan,Pricing *pt_pricing, Option *opt,Model *mod,PricingMethod *met,DynamicTest **test);

void BuildGnuStuff(Planning *plan,Model* pt_model,Option *pt_option,Pricing *pt_pricing,PricingMethod **pt_method);
int  make_titles_file(char *name,Model* pt_model,Option* pt_option, Pricing *pt_pricing, PricingMethod **pt_method, int metno);
void begin_tex_file(FILE *f_tex);
void end_tex_file(FILE *f_tex);
void begin_gnu(FILE *fp);
void FurtherMsg(void);

void WellcomeMsg(int user);
int NextSession(Planning *pt_plan,char action,int user);


void FreeTest(DynamicTest *test);
void BuildGnuStuffTest(Model* model,Option* option, Pricing* pricing, PricingMethod* method ,DynamicTest* test);
int Premia_match_model_option (Model *pt_model, Option *pt_opt, Pricing **pricing);

void get_model_helpfile (Model *mod, char *helpfile);
void get_option_helpfile (Option *opt, char *helpfile);
void get_method_helpfile (Pricing *Pr, PricingMethod *Met, char *helpfile);
void get_method_helpfile_with_ids (PricingMethod *Met, const char *mod_id, const char *opt_id, char *helpfile);

#endif

