#pragma once

extern "C" {
#include "error_msg.h"
#include "optype.h"
#include "ftools.h"
//#include "random.h"
#include "var.h"
//#include "math\bsnd_math\basis.h"
#include "..\premia_obj.h"
//#include "premia_path.h"
//#include "matrix.h"

    extern Model BSND_model;
    extern Pricing BSND_STDND_pricing;
	extern Family STDND_family;
    extern FILE * out_stream;
    extern VAR g_printvararray[255];
    extern int g_printvararray_size;
    extern int g_dup_printf;
    extern FILE * g_dup_file;
    extern char premiamandir[256];
    extern int * true_typeV;

    extern PremiaAsset premia_assets[];

    int InitErrorMsg();
    int FGetMethod(char **InputFile,int user,Planning *pt_plan,Pricing *Pr,PricingMethod *Met);

    
}

#include <comutil.h>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <sstream>

#include "registry.h"

namespace prxl {
    /// contains premialib related stuff

namespace lib {
	
	typedef int Mode;

/// \param mode - a mode for which models_e are queried
/// \return a pointer to an array of models_e belonging to the mode
inline Model ** getModels(Mode mode)
{
    return premia_assets[mode].models;
}

/// \param mode - a mode for which pricings_e are queried
/// \return a pointer to an array of pricings_e belonging to the mode
inline Pricing ** getPricings(Mode mode)
{
    return premia_assets[mode].pricings;
}

/// \param mode - a mode for which families_e are queried
/// \return a pointer to an array of families_e belonging to the mode
inline Family ** getFamilies(Mode mode)
{
    return premia_assets[mode].families;
}

/// Searches for a model with given name
/// \param mode - a mode where models_e should be searched
/// \param name - a name for a model being searched
/// \return a pointer to a model instance if found, otherwise 0
inline Model * getModel(Mode mode, _bstr_t const & name)
{
    for (Model ** m_it = getModels(mode); *m_it; ++m_it)
        if (_bstr_t((**m_it).Name) == name)
            return *m_it;

    return 0;
}

/// Searches for a family with a given name
/// \param mode - a mode which families_e should be searched
/// \param familyName - a name for a family being searched
/// \return a pointer to a family instance if found, otherwise 0
inline Family * getFamily(Mode mode, _bstr_t const & familyName)
{
    for (Family ** f_it = getFamilies(mode); *f_it; ++f_it)
        if (_bstr_t((**f_it)[0]->ID) == familyName)
            return *f_it;

    return 0;
}

/// Searches for an option with specific name in a given family
/// \param optionName - name of an option being searched
/// \param family - a family where the option is searched
/// \return a pointer to an option instance if found, otherwise 0
inline Option * getOption(_bstr_t const & optionName, Family * family)
{
    for (Option ** opt_it = *family; *opt_it; ++opt_it)
        if (_bstr_t((**opt_it).Name) == optionName)
            return *opt_it;

    return 0;
}

/// Searches for an option with specific name in a given mode
/// \param mode - a mode where an option is searched
/// \param optionName - name of an option being searched
/// \param pFamily - output parameter for putting a pointer to a family of option found
/// \return a pointer to an option instance if found, otherwise 0
inline Option * getOption(Mode mode, _bstr_t const & optionName, Family ** pFamily = 0)
{
    for (Family ** f_it = getFamilies(mode); *f_it; ++f_it)
        if (Option * opt = getOption(optionName, *f_it))
        {
            if (pFamily)
                *pFamily = *f_it;
            return opt;
        }

        return 0;
}

/// Searches for an option with specific name in a given mode
/// \param mode - a mode where an option is searched
/// \param optionName - name of an option being searched
/// \return a pointer to an option instance if found, otherwise 0
inline Option * getOption_2(Mode mode, _bstr_t const & optionName)
{
    return getOption(mode, optionName);
}

/// searches for a family that contains the option
/// \param mode - a mode where to search
/// \param opt - an option
/// \return a pointer to a family instance if found, otherwise 0
inline Family * getFamily(Mode mode, Option * opt)
{
    for (Family ** f_it = getFamilies(mode); *f_it; ++f_it)
        if (0 == strcmp(opt->ID, (**f_it)[0]->ID))
            return *f_it;

    return 0;
}

/// returns a string containing a name of a pricing for given model and family
/// \param modelName - a model ID
/// \param familyName - a family ID
/// \return a pricing ID 
inline _bstr_t getPricingName(_bstr_t const & modelName, _bstr_t const & familyName)
{
    _bstr_t name (modelName);

    name += "_";
    name += familyName;

    return name;
}

/// returns a string containing a name of a pricing for given model and option
/// \param model - a model 
/// \param opt - an option
/// \return a pricing ID 
inline _bstr_t getPricingName(Model * model, Option * opt)
{
    return getPricingName(model->ID, opt->ID);
}

/// Searches for a pricing with given name
/// \param mode - a mode where a pricing is searched
/// \param name - a ID for a pricing to be found
/// \return a pointer to a pricing if found, otherwise 0
inline Pricing * getPricing(Mode mode, _bstr_t const & name)
{
    for (Pricing ** p_it = getPricings(mode); *p_it; ++p_it)
        if (_bstr_t((**p_it).ID) == name)
            return *p_it;

    return 0;
}

/// Searches for a pricing for a given model and option
/// \param mode - a mode where a pricing is searched
/// \param mod - a model
/// \param opt - an option
/// \return a pointer to a pricing if found, otherwise 0
inline Pricing * getPricing(Mode mode, Model * mod, Option * opt)
{
    return getPricing(mode, getPricingName(mod, opt));
}

/// Searches for a pricing containing the given pricing method
/// \param mode - a mode where a pricing is searched
/// \param method - a method of a pricing being searched
/// \return a pointer to a pricing if found, otherwise 0
inline Pricing * getPricing(Mode mode, PricingMethod * method)
{
    for (Pricing ** p_it = getPricings(mode); *p_it; ++p_it)
        for (PricingMethod ** pm_it = (**p_it).Methods; *pm_it; ++pm_it)
            if (*pm_it == method)
                return *p_it;

    return 0;
}

/// Searches for a pricing method with a given name compatible with a model and an option
/// \param mode - a mode where the method is searched
/// \param model - a model for the method
/// \param option - an option for the method
/// \param methodName - a name of a method being searched
/// \param pFoundButNotCompatable - output parameter which is set to true iff a method with a given name is found but not compatible with the model and the option
/// \return a pointer to a pricing if found, otherwise 0
inline PricingMethod * getMethod(Mode mode, Model * model, Option * option, _bstr_t const & methodName, 
                                 bool * pFoundButNotCompatable = 0)
{
    pFoundButNotCompatable && (*pFoundButNotCompatable = false);

    Pricing * pricing = getPricing(mode, model, option);

    model->Init(model);
    option->Init(option, model);

    for (PricingMethod ** pm_it = pricing->Methods; *pm_it; ++pm_it)
    {
        PricingMethod * pm = *pm_it;
        pm->Init(pm, option);

        if (_bstr_t(pm->Name) == methodName)
            if (pm->CheckOpt(option, model) == OK)
                return *pm_it;
            else
                pFoundButNotCompatable && (*pFoundButNotCompatable = true);
    }  

    return 0;
}

/// Searches for a pricing method with a given name compatible with a model and an option
/// \param mode - a mode where the method is searched
/// \param modelName - a name of a model for the method
/// \param optionName - a name of an option for the method
/// \param methodName - a name of a method being searched
/// \return a pointer to a pricing if found, otherwise 0
inline PricingMethod * getMethod(Mode mode, _bstr_t const & modelName, _bstr_t const & optionName, _bstr_t const & methodName)
{
    Model * model = getModel(mode, modelName);
    Option *option = getOption(mode, optionName);

    return getMethod(mode, model, option, methodName);
}

/// /return number of elements in an array of VARs
inline int getParArraySize(VAR * array)
{
    int i = 0;

    for (; array->Vtype != PREMIA_NULLTYPE; ++array)
        ++i;

    return i;
}

/// Retrieves a vector of entity parameters
/// \param T - entity type
/// \param e - an entity
/// \return a vector of entity parameters

    inline void getParameters(VAR * v, int stopper, std::vector<VAR> & result)
    {
        for (; stopper-- && v->Vtype != PREMIA_NULLTYPE; ++v)
        {
            switch (v->Vtype)
            {
            case NUMFUNC_1: 
                getParameters(v->Val.V_NUMFUNC_1->Par, MAX_PAR, result);
                break;
            case NUMFUNC_2:
                getParameters(v->Val.V_NUMFUNC_2->Par, MAX_PAR, result);
                break;
            case NUMFUNC_ND:
                getParameters(v->Val.V_NUMFUNC_ND->Par, MAX_PAR, result);
                break;
            default:
                if (v->Vsetable == SETABLE)
                    result.push_back(*v);
                break;
            }
        }
    }

    inline std::vector<VAR> getParameters(Model * model)
    {
        std::vector<VAR> v;
        getParameters(reinterpret_cast<VAR*>(model->TypeModel), model->nvar, v);
        return v;
    }

    inline std::vector<VAR> getParameters(Option * option)
    {
        std::vector<VAR> v;
        getParameters(reinterpret_cast<VAR*>(option->TypeOpt), option->nvar_setable, v);
        return v;
    }

    /// Retrieves a vector of a pricing method parameters
    /// \param e - a pricing method
    /// \return a vector of a pricing method parameters
    inline std::vector<VAR>    getParameters(PricingMethod * e)
    {
        return std::vector<VAR>(e->Par, e->Par + getParArraySize(e->Par));
    }

    inline bool isEnumWithParameters(VAR const & V)
    {
        if (V.Vtype == ENUM)
        {
            PremiaEnum * e = V.Val.V_ENUM.members;

            for (PremiaEnumMember * em = e->members; em->label; 
                em = reinterpret_cast<PremiaEnumMember*>(reinterpret_cast<char*>(em) + e->size))
            {
                if (em->nvar > 0)
                    return true;
            }
        }

        return false;
    }



/// Tag class used to differentiate between input and output parameters of a pricing method
/// \param T == PricingMethod
template <class T> struct ResultOf
{
    typedef T   pricing_method_t;

    /// wraps a pricing method
    ResultOf(pricing_method_t * ptr) : ptr_(ptr) {}

    /// \return the underlying pricing method
    pricing_method_t * get() const { return ptr_; }

private:
    pricing_method_t    * ptr_;
};

/// Retrieves a vector of a pricing method output parameters
/// \param e - a pricing method result wrapper
/// \return a vector of a pricing method output parameters
inline std::vector<VAR> getParameters(ResultOf<PricingMethod> * e)
{
    return std::vector<VAR>(e->get()->Res, e->get()->Res + getParArraySize(e->get()->Res));
}

inline bool init(Model * mod)
{
    mod->Init(mod);
    return true;
}

inline bool init(Option * opt, Model * mod)
{
    (opt->Init)(opt, mod);
    return true;
}

inline bool init(PricingMethod * met, Model * mod)
{
    met->Init(met, reinterpret_cast<Option*>(mod));
    return true;
}


/// Tests parameter setability 
/// \param v -- a parameter to test
/// \return true iff the parameter is setable
inline bool notsetable(VAR const & v)
{
    return strncmp(v.Vname, "-->", 3) == 0;
}

/// removes calculable fields from a parameter array
/// \param pars - a parameter array to be filtered
inline void removeCalculableFields(std::vector<VAR> & pars)
{
    pars.erase(std::remove_if(pars.begin(), pars.end(), notsetable), pars.end());
}

/// Predicate comparing VARs by their names
struct EqualNames
{
    bool operator ()(VAR const& a, VAR const& b) const
    {
        return strcmp(a.Vname, b.Vname) > 0;
    }
};

/// Accumulates parameters for pricing methods: PricingMethods with same name --> set of their parameters
struct MethodParameters : std::map<std::string, std::set<VAR, EqualNames> >   
{
    /// Constructor which initializes the map
    MethodParameters()
    {
        // for each pricing 
        for (Pricing ** pr_it = pricings_e; *pr_it; ++pr_it)
        {
            if (*pr_it == &BSND_STDND_pricing)
                continue;

            // iterate its methods
            for (PricingMethod ** m_it = (*pr_it)->Methods; *m_it; ++m_it)
            {
                PricingMethod * met = *m_it;

                (met->Init)(met, 0);

                // and for each method add all parameters to the map
				for (VAR * v_it = met->Par; v_it->Vtype != PREMIA_NULLTYPE; ++v_it)
                {
                    (*this)[met->Name].insert(*v_it);
                }
            }
        }
    }    
};

/// Singleton for the map: PricingMethods with same name --> set of their parameters
inline MethodParameters & methodParameters()
{
    static MethodParameters m;
    return m;
}

/// Tests whether a string is a number.
/// \param t - a string to be tested
/// \return true iff the string is a number
inline bool is_number(_bstr_t const & t)
{
    try { (double)_variant_t(t); }      // very stupid....
    catch (...)
    {
        return false;
    }
    return true;
}

/// Tests that a parameter value is a number.
/// If not shows a message box with an error.
/// \param key - a parameter key
/// \param value - a parameter value
/// \return true iff the value is a number
inline bool CheckNumericParameter(_bstr_t const & key, _bstr_t const & value)
{
    if (!is_number(value))
    {
        ::MessageBox(0, L"Parameter '" + key + "' should contain a number but contains " + value, L"Premia", MB_OK);
        return false;
    }

    return true;
}


/// Looks up a parameter by given key in an array of VARs.
/// Searching stops either when the parameter is found,
/// or end of array mark is reached, or stopper VARs are examined.
/// \param vars -- array of VARS where to perform the search
/// \param stopper -- maximum number of elements in vars to be examined
/// \param key -- a name of a parameter being searched
/// \return a pointer to parameter if found, null pointer otherwise
inline VAR * findParameter(VAR * vars, int stopper, std::string const & key)
{
    while (stopper-- && vars->Vtype != PREMIA_NULLTYPE)
    {
        if (vars->Vtype == NUMFUNC_1)
            if (VAR * ret = findParameter(vars->Val.V_NUMFUNC_1->Par, MAX_PAR, key))
                return ret;

        if (vars->Vtype == NUMFUNC_2)
            if (VAR * ret = findParameter(vars->Val.V_NUMFUNC_2->Par, MAX_PAR, key))
                return ret;

        if (key == vars->Vname)
            return vars;

        vars++;
    }

    return 0;
}

/// Looks up a parameter in a model's parameters.
/// \param m - a model where to look for a parameter
/// \param key - a name of a parameter being searched
/// \return a pointer to parameter if found, null pointer otherwise
inline VAR * findParameter(Model * m, const char * key)
{
    return findParameter(reinterpret_cast<VAR*>(m->TypeModel), m->nvar, key);
}

/// Looks up a parameter in an option's parameters.
/// \param o - an option where to look for a parameter
/// \param key - a name of a parameter being searched
/// \return a pointer to parameter if found, null pointer otherwise
inline VAR * findParameter(Option * o, const char * key)
{
    return findParameter(reinterpret_cast<VAR*>(o->TypeOpt), o->nvar_setable, key);
}

/// Looks up a parameter in a pricing method's parameters.
/// \param m - a pricing method where to look for a parameter
/// \param key - a name of a parameter being searched
/// \return a pointer to parameter if found, null pointer otherwise
inline VAR * findParameter(PricingMethod * m, const char * key)
{
    return findParameter(m->Par, MAX_PAR, key);
}

/// Loads a value to a given parameter.
/// \param V - VAR to be loaded
/// \param value - a string in Excel number format (uses commas as separator)
/// \return true iff the value is loaded successfully
inline bool loadParameter(VAR & V, const char * value)
{
    std::string ss = value;
    std::replace(ss.begin(), ss.end(), ',', '.');

    std::istringstream   s(ss);

	if (V.Vtype < FIRSTLEVEL && V.Vtype != ENUM && V.Vtype != FILENAME)
	{ 
		switch (true_typeV[V.Vtype])
		{
		case INT:
			s >> V.Val.V_INT; return s != 0;
		case 3 /*LONG*/:
			s >> V.Val.V_LONG; return s != 0;
		case DOUBLE:
			s >> V.Val.V_DOUBLE; return s != 0;
		default:
            return false;
		}
	}
    switch (V.Vtype)
    {
    case ENUM:
        s >> V.Val.V_ENUM.value; return s != 0;
    }

    return false;
}

/// Reads a parameter to an entity.
/// \param T -- entity type
/// \param target - an entity
/// \param key - a parameter key
/// \param value - a parameter value
template <class T>
bool putParameter(T * target, _bstr_t key, _bstr_t value)
{
    if (VAR * v = findParameter(target, key))
    {
        if (loadParameter(*v, value)) 
            return true;
        else
            ::MessageBox(0, L"Unable to load value = {" + value + L"} to parameter '" + key + L"'", L"Premia 10", MB_OK);
    }

    return false;
}


}}