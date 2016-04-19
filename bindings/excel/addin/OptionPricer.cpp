// OptionPricer.cpp : implémentation de COptionPricer

#include "stdafx.h"
#include "OptionPricer.h"

#include "storage.h"
#include "readingcells.h"

namespace prxl {

/// COptionPricer implementation members
struct COptionPricer::Impl
{
    /// problem mode
    lib::Mode    mode;
    /// problem model
    Model * model;
    /// problem option
    Option* option;
    /// problem method
    PricingMethod * method;
    /// problem iterable parameters
    fmt::iterators_t  iters_;
    /// a name of a requested field
    _bstr_t requested_field_;
};

/// Extracts a field from a result 
/// \param x - a result of a pricing method
/// \param name - a name of the field
/// \param pRes - output parameter where the extracted field should be put into
/// \return true iff a requested parameter is found
bool extractFieldFromResult(VAR * x, _bstr_t const & name, double *pRes)
{
    int count = lib::getParArraySize(x);

    for (VAR * v = x; count; count--, v++)
    {
        if (_stricmp(v->Vname, name) == 0)
        {
            *pRes = v->Val.V_DOUBLE;
            return true;
        }
    }

    ::MessageBox(0, "There is no such output parameter: '" + name + "'", L"Premia", MB_OK);

    return false;
}

/// Tests whether a string contains a formula 
/// \param s - a string to test
/// \return true iff the string is a formula
bool isFormula(BSTR s)
{
    return std::string(_bstr_t(s)).find('(') != std::string::npos;
}

/// Parses a formula strings for parameters of an entity and puts them to a premia structure.
/// \param T - entity type
/// \param target - an entity to be filled with parameters read
/// \param v - an array of entity parameters (for a pricing methods with same name their parameters are contained in a single array)
/// \param s - a formula string
/// \param q - a position in s where to start parsing
/// \return true iff there are no errors
template <class T>
bool parseAndLoadParameters(T * target, std::vector<VAR> const & v, 
                            std::string const & s, std::string::size_type q)
{
    // parameters read so far
    int no = 0;

    do 
    {
        // at begin q points to an opening parenthesis; after q points to a semicolon
        std::string::size_type p = q + 1;

        // looking for the end of the parameter
        q = s.find_first_of(";)", p);

        // if end of the string reached there is an error
        if (q == std::string::npos)
        {
            ::MessageBox(0, L"Formula should ends with ')'", L"Premia", MB_OK);

            return false;
        }

        // if the substring is not empty
        if (q > p)
        {
            // extract the parameter string
            std::string param = s.substr(p, q - p); 

            // check its validity
            lib::CheckNumericParameter(v.at(no).Vname, param.c_str());

            lib::putParameter(target, v.at(no).Vname, param.c_str());
        }

        ++no;

    } while(s[q] != ')');

    return true;
}

/// Extracts name of a function from a string.
/// \param s - a string which is a formula
/// \return name of a function 
std::string getFunctionName(std::string const & s)
{
    std::string::size_type q = s.find('(');

    return s.substr(0, q);
}

/// parses a formula that contains an option description
/// \param strOption - an option formula
/// \return a pointer to a created option, 0 if there are some errors
Option * parseOptionFormula(BSTR strOption)
{
    std::string s = _bstr_t(strOption);

    std::string option_name = getFunctionName(s);

    // position of the opening parenthesis
    std::string::size_type q = option_name.size();

    // look up the option
	if (Option * opt = lib::getOption(0, option_name.c_str()))
    {
        // extract parameters from the formula
        std::vector<VAR> v = lib::getParameters(opt);

        lib::removeCalculableFields(v);

        // load parameters to the option
        parseAndLoadParameters(opt, v, s, q);

        return opt;
    }
    else
    {
        ::MessageBox(0, "Cannot find an option named " + _bstr_t(option_name.c_str()), L"Premia", MB_OK);

        return 0;
    }

}

/// parses a formula that contains an model description
/// \param strModel - an option formula
/// \return a pointer to a created model, 0 if there are some errors
Model * parseModelFormula(BSTR strModel)
{
    std::string s = _bstr_t(strModel);

    std::string model_name = getFunctionName(s);

    // position of the opening parenthesis
    std::string::size_type q = model_name.size();

    // look up the model
	if (Model * mod = lib::getModel(0, model_name.c_str()))
    {
        // extract parameters from the formula
        std::vector<VAR> v = lib::getParameters(mod);

        lib::removeCalculableFields(v);

        // load parameters to the model
        parseAndLoadParameters(mod, v, s, q);

        return mod;
    }
    else
    {
        ::MessageBox(0, "Cannot find a model named " + _bstr_t(model_name.c_str()), L"Premia", MB_OK);

        return 0;
    }
}

/// parses a formula that contains a pricing method description
/// \param strMethod - a pricing method formula
/// \param mod - a pointer to a model for the pricing method
/// \param opt - a pointer to an option for the pricing method
/// \return a pointer to a created pricing method, 0 if there are some errors
PricingMethod * parseMethodFormula(BSTR strMethod, Model * mod, Option * opt)
{
    std::string s = _bstr_t(strMethod);

    std::string original_name = getFunctionName(s);

    // position of the opening parenthesis
    std::string::size_type q = original_name.size();

    // the method found but it is not compatible with mod and opt
    bool found_but_not_compatable;

    // lookup for the pricing method
	PricingMethod * met = lib::getMethod(0, mod, opt, original_name.c_str(), &found_but_not_compatable);

    // if not found
    if (met == 0)
    {
        bool b;
        // replace underscores in the method name by spaces
        std::replace(original_name.begin(), original_name.end(), '_', ' ');
        // try one more time
		met = lib::getMethod(0, mod, opt, original_name.c_str(), &b);
        found_but_not_compatable = found_but_not_compatable || b;
    }

    if (met)
    {
        // get an array of parameters of all methods with name == original_name
        std::set<VAR, lib::EqualNames> const & vs = lib::methodParameters()[original_name];

        std::vector<VAR> v(vs.begin(), vs.end());

        parseAndLoadParameters(met, v, s, q);

        return met;    
    }
    else
    {
        ::MessageBox(0, "Method " + _bstr_t(strMethod) + 
            (found_but_not_compatable 
            ? _bstr_t(" found but not compatible with model ") + mod->Name + " and option " + opt->Name 
            : _bstr_t(" not found in pricing ") + lib::getPricingName(mod,opt)), L"Premia", MB_OK);

        return 0;
    }
}

/// Loads a model from a string
/// \param model_id - a string that contains either a model instance label either a model formula
/// \param iters - output parameter for storing information about iterable parameter if any
/// \param pMode - output parameter for storing mode of the model
/// \return a pointer to a model created
Model * loadModel(BSTR model_id, fmt::iterators_t & iters, lib::Mode * pMode)
{
    if (isFormula(model_id))
    {
		*pMode = 0;

        return parseModelFormula(model_id);
    }
    else    // model_id is a model instance label
    {
        if (db::Entry<Model> const * mod = db::theModels().lookup(model_id))
        {
            *pMode = mod->getMode();

            return fmt::load(*mod, iters) ? mod->getData() : 0;
        }
        else
        {
            ::MessageBox(0, "There is no model instance with ID = " + _bstr_t(model_id), L"Premia", MB_OK);

            return 0;
        }
    }
}

/// Loads an option from a string
/// \param option_id - a string that contains either an option instance label either an option formula
/// \param iters - output parameter for storing information about iterable parameter if any
/// \return a pointer to an option created
Option * loadOption(BSTR option_id, fmt::iterators_t & iters)
{
    if (isFormula(option_id))
    {
        return parseOptionFormula(option_id);
    }
    else
    {
        if (db::Entry<Option> const * opt = db::theOptions().lookup(option_id))
        {
            return fmt::load(*opt, iters) ? opt->getData() : 0;
        }
        else
        {
            ::MessageBox(0, "There is no option instance with ID = " + _bstr_t(option_id), L"Premia", MB_OK);

            return 0;
        }
    }
}
/// Loads a pricing method from a string
/// \param method_id - a string that contains either a pricing method instance label either a pricing method formula
/// \param iters - output parameter for storing information about iterable parameter if any
/// \param mod - a pointer to a model for the method
/// \param opt - a pointer to an option for the method
/// \return a pointer to a pricing method created
PricingMethod * loadMethod(BSTR method_id, fmt::iterators_t & iters, Model * mod, Option * opt)
{
    if (isFormula(method_id))
    {
        return parseMethodFormula(method_id, mod, opt);
    }
    else
    {
        if (db::Entry<PricingMethod> const * met = db::theMethods().lookup(method_id))
        {
            return fmt::load(*met, iters) ? met->getData() : 0;
        }
        else
        {
            ::MessageBox(0, "There is no pricing method instance with ID = " + _bstr_t(method_id), L"Premia", MB_OK);

            return 0;
        }
    }
}


STDMETHODIMP COptionPricer::Init(BSTR requested_param, BSTR model_id, BSTR option_id, BSTR method_id)
{
    // reset pimpl
    impl_.reset(new Impl());

    impl_->requested_field_ = requested_param;

    if ((impl_->model = loadModel(model_id, impl_->iters_, &impl_->mode)) && 
        (impl_->option = loadOption(option_id, impl_->iters_)) &&
        (impl_->method = loadMethod(method_id, impl_->iters_, impl_->model, impl_->option)))
    {
        return S_OK;
    }

    impl_.reset(0);
    return E_FAIL;
}

/// Returns number of elements in an input values array
/// \param from - the first cell of the input array
/// \return number of elements in an input values array
int iteration_count(Excel::RangePtr from)
{
    int c = 0;

    for (; from->Value2 != _variant_t(); ++c, from = Down(from));

    return c;
}

STDMETHODIMP COptionPricer::GetExtents(long* pSize_1, long* pSize_2)
{
    if (impl_.get())
    {
        *pSize_1 = impl_->iters_.size() >= 1 ? iteration_count(fmt::derefCell(Right(impl_->iters_.front().getReferenceToKey()))) : 1;
        *pSize_2 = impl_->iters_.size() >= 2 ? iteration_count(fmt::derefCell(Right((*++impl_->iters_.begin()).getReferenceToKey()))) : 1;
        return S_OK;
    }
    
    return E_FAIL;
}

STDMETHODIMP COptionPricer::Run(long it_1, long it_2, double* pRes)
{
    if (impl_.get() == 0)
        return E_FAIL;

    switch (impl_->iters_.size())
    {
    case 2:
        {
            fmt::IterableParam &iter_2 = *++impl_->iters_.begin();
            Excel::RangePtr input_2 = fmt::derefCell(Right(iter_2.getReferenceToKey()))->GetOffset(it_2, 0);
            iter_2.update((_bstr_t)input_2->Value2);
        }
    case 1:
        {
            fmt::IterableParam &iter_1 = *impl_->iters_.begin();
            Excel::RangePtr input_1 = fmt::derefCell(Right(iter_1.getReferenceToKey()))->GetOffset(it_1, 0);
            iter_1.update((_bstr_t)input_1->Value2);
        }
    case 0:
        {
            if (lib::getPricing(impl_->mode, impl_->model, impl_->option) == 
                lib::getPricing(impl_->mode, impl_->method))
            {
                if (OK == (impl_->method->CheckOpt)(impl_->option, impl_->model))
                {
                    if (lib::executeMethod(impl_->model, impl_->option, impl_->method, lib::getPricing(impl_->mode, impl_->method)))
                    {
                        if (extractFieldFromResult(impl_->method->Res, impl_->requested_field_, pRes))
                            return S_OK;
                    }    
                }
                else
                {
                    ::MessageBox(0, _bstr_t("The method ") + impl_->method->Name + " is not compatible with"
                        " model " + impl_->model->Name + " and option " + impl_->option->Name, L"Premia", MB_OK);
                }
            }
            else
                ::MessageBox(0, _bstr_t("Pricing for model ") + impl_->model->Name + 
                " and option " + impl_->option->Name + " doesn't contain method " + impl_->method->Name, 
                L"Premia", MB_OK);
        }
        break;
    }           

    return S_OK;
}
}