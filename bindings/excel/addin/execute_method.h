#pragma once

namespace prxl {

    namespace lib {
/// Launches premia's mechanism for checking an entity's parameters
/// \param Target - an entity type
/// \param target - an entity 
/// \return true iff there are no errors
template <class Target>
    bool are_there_any_errors_in_params(Target * target)
{
    g_dup_printf = true;
    fopen_s(&g_dup_file, "checking_premia_input", "w");

    Planning plan;
    plan.Action = 'p';
    plan.VarNumber = 0;

    int result = target->Check(TOSCREEN, &plan, target);

    fclose(g_dup_file);
    g_dup_printf = false;

    if (result != OK)
    {
        std::ifstream infile("checking_premia_input");
        char buf[256];

        std::string err_msg;

        while (infile)
        {
            infile.getline(buf, 256);  
            std::string s = buf;

            if (s == "" || s == "All Right (ok: Return, no: n, h for Help) ? 	");
            else
                if (s == "Bad value:")
                {
                    infile.getline(buf, 256);
                    err_msg += buf;

                    infile.getline(buf, 256);
                    err_msg += ". ";
                    err_msg += buf;
                    err_msg += "\n";
                }
                else
                {
                    err_msg += buf;
                    err_msg += "\n";
                }
        }


        if (err_msg.size() > 0)
        {
            ::MessageBoxA(0, err_msg.c_str(), "Premia input data error", MB_OK);
            return false;
        }
    }

    return  true;
}


/// Checks whether there are any incompatibilities in parameter values of a model and an option
/// \param option - an option
/// \param model - a model
/// \param pricing - a pricing for the option and the model
/// \return true iff there are no errors
inline bool checkMixing(Option * option, Model * model, Pricing * pricing)
{
    // todo: check that pricing corresponds to the model and option
    g_dup_printf = true;
    fopen_s(&g_dup_file, "checking_premia_input", "w");

    int res = (*pricing->CheckMixing)(option, model);

    fclose(g_dup_file);
    g_dup_printf = false;

    if (res != OK)
    {
        std::ifstream infile("checking_premia_input");

        std::string s;
        char buf[256];

        while (infile)
        {
            infile.getline(buf, 256);
            s += buf;
            s += "\n";
        }

        ::MessageBoxA(0, s.c_str(), "Premia input data validation", MB_OK);

        return false;
    }

    return true;
}

/// Launches premia pricing method
/// \param model - a model 
/// \param option - an option
/// \param method - a pricing method to run
/// \param pricing - a pricing which the method belongs to
/// \return true iff there are no errors
inline bool executeMethod(Model * model, Option * option, PricingMethod * method, Pricing * pricing)
{
    if (!checkMixing(option, model, pricing))
        return S_OK;

    int result = method->Compute(option->TypeOpt, model->TypeModel, method);

    if (result != OK)
    {
        InitErrorMsg();

        ::MessageBoxA(0, error_msg[result], "Premia input data validation", MB_OK);

    }

    return result == OK;
}
}
}