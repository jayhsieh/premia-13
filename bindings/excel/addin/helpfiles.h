#pragma once

#include <boost/algorithm/string.hpp>
#include "registry.h"
#include "premia.h"

namespace prxl {
    
    /// Contains routines for accessing Premia help files
    namespace hlp {

    /// Corrects model ID for forming a correct path to the model helpfile
    /// \param mode - a model mode
    /// \param n - a model ID
    /// \return corrected ID
    inline const char * correctModelId(long mode, Model * model)
    {
		if (const char * corrected = model->HelpFilenameHint)
		{
			return corrected;
		}

        std::string name = model->ID;

        return 
//             name == "DUPIRE1D" ? "DUP1D" :
//            name == "MERTON1D" ? "MER1D" :
//            name == "HESTON1D" ? "HES1D" : 
//            name == "FOUQUEPAPANICOLAUSIRCAR1D" ? "FPS1D" :
//            name == "FOUQUEPAPANICOLAUSIRCAR2D" ? "FPS2D" :
//             name == "MERTONHESTON1D" ? "MERHES1D" : 
//             name == "CirPlus1D" ? "CIRPP1D" :
//             name == "SquaredGaussian1D" ? "SG1D" :
//             name == "LiRitchkenSankarasubramanian1D" ? "LRSHJM1D" : 
//             name == "HuntKennedyPelsser1D" ? "HK1D" : 
// 			name == "HULLWHITE1D" && mode == lib::equity ? "HW1D" :
            model->ID;
    }

    /// Corrects option ID for forming a correct path to the option helpfile
    /// \param mode - a option mode
    /// \param n - a option ID
    /// \return corrected ID
    inline const char * correctOptionName(Option * opt)
    {
		if (const char * corrected = opt->HelpFilenameHint)
		{
			return corrected;
		}

        std::string name = opt->Name;

        return 
//            name == "PayerSwaption" ? "PayerSwaptions" : 
//            name == "ReceiverBermudanSwaption" ? "ReceiverBermudanSwaptions" : 
//             name == "ReceiverSwaption" ? "ReceiverSwaptions" : 
//             name == "ZeroCouponBond" ? "zcbond" : 
//             name == "ZeroCouponCallBondAmer" ? "zccallbondamer" :
//            name == "ZeroCouponCallBondEuro" ? "zccallbondeuro" : 
//            name == "ZeroCouponPutBondAmer" ? "zcputbondamer" : 
//            name == "ZeroCouponPutBondEuro" ? "zcputbondeuro" : 
//             name == "CallDownOutDiscEuro" ? "CallDownOutEuro" :
            opt->Name;
    }

    /// Corrects pricing method ID for forming a correct path to the pricing method helpfile
    /// \param n - a pricing method ID
    /// \param modelName - a model name the method applies to 
    /// \return corrected ID
	inline const char * correctMethodName(PricingMethod * method, std::string corrected_model_name)
    {
		if (const char * corrected_name = method->HelpFilenameHint)
		{
			return corrected_name;
		}

        std::string name = method->Name;

        return 
//             name == "AP_Whaley" ? "AP_Waley" :
/*
                        name == "AP_CarrMadan" ? 
                        (corrected_model_name == "MER1D" ? "AP_Carr_Mer" :
                        corrected_model_name == "VARIANCEGAMMA1D" ? "AP_Carr_VG" : 
                        corrected_model_name == "TEMPEREDSTABLE1D" ? "AP_Carr_TS" :
                        corrected_model_name == "NIG1D" ? "AP_Carr_NIG" :
                        corrected_model_name == "KOU1D" ? "AP_Carr_Kou" :
                        method->Name) :
                    name == "MC_RobbinsMonro" ? 
                        (corrected_model_name == "HW1D" ? "MC_RobbinsMoro_HW" :
                        corrected_model_name == "HES1D" ? "MC_RobbinsMoro_Hes" :
                        method->Name) :
                    name == "MC_FixedAsian_RobbinsMonro" ? "MC_FixedAsian_RobbinsMoro" :
                        name == "FD_FixedAsian_LelievreDubois" ? "FD_FixedAsian_RodgerShi2" 
                        :*/
             method->Name;        
    }

	inline std::string slash(std::string s)
	{
		boost::algorithm::replace_all(s, "\\", "/");
		return s;
	}

    /// Returns a relative path to a model help file
    /// \param mode - a model mode
    /// \param model - a model
    /// \return a relative path to a model help file
    inline _bstr_t getHelpFileName(lib::Mode mode, Model * model, std::ostream & err_log)
    {
		(model->Init)(model);
        const char * modelName = correctModelId(mode, model);
        _bstr_t path = env::getManPath() + _bstr_t("\\mod\\") + modelName + "\\" + modelName + "_doc.pdf";

        if (!std::ifstream((const char *)path))
        {
            err_log 
				<<  "Unable to open help file '" << slash((const char*)path) << "' for model \n\t" 
                << "{ \n\t\tName = " << model->Name << ", \n\t\t"
                << "\n\t\tID = " << model->ID << ", "
                << "\n\t\tcorrected name = " << modelName 
                << "\n\t}\n" << std::endl;

            return "";
        }

        return path;
    }

    /// Returns a relative path to an option help file
    /// \param mode - an option mode
    /// \param model - an option
    /// \return a relative path to an option help file
    inline _bstr_t getHelpFileName(lib::Mode mode, Option * option, std::ostream & err_log)
    {
		Model * model = strcmp(option->ID, "STDND") == 0 || strcmp(option->ID, "STDNDc") == 0 ? &BSND_model : 0;

		(*option->Init)(option, model);

        _bstr_t path = env::getManPath() + _bstr_t("\\opt\\") + option->ID + "\\" + correctOptionName(option)+ "_doc.pdf";

        if (!std::ifstream((const char *)path))
        {
            err_log 
                <<  "Unable to open help file '" << slash((const char*)path) << "' for option "
                << "\n\t{ \n\t\tName = " << option->Name << ", "
                << "\n\t\tID = " << option->ID << ", "
                << "\n\t\tcorrected name = " << correctOptionName(option) 
                << "\n\t}\n" << std::endl;

            return "";
        }

        return path;
    }

    /// Returns a relative path to a pricing method help file
    /// \param mode - a pricing method mode
    /// \param model - a pricing method 
    /// \return a relative path to a pricing method file
    inline _bstr_t getHelpFileName(lib::Mode mode, PricingMethod * method, std::ostream & err_log)
    {
        Pricing * p = lib::getPricing(mode, method);

        std::string pr_name = p->ID;

        std::string::size_type pos = pr_name.rfind('_');

        std::string model_name = pr_name.substr(0, pos);
        std::string family_name = pr_name.substr(pos + 1);
        std::string corrected_family_name = /*family_name == "STDg" ? "STD" :*/ family_name;

		Model * model = lib::getModel(mode, model_name.c_str());

		std::string corrected_model_name = model ? correctModelId(mode, model) : model_name;

		Family *f = lib::getFamily(mode, _bstr_t(corrected_family_name.c_str()));

		(method->Init)(method, (*f)[0]);
        std::string corrected_method_name = correctMethodName(method, corrected_model_name);

        std::string path = env::getManPath() + "\\mod\\";
        path += corrected_model_name;
        path += "\\";
        path += corrected_model_name;
        path += "_";
        path += corrected_family_name;
        path += "\\";
        path += corrected_method_name;
        path += "_doc.pdf";

        if (!std::ifstream(path.c_str()))
        {

            err_log 
                <<  "Unable to open help file '" << slash(path) << "' for pricing method "
                << "\n\t{\n\t\tName = " << method->Name << ", "
                << "\n\t\tpricing_name = " << p->ID << ", "
                << "\n\t\tmodel_name = " << model_name << ", "
                << "\n\t\tfamily_name = " << family_name << ", "
                << "\n\t\tcorrected_model_name = " << corrected_model_name << ", "
                << "\n\t\tcorrected_family_name = " << corrected_family_name << ", "
                << "\n\t\tcorrected_method_name = " << corrected_method_name
                << "\n\t}\n" << std::endl;

            return "";
        }

        return path.c_str();
    }
 
} }
