#pragma once

#include <iostream>
#include <premia/generator/api/model.h>
#include <premia/generator/api/option.h>

/// \brief extract from a pricing name a model and a family names the pricing corresponds to
/// \param pricing_name a pricing name to parse
/// \param model_name [out] extracted model name
/// \param family_name [out] extracted family name
inline void parsePricingName(std::string const & pricing_name, std::string & model_name, std::string & family_name)
{
	// since some model IDs contain '_', we should look for the last one
	std::string::size_type pos = pricing_name.find_last_of('_');

	if (pos == std::string::npos)
		throw std::string("there is no '_' in pricing name: ") + pricing_name;

	model_name = pricing_name.substr(0, pos);
	family_name = pricing_name.substr(pos + 1);
}

namespace premia {
namespace pygen  {
namespace api    {

	struct PricingMethod;

	/// \brief A wrapper over a Premia pricing
	struct Pricing : boost::noncopyable
	{
		/// a pointer to a corresponding model wrapper
		Model const  * model;
		/// a pointer to a corresponding family wrapper
		Family		 * family;
		/// a pointer to the native pricing
		::Pricing **   source;
		/// a reference to the parent asset
		Asset const&   asset;

		/// a list of all methods
		boost::ptr_list<PricingMethod>	methods;

		std::map<Option const*, std::list<PricingMethod*> >	methods_for_options;

		/// collects all methods in the pricing
		/// \brief pricing a pointer to the native pricing
		/// \param asset a reference to the parent asset
		template <class Ctx>
   		     Pricing(Ctx & ctx, ::Pricing** pricing, Asset const & asset);
	};

	/// \brief gets the ordinal number of the pricing
	inline int id(Pricing const & pricing)
	{
		return pricing.source - pricing.asset.source->pricings;
	}

    struct empty_method_exception {};

	/// \brief A wrapper over a Premia pricing method
	struct PricingMethod : boost::noncopyable
	{
		/// \brief initializes the method
		/// \param source a pointer to the native pricing method
		/// \pricing a reference to the parent pricing
		template <class Ctx>
		PricingMethod(Ctx & ctx, ::PricingMethod **source, Pricing & pricing)
			:	source(source), pricing(pricing), name((*source)->Name)
		{
            bool initialized = false;

			// obtaining a list of options compatible with the method
			BOOST_FOREACH(Option & opt, pricing.family->options)
			{
				if ((*(*source)->CheckOpt)(*opt.source, *pricing.model->source) == OK)
				{
                    if (!initialized)
                    {
                        // initializing the native method
                        (*(*source)->Init)(*source, *opt.source);
                        initialized = true;
                    }

					compatible_options.push_back(&opt);

					opt.methods_for_models[pricing.model].insert(this);
				}
			}

            if (!initialized)
            {
                // it is very bad
                (*(*source)->Init)(*source, *pricing.family->options.front().source);

                std::cout << "Pricing method " << name << " has no compatible options" << std::endl;

                throw empty_method_exception();
            }

            // getting variables
            getVars(ctx, reinterpret_cast<VAR*>((*source)->Par), MAX_PAR, members);

            // getting method results
			getVars(ctx, reinterpret_cast<VAR*>((*source)->Res), MAX_PAR, results);
		}

		/// a pointer to the native pricing method
		::PricingMethod **  source;
		/// a reference to the parent pricing
		Pricing const &		pricing;
		/// corrected name of the method
		std::string const   name;
		/// linearized list of the method members
		VarList				members;
		/// linearized list of the method result parameters
		VarList				results;
		/// list of the options compatible with the method
		std::list<Option const*>	compatible_options;
	};

    inline const char * correctedFilename(PricingMethod const & m)
    {
        if (const char * corrected = (*m.source)->HelpFilenameHint)
        {
            return corrected;
        }
        return (*m.source)->Name;
    }

    inline fs::path relativeDocPath(PricingMethod const &method)
    {
        std::string model_name = correctedFilename(*method.pricing.model);
        std::string method_name = correctedFilename(method);
        std::string family_name = method.pricing.family->name;
        boost::to_lower(model_name);
        boost::to_lower(method_name);
        boost::to_lower(family_name); 
        std::string pricing_name = model_name + "_" + family_name;

        return fs::path("mod") / model_name / pricing_name / method_name;
    }

    inline fs::path relativeHtmlPath(PricingMethod const &method)
    {
        std::string model_name = correctedFilename(*method.pricing.model);
        std::string method_name = correctedFilename(method);
        std::string family_name = method.pricing.family->name;
        boost::to_lower(model_name);
        boost::to_lower(method_name);
        boost::to_lower(family_name); 
        std::string pricing_name = model_name + "_" + family_name;

        return fs::path("mod") / model_name / pricing_name / (method_name + "_doc") / (method_name + "_doc.html");
    }


	/// \brief gets the ordinal number of the method in the pricing
	inline int id(PricingMethod const & m)
	{
		return m.source - (*m.pricing.source)->Methods;
	}

    struct empty_pricing_exception {};

	/// \brief holds all the Premia pricings 
	struct Pricings
	{
		/// \brief collects the pricings
		template <class Ctx>
		Pricings(Ctx & ctx) 
		{
			std::set< ::Pricing*> pricings_processed;

			for (PremiaAsset * asset = premia_assets; asset->name; ++asset)
			{
				for (::Pricing ** pricing = asset->pricings; *pricing; ++pricing)
				{
					if (pricings_processed.find(*pricing) == pricings_processed.end())
					{
                        try {
    						Pricing * p = new Pricing(ctx, pricing, ctx.assets.lookup(asset));
    						pricings.push_back(p);
                        } catch (empty_pricing_exception) {
                            std::cout << "Pricing " << (*pricing)->ID << " is empty" << std::endl;
                        }
						pricings_processed.insert(*pricing);
					}
				}
			}
		}

		/// list of the pricings
		boost::ptr_list<Pricing>	pricings;
	};

		template <class Ctx>
		Pricing::Pricing(Ctx & ctx, ::Pricing** pricing, Asset const & asset)
			:	source(pricing), asset(asset)
		{
			std::string model_name, family_name;

			parsePricingName((*pricing)->ID, model_name, family_name);

			// corresponding model
			Family * family = ctx.families.lookup(family_name);
			// corresponding family
			Model * model = ctx.models.lookup(model_name);

			this->family = family;
			this->model = model;

			for (::PricingMethod **method = (*source)->Methods; *method; ++method)
			{
                try {
    				PricingMethod* pm = new PricingMethod(ctx, method, *this);
    				methods.push_back(pm);
    				BOOST_FOREACH(Option const * opt, pm->compatible_options)
    				{
    					methods_for_options[opt].push_back(pm);
    				}
                } catch (empty_method_exception) {}
			}

            if (methods_for_options.empty())
                throw empty_pricing_exception();

            // register the pricing in its model
            model->pricings.push_back(this);
		}

}

using api::Pricing;
using api::PricingMethod;
using api::Pricings;
using api::id;

}}
