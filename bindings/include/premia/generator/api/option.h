#pragma once

#include <string>
#include <map>

#include <boost/noncopyable.hpp>
#include <boost/ptr_container/ptr_list.hpp>

#include <premia/generator/api/var_conversion.h>

namespace premia {
namespace pygen  {
namespace api    {

	struct Family;
	struct PricingMethod;

	/// \brief A wrapper over a Premia option 
	struct Option : boost::noncopyable
	{
		/// \brief constructs a wrapper
		/// \param option a pointer to native option
		/// \param family a reference to family the option belongs to
		/// \param name corrected name of the option 
		/// \warning models should be initialized before 
		///	 since n-dimensional options take number of dimensions from corresponding models
		template <class Ctx>
		Option(Ctx & ctx, ::Option ** option, Family const & family, std::string const & name)
			:	family(family), name(name), source(option)
		{
			// if the option is n-dimensional, we need a model instance to correctly initialize it
			::Model * model = 
                strcmp((*option)->ID, "STDND") == 0  ? &BSND_model :
                strcmp((*option)->ID, "STDNDc") == 0 ? &COPULA_model :
                0;

			(*(*option)->Init)(*option, model);

			// getting the option variables
			getVars(ctx, reinterpret_cast<VAR*>((*option)->TypeOpt), (*option)->nvar, vars);
		}

		typedef std::map<Model const *, std::set<PricingMethod const *> >	AvailableMethods;

		/// family the option belongs to
		Family const & family;
		/// linearized list of the option variables
		VarList	  vars;
		/// corrected name of the option 
		std::string const name;
		/// pointer to the native option
		::Option ** source;

		AvailableMethods	methods_for_models;
	};

	/// \brief return the name of the native family (e.g. LIM, STD etc.)
	inline const char * getName(::Family * f)
	{
		return (**f)[0].ID;
	}

    inline const char * correctedFilename(Option const & opt)
    {
        if (const char * corrected = (*opt.source)->HelpFilenameHint)
        {
            return corrected;
        }
        return (*opt.source)->Name;
    }

	/// \brief A wrapper over an option family
	struct Family : boost::noncopyable
	{
		/// \brief collects all options for the family
		/// \param f a pointer to the native family
		/// \param asset an asset the family belongs to
		template <class Ctx>
		Family(Ctx& ctx, ::Family **f, Asset const& asset)
			: asset(asset), name(getName(*f)), source(f)
		{
			for (::Option ** option = **f; *option; ++option)
			{
				UsedIds used_names;

				options.push_back(new Option(ctx, option, *this, 
					correctIdentifierName((*option)->Name, used_names)));
			}
		}

		typedef boost::ptr_list<Option> Options;

		/// asset the family belongs to
		Asset const &			asset;
		/// 'id' of the family (STD, LIM etc.)
		std::string				name;
		/// list of the family options
		Options 				options;
		/// pointer to the native family
		::Family**				source;
	};

    inline fs::path relativeDocPath(Family const &f)
    {
        std::string familyName = f.name;
        boost::to_lower(familyName);
        return fs::path("opt") / familyName / familyName;
    }

    inline fs::path relativeHtmlPath(Family const &f)
    {
        std::string familyName = f.name;
        boost::to_lower(familyName);
        return fs::path("opt") / familyName / (familyName + "_doc") / (familyName + "_doc.html") ;
    }

    inline fs::path relativeDocPath(Option const &opt)
    {
        std::string optionName = correctedFilename(opt);
        boost::to_lower(optionName);
        return fs::path("opt") / boost::to_lower_copy(opt.family.name) / optionName;
    }

    inline fs::path relativeHtmlPath(Option const &opt)
    {
        std::string optionName = correctedFilename(opt);
        boost::to_lower(optionName);
        return fs::path("opt") / boost::to_lower_copy(opt.family.name) / (optionName + "_doc") / (optionName + "_doc.html");
    }

	/// \brief gets the ordinal number of the family
	inline int id(Family const & f)
	{
		return f.source - f.asset.source->families;
	}

	/// \brief gets the ordinal number of the option
	inline int id(Option const & opt)
	{
		return opt.source - **opt.family.source;
	}


	/// \brief represents all option family of Premia
	struct Families : boost::noncopyable
	{
		/// \brief collects all families
		template <class Ctx>
		Families(Ctx & ctx)
		{
			for (PremiaAsset * asset = premia_assets; asset->name; ++asset)
			{
				for (::Family ** f = asset->families; *f; ++f)
				{
					std::string name = getName(*f);

					if (!exist(name))
					{
						families.push_back(names[name] = new Family(ctx, f, ctx.assets.lookup(asset)));
					}
				}
			}
		}

		/// \brief checks whether a family with given name has been already processed
		bool exist(std::string const & family_name) const 
		{
			return names.find(family_name) != names.end();
		}
		
		typedef std::map<std::string, Family*> Names;

		/// \brief looks for a family with given name
		/// \returns pointer to the wrapper if found; throw an exception unless
		Family * lookup(std::string const & family_name) const 
		{
			Names::const_iterator it = names.find(family_name);

			if (it == names.end())
				throw std::string("cannot find a family with name ") + family_name;

			return it->second;
		}

		/// collected families
		boost::ptr_list<Family>			families;
		/// maps family names to pointers to families
		Names	names;
	};
}

using api::Option;
using api::Family;
using api::Families;
using api::id;

}}
