#ifndef _premia_generator_api_enum_h_included_
#define _premia_generator_api_enum_h_included_

#include <map>
#include <string>
#include <boost/format.hpp>
#include <premia/exception.h>
#include <premia/generator/utils/symbols.h>
#include <premia/generator/api/enum.h>

namespace premia {
namespace pygen {
namespace api {

	//template <class Ctx>
    //	Var convertToVar(Ctx & ctx, VAR const & v);

	template <class Ctx>
    	int getVars(Ctx & ctx, VAR const* vars, int stopper, UsedIds & used_pars, VarList & var_list);

	/// \brief an enum instance: ordinal value + pointer to enumeration type
	struct EnumValue;
	
	inline std::string quote(const char * src)
	{
	    std::string s;
	    for (; *src; ++src)
	    {
	        if (*src == '\'')
	            s.push_back('\\');
	        s.push_back(*src);
	    }
	    return s;
	}

	/// \brief represents a premia enumeration type
	struct Enum 
	{
	    /// \brief represents a member of an enumeration
	    struct Member 
	    {
	         PremiaEnumMember const* src; 
	         std::string     label;  //!< corrected label of the member
	         std::string     quoted_original_label;
	         VarList         params; //!< optional parameters for the member
	    };
	    
		/// mapping from enumeration ordinal numbers to corresponding labels
		typedef std::map<int, Member> Members;

		/// constructs a wrapper over a premia enumeration
		/// \param em a pointer to the original enumeration
		/// \param correctedLabel unique name for the enumeration
		template <class Ctx>
		Enum(Ctx &ctx, PremiaEnum const * em, std::string const & correctedLabel)
			:	label(correctedLabel)
		{
			UsedIds used_members;

			/// iterating the enumeration members
			for (PremiaEnumMember const * eit = em->members; 
				eit->label != NULL; 
				eit = reinterpret_cast<PremiaEnumMember const *>(reinterpret_cast<char const *>(eit) + em->size))
			{
			    if (members.find(eit->key) != members.end())
			        throw premia::api::exception((boost::format("Members with same key (%0%) in enumeration %1%: %2% and %3%") % eit->key % label % members[eit->key].src->label % eit->label).str());
			        
				Member& m = members[eit->key];
				m.src = eit;
				m.quoted_original_label = quote(m.src->label);
				m.label = correctIdentifierName(eit->label, used_members);
				
				UsedIds used_params;
				getVars(ctx, eit->Par, eit->nvar, used_params, m.params);
				//for (int i = 0; i != eit->nvar; ++i)
				//    m.params.push_back(convertToVar(ctx, eit->Par[i]));
			}
		}

		std::string	const			label;	//!< enum unique name
		Members                     members;
	};
	
	/// \brief an enum instance: ordinal value + pointer to enumeration type
	struct EnumValue 
	{
		int			 value;				//!< ordinal value of the enum
		Enum const * type;				//!< pointer to the enum type

		/// constructs enum value
		/// \param v ordinal value of the enum
		/// \param t pointer to the enum type
		EnumValue(int v, Enum const * t) : value(v), type(t) {}
	};
	

	/// \brief represents all Premia enumeration types
	struct Enums 
	{
		/// mapping from source enumerations to wrappers over them
		typedef std::map<PremiaEnum const *, Enum> EnumsToLabels;

		/// registers an enumeration type
		/// \param e pointer to the original enumeration type
		/// \return a reference to a wrapper over it
		template <class Ctx>
		Enum const& insert(Ctx &ctx, PremiaEnum const * e)
		{
			// check if e has been registered
			EnumsToLabels::const_iterator it = enums_.find(e);

			// if no,
			if (it == enums_.end())
			{
				// generate unique name for it
				std::string corrected = correctIdentifierName(e->label, used_);

				// put to the map
				return enums_.insert(std::make_pair(e, Enum(ctx, e, corrected))).first->second;
			}

			return it->second;
		}

		/// \return a reference to the mapping from original enumerations to wrappers over them
		EnumsToLabels const & enums() const { return enums_; }

	private:
		/// mapping from original enumerations to wrappers over them
		EnumsToLabels	enums_;
		/// set of used identifiers for enumerations
		UsedIds			used_;
	};
}

using api::Enum;
using api::Enums;
using api::id;

}}

#endif
