#pragma once

#include <list>
#include <string>
#include <boost/variant.hpp>
#include <premia/generator/api/var_constraint.h>

namespace premia {
namespace pygen {
namespace api   {

	/// \brief represents a value of a VAR of type int, long or double with a constraint if any
	/// \param Scalar int, long or double
	template <class Scalar>
		struct Numeric
		{
			/// default value of the VAR
			Scalar				  value;
			/// pointer to the constraint wrapper for the VAR if any, 0 otherwise
			Range<Scalar> const * constraint;

			/// \param value default value of the VAR
			/// \param pointer to the constraint wrapper for the VAR if any, 0 otherwise
			Numeric(Scalar value, Range<Scalar> const * constraint = 0)
				:	value(value), constraint(constraint)
			{}
		};

	/// \brief creates a wrapper for a VAR of numeric type
	/// \param value default value of the VAR
	/// \param Vtype VAR's type
	template <class Scalar>
		Numeric<Scalar>	numeric(Scalar value, int Vtype)
		{
			return Numeric<Scalar>(value, getRangeConstraint<Scalar>(Vtype));
		}

    struct EnumValue; 
    
	/// \brief this discriminated union represents all types that can be held in a VAR
	typedef 
		boost::variant<
			//bool,			// BOOL
			Numeric<int>,	// INT
			Numeric<long>,	// LONG
			Numeric<double>,	// DOUBLE
			std::string,		// FILENAME
			std::vector<double>,	// PNLVECT, PNLVECTCOMPACT
			boost::recursive_wrapper<EnumValue> // ENUM
		>	
		Var;
		
	/// \brief base class for visitors for Var
	/// \param Derived actual visitor class
	/// \param ResultT visitor return type
	template <class Derived, class ResultT = void>
		struct var_visitor : boost::static_visitor<ResultT>
	{
		/// applies visitor for the Var
		ResultT apply(Var const & v) const
		{
			return boost::apply_visitor(static_cast<Derived const &>(*this), v);
		}

		/// applies visitor for the Var
		ResultT apply(Var const & v)
		{
			return boost::apply_visitor(static_cast<Derived &>(*this), v);
		}
	};


	/// \brief VAR wrapper 
	/// abstract representation of a member of type VAR
	struct NamedVar
	{
		VAR const * src;	//!<  pointer to the source VAR 
		std::string name;	//!<  "corrected" name of the member
		Var			value;	//!<  wraps a default value, type and constraints for the field value
		bool        has_setter; //!< true iff the field has a non-trivial setter
		bool        iterable; //!< true iff the field allows iteration

		/// \param src pointer to the source VAR 
		/// \param name "corrected" name of the member
		/// \param value a wrapper holding the default value, type and constraints for the field value
		NamedVar(VAR const * src, std::string name, Var value, bool has_setter, bool iterable)
			: src(src), name(name), value(value), has_setter(has_setter), iterable(iterable)
		{}
	};

	/// list representing members of a model, an option or a pricing method
	typedef std::list<NamedVar>  VarList;
}

using api::Numeric;
using api::var_visitor;
using api::NamedVar;
using api::VarList;

}}
