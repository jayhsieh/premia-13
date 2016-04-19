#pragma once

#include <boost/format.hpp>
#include <boost/foreach.hpp>
#include <premia/generator/api/var_wrapper.h>
#include <premia/generator/utils/formatter_dsl.h>

namespace premia {
namespace pygen {
namespace fsharp {

    using api::EnumValue;

	namespace print {

		/// \brief returns F# type for a VAR
		struct fs_type_getter : var_visitor<fs_type_getter, const char *>
		{
			const char * operator () (Numeric<int> const &x) const		{ return "int"; }
			const char * operator () (Numeric<long> const &x) const		{ return "int64"; }
			const char * operator () (Numeric<double> const &x) const	{ return "float"; }
			const char * operator () (std::string const &x) const		{ return "string"; }
			const char * operator () (EnumValue const &x) const		{ return x.type->label.c_str(); }
			const char * operator () (std::vector<double> const &) const{ return "seq<double>"; }
		};

		/// prints a member declaration corresponding to the given VAR
		inline Formatter& memberDecl(Formatter & out, NamedVar const & vr)
		{
			out ("FLD_NAME",  vr.name)
				("FLD_TYPE",  fs_type_getter().apply(vr.value)) << "%FLD_NAME% : %FLD_TYPE%";

			return out;
		}

		/// \brief returns a string representation of a VAR value
		struct var_value_printer : var_visitor<var_value_printer, std::string>
		{
			/// for ENUM
			std::string operator () (EnumValue const & x) const
			{
				return (boost::format("%1%.%2%") % x.type->label % x.type->members.find(x.value)->second.label).str();
			}

			/// for FILENAME
			std::string operator () (std::string const & x) const
			{
				return (boost::format("@\"%1%\"") % x).str();
			}

			/// for INT and DOUBLE based types
			template <class T>
			std::string operator () (Numeric<T> const & i) const
			{
				std::stringstream out;
				out.setf(std::ios::showpoint);
				out << i.value;
				return out.str();
			}

			/// for LONG based types
			std::string operator () (Numeric<long> const & i) const
			{
				return str(i.value) + "L";
			}

			/// for PNLVECT and PNLVECTCOMPACT
			std::string operator () (std::vector<double> const & i) const
			{
				std::stringstream out;
				out.setf(std::ios::showpoint);
				out << "[";
				BOOST_FOREACH(double x, i)
					out << x << ";";
				out << "]";
				return out.str();
			}		
		};

		/// prints initialization clause for the VAR
		inline Formatter& memberIni(Formatter & out, NamedVar const & vr)
		{
			out ("FLD_NAME", vr.name)
				("FLD_VALUE", var_value_printer().apply(vr.value))
				<< "%FLD_NAME% = %FLD_VALUE%";

			return out;
		}

		/// \brief returns code sending VAR content to Premia used at Premia data structures initialization
		struct param_writer : var_visitor<param_writer, const char *>
		{
			const char * operator () (Numeric<int> const &x) const		{ return "write_int(x.%PROP_NAME%)"; }
			const char * operator () (Numeric<long> const &x) const		{ return "write_long(x.%PROP_NAME%)"; }
			const char * operator () (Numeric<double> const &x) const	{ return "write_double(x.%PROP_NAME%)"; }
			const char * operator () (std::string const &x) const		{ return "write_filename(x.%PROP_NAME%)"; }
			const char * operator () (EnumValue const &x) const		{ return "write_enum(x.%PROP_NAME% |> int)"; }
			const char * operator () (std::vector<double> const&) const { return "writearray(x.%PROP_NAME%)"; }
		};

		/// prints code sending VAR content to Premia used at Premia data structures initialization
		inline void copy_param_fs(Formatter & out, VarList::const_reference vr)
		{
			out("PROP_NAME", vr.name) << param_writer().apply(vr.value);
		}
	}


}}}
