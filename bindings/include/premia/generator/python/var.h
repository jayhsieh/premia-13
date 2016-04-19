#pragma once

#include <premia/generator/api/var_wrapper.h>

namespace premia {
namespace pygen {
namespace python {

    using api::EnumValue;

	namespace print {

		/// creates a string with a python expression representing the range constraint
		/// \param c pointer to the constraint if exist, 0 otherwise
		/// \return pointer to string with text representation, 0 otherwise
		template <class Scalar>
			const char * tostr(Range<Scalar> const * c)
		{
			// since number of different ranges is small, let's memoize their representations
			typedef std::map<Range<Scalar> const *, std::string> Cache;

			static Cache	cache;

			if (c == 0) return 0;

			// looking in the cache
			typename Cache::const_iterator it = cache.find(c);

			if (it != cache.end())
				return it->second.c_str();
			else	// if not in the cache
			{			
				std::string lhs = "", rhs; 

				// format the lower bound of the range if exists
				if (c->has_low())
					lhs += (boost::format("%1% %2% x") % c->low % (c->low_inclusive ? "<=" : "<")).str();

				// format the upper bound of the range if exists
				if (c->has_hi())
					rhs = (boost::format("x %2% %1%") % c->hi % (c->hi_inclusive ? "<=" : "<")).str();

				std::string res = c->has_low() && c->has_hi() ? lhs + " and " + rhs : lhs + rhs;

				cache[c] = res;

				return cache[c].c_str();
			}
		}

		template <class Scalar> const char* symbol() { return "Z"; }
		template <> const char* symbol<double>() { return "R"; }

		template <class Scalar>
			const char * tometa(Range<Scalar> const * c)
		{
			// since number of different ranges is small, let's memoize their representations
			typedef std::map<Range<Scalar> const *, std::string> Cache;

			static Cache	cache;

			// looking in the cache
			typename Cache::const_iterator it = cache.find(c);

			if (it != cache.end())
				return it->second.c_str();
			else	// if not in the cache
			{			
				std::string lhs, rhs; 

				// format the lower bound of the range if exists
				if (c && c->has_low())
					lhs = (boost::format("%1%, %2%") % (c->low_inclusive ? "'I'" : "'S'") % c->low).str();

				// format the upper bound of the range if exists
				if (c && c->has_hi())
					rhs = (boost::format("%1%, %2%") % (c->hi_inclusive ? "'I'" : "'S'") % c->hi).str();

				std::string res = (boost::format("['%1%', [%2%], [%3%]]") % symbol<Scalar>() % lhs % rhs).str();

				cache[c] = res;

				return cache[c].c_str();
			}
		}

		/// \brief prints Python assertion(s) for the given Premia VAR type
		struct property_value_assert : var_visitor<property_value_assert>
		{
			property_value_assert(Formatter & out, VAR const * src) : out(out), src(src) {}

			void scalar_iteration(const char * cast) {
				out("CAST", cast) 
				<< (seq, 
					"if type(x)==list and x[0]=='iterate':",
					"   def setter(v): self.%FLD_NAME%=v",
					"   iterables(r'%PRINTABLE%', setter, %CAST%(x[1]), %CAST%(x[2]), x[3])",
					"   x = x[1]");
			}


			/// assert for INT based types
			void operator () (Numeric<int> const & i) 
			{
				scalar_iteration("int");
				print_asserts(i, "assert(type(x) == int)");
			}

			/// assert for LONG based types
			void operator () (Numeric<long> const & i) 
			{
				scalar_iteration("long");
				print_asserts(i, "assert(type(x) == int or type(x) == long)");
			}

			/// assert for DOUBLE based types
			void operator () (Numeric<double> const & i) 
			{
				scalar_iteration("float");
				print_asserts(i, "assert(type(x) == int or type(x) == long or type(x) == float)");
			}

			/// assert for FILENAME
			void operator () (std::string const & i)  
			{
				out << (seq, 
					"if type(x)==unicode: x = str(x)",
					"assert(type(x) == str)");
			}

			/// assert for PNLVECT and PNLVECTCOMPACT
			void operator () (std::vector<double> const & i) 
			{
				if (src && src->Vtype==PNLVECTCOMPACT)
					out << (seq, 
							"x = list(x)", 
							"for i in range(len(x)):",
							"   if type(x[i])==list and x[i][0]=='iterate':",
							"      if len(x)==1:",
							"         def setter(v): self.%FLD_NAME% = len(self.%FLD_NAME%) * [v]",
							"      else:",
							"         def setter(v): self.%FLD_NAME%[i] = v",
							"      iterables(r'%PRINTABLE%', setter, float(x[i][1]), float(x[i][2]), x[i][3])",
							"      x[i] = x[i][1]",
							"if len(x)==1:",
							"   x = len(self.%FLD_NAME%) * x",
							"assert(len(x) == 0 or type(x[0]) == int or type(x[0]) == long or type(x[0]) == float)"); 
				else 
					out << (seq, 
							"x = list(x)", 
							"for i in range(len(x)):",
							"   if type(x[i])==list and x[i][0]=='iterate':",
							"      def setter(v): self.%FLD_NAME%[i] = v",
							"      iterables(r'%PRINTABLE%', setter, float(x[i][1]), float(x[i][2]), x[i][3])",
							"      x[i] = x[i][1]",
							"assert(len(x) == 0 or type(x[0]) == int or type(x[0]) == long or type(x[0]) == float)"); 
			}

			/// assert for ENUM
			void operator() (EnumValue const & e) 
			{
				out << (seq, 
					"if type(x) == list:",
					"  self.%FLD_NAME%.assign(x[0], x[1], iterables)",
					"  return",
					"assert(x.__class__ == self.%FLD_NAME%.__class__)");
			}
		private:

			/// prints assert for a numeric type and assert for the associated constraint of any
			template <class T> void print_asserts(Numeric<T> const & x, const char * pytype_constraint)
			{
				out << pytype_constraint;

				if (const char * constraint = tostr(x.constraint))
					out << str("assert ") + constraint + ", \"" + constraint + " where x is \" +str(x)";
			}
		private:
			Formatter & out;
			VAR const * src;
		};

		void copy_or_ignore_param(NamedVar const *to_copy, Formatter & out, NamedVar const & vr);
		void load_param(Formatter & out, NamedVar const & vr);

		/// prints body of the generated property setter for the variable
		inline void propertySetter(Formatter & out, NamedVar const & vr, VarList const & vars)
		{
			// prints type and constraint assertions
			property_value_assert(out, vr.src).apply(vr.value);
			
			if (vr.has_setter)
    			out << "self.makeCurrent()";
    			
			out << "self.%FLD_NAME% = x";

			if (vr.has_setter)
			{
			    out << (seq,
					"from premia import interop",
					"interop.setCurrentAsset(%ASSET_ID%)",
					"interop.setCurrentModel(%MODEL_ID%)",
					foreach_x(vars, boost::bind(print::copy_or_ignore_param, &vr, _1, _2)),
					"interop.stopWriteParameters()",
					"interop.readCurrentModel()",
					foreach_x(vars, print::load_param),
					"interop.stopReadParameters()")
					;
			}
		}

		/// python fields will be prefixed by underscore
		inline std::string fieldName(std::string const & s)
		{
			return "_" + s;
		}

		/// prints a property with a getter and a setter for a variable
		inline void PropertyEx(Formatter & out, NamedVar const& vr, VarList const & vars)
		{
			out("PROP_NAME", vr.name)
			   ("PRINTABLE", vr.src->Vname)
			   ("FLD_NAME", fieldName(vr.name))  << (seq, 
					"def get_%PROP_NAME%(self): return self.%FLD_NAME%", 
					"def set_%PROP_NAME%(self,x,iterables=None):", +
						call(boost::bind(print::propertySetter, _1, vr, boost::cref(vars))), 
					"%PROP_NAME% = property(get_%PROP_NAME%, set_%PROP_NAME%)", ""
			);
		}

		/// \brief prints a Python expression for the value of the given variable
		struct var_value_printer : var_visitor<var_value_printer, std::string>
		{
			/// for numeric types
			template <class T>	std::string operator () (Numeric<T> const & i) 
			{
				return str(i.value);
			}

			/// for FILENAME
			std::string operator () (std::string const & i)  
			{
				return (boost::format("r'%1%'") % i).str(); 
			}

			/// for PNLVECT and PNLVECTCOMPACT
			std::string operator () (std::vector<double> const & i) 
			{
				std::stringstream out;
				out << "[";
				BOOST_FOREACH(double x, i)
					 out << x << ",";
				out << "]";
				return out.str();
			}

			/// for ENUM
			std::string operator() (EnumValue const & e) 
			{
    			std::string enumname = e.type->label;
			    std::string label = e.type->members.find(e.value)->second.label;
			    return enumname + "(" + label + "())";
				//return (boost::format("%%ENUM_ID%%(%1%)") % e.value).str();
			}
		};

		/// prints an initializer by default value for a variable in the constructor of a Python class
		inline Formatter& initializer(Formatter & out, NamedVar const & vr)
		{
			// if the VAR is ENUM
			if (EnumValue const * e = boost::get<EnumValue>(&vr.value))
			{
                std::string label = e->type->members.find(e->value)->second.label;
                			
				// print import statements for the enum type
				out("PREMIA_LIB", premia_lib_name)
				   ("ENUM_ID", e->type->label)
				   ("LABEL", label) 
				   <<  (seq,
				        "from %PREMIA_LIB%.enum.%ENUM_ID% import %ENUM_ID%", 
				        "from %PREMIA_LIB%.enum.%ENUM_ID% import %LABEL%");
			}

			out("FLD_NAME",  fieldName(vr.name)) << "self.%FLD_NAME% = " + var_value_printer().apply(vr.value) ;

			return out;
		}

		/// prints constructor body of a Python class
		inline Formatter & Initializers (Formatter & out, VarList const & vars)
		{
			if (vars.empty())	out << "pass";
			else				out << foreach_x(vars, print::initializer);

			return out;
		}

		/// prints name of a property 
		inline void member(Formatter &out, NamedVar const & vr)
		{
			out("PROP_NAME", vr.name) << "'%PROP_NAME%',";
		}

		struct meta_writer : var_visitor<meta_writer>
		{
			meta_writer(Formatter & out, VAR const *src) : out(out), src(src) {}

			/// assert for INT based types
			template <class Scalar>
				void operator () (Numeric<Scalar> const & x) 
			{
				out("VAL", x.value)
					("CONSTRAINT", tometa(x.constraint)) 
					<< 
					"['%PROP_NAME%', 0, self.%FLD_NAME%, %CONSTRAINT%, %ITERABLE%, '%SETTER%'],";
			}

			/// assert for FILENAME
			void operator () (std::string const & i)  
			{
				out("VAL", i) << "['%PROP_NAME%', 3, self.%FLD_NAME%],";
			}

			/// assert for PNLVECT and PNLVECTCOMPACT
			void operator () (std::vector<double> const & i) 
			{
				std::stringstream acc;
				acc << "[";
				BOOST_FOREACH(double x, i)
					 acc << x << ",";
				acc << "]";

		      out("VAL", acc.str())
		      	("KEY", (src && src->Vtype==PNLVECTCOMPACT) ? 2 : 1) 
		      	<< 
			   		"['%PROP_NAME%', %KEY%, self.%FLD_NAME%],";
			}

			/// assert for ENUM
			void operator() (EnumValue const & e) 
			{
    			out("KEY", e.type->label)
    				("VAL", e.type->members.find(e.value)->second.quoted_original_label)
    				<< "['%PROP_NAME%', self.%FLD_NAME%.meta(), self.%FLD_NAME%._labels[self.%FLD_NAME%._value.key()], self.%FLD_NAME%._value.meta()],";
			}
		private:
			Formatter & out;
			VAR const * src;
		};


		/// prints meta information about the property 
		inline void meta(Formatter &out, NamedVar const & vr)
		{
			meta_writer(
				out("PROP_NAME", vr.src->Vname)
					("FLD_NAME", vr.name)
					("ITERABLE", vr.iterable)
					("SETTER", vr.has_setter ? vr.name : ""), 
				vr.src
				).apply(vr.value);
		}

		/// \brief returns a python statement copying a parameter from a Python class to Premia runtime
		struct param_writer : var_visitor<param_writer, const char *>
		{
			const char * operator () (Numeric<int> const &x) const		{ return "interop.write_int(self._%PROP_NAME%)"; }
			const char * operator () (Numeric<long> const &x) const		{ return "interop.write_long(self._%PROP_NAME%)"; }
			const char * operator () (Numeric<double> const &x) const	{ return "interop.write_double(self._%PROP_NAME%)"; }
			const char * operator () (std::string const &x) const		{ return "interop.write_filename(self._%PROP_NAME%)"; }
			const char * operator () (EnumValue const &x) const		{ return "self._%PROP_NAME%.export()"; }
			const char * operator () (std::vector<double> const&) const { return "writearray(self._%PROP_NAME%)"; }
		};

		/// prints a python statement copying a parameter from a Python class to Premia runtime
		inline void copy_param(Formatter & out, NamedVar const & vr)
		{
			out("PROP_NAME", vr.name) << param_writer().apply(vr.value);
		}

		struct param_ignorer : var_visitor<param_ignorer, const char *>
		{
			const char * operator () (Numeric<int> const &x) const		{ return "interop.ignore_int(self._%PROP_NAME%)"; }
			const char * operator () (Numeric<long> const &x) const		{ return "interop.ignore_long(self._%PROP_NAME%)"; }
			const char * operator () (Numeric<double> const &x) const	{ return "interop.ignore_double(self._%PROP_NAME%)"; }
			const char * operator () (std::string const &x) const		{ return "interop.ignore_filename(self._%PROP_NAME%)"; }
			const char * operator () (EnumValue const &x) const		{ return "self._%PROP_NAME%.export_ignoring()"; }
			const char * operator () (std::vector<double> const&) const { return "ignorearray(self._%PROP_NAME%)"; }
		};
		
		inline void copy_or_ignore_param(NamedVar const *to_copy, Formatter & out, NamedVar const & vr)
		{
		    if (vr.name == to_copy->name)
    			out("PROP_NAME", vr.name) << param_writer().apply(vr.value);
    		else
    			out("PROP_NAME", vr.name) << param_ignorer().apply(vr.value);
		}

		inline void ignore_param(Formatter & out, NamedVar const & vr)
		{
			out("PROP_NAME", vr.name) << param_ignorer().apply(vr.value);
		}
		
		struct param_loader : var_visitor<param_loader, const char *>
		{
			const char * operator () (Numeric<int> const &x) const		{ return "self._%PROP_NAME% = interop.read_int()"; }
			const char * operator () (Numeric<long> const &x) const		{ return "self._%PROP_NAME% = interop.read_long()"; }
			const char * operator () (Numeric<double> const &x) const	{ return "self._%PROP_NAME% = interop.read_double()"; }
			const char * operator () (std::string const &x) const		{ return "self._%PROP_NAME% = interop.read_filename()"; }
			const char * operator () (EnumValue const &x) const		{ return "self._%PROP_NAME%.load()"; }
			const char * operator () (std::vector<double> const&) const { return "self._%PROP_NAME% = readarray()"; }
		};

		inline void load_param(Formatter & out, NamedVar const & vr)
		{
			out("PROP_NAME", vr.name) << param_loader().apply(vr.value);
		}

		inline void assign_param(Formatter &out, NamedVar const &vr)
		{
			out("PROP_NAME", vr.name) << (seq, 
				"self.set_%PROP_NAME%(it.next(), iterables)");		
		}
	}

}}}
