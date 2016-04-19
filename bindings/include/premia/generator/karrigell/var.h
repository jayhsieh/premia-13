#pragma once

#include <premia/generator/utils/formatter_dsl.h>
#include <premia/generator/utils/symbols.h>
#include <premia/generator/karrigell/ctx.h>
#include <premia/generator/api/var_wrapper.h>

namespace premia	{
namespace pygen		{
namespace karrigell {

    using api::EnumValue;

	namespace print 
	{
        inline void FieldEx(Formatter & out, NamedVar const & vr);
	    
        inline void enumFields(Formatter &out, NamedVar const &vr)
        {
            FieldEx(out("OBJ", "pmem.get_%LABEL%()")("PREFIX", "self.fullName+'_get_%LABEL%_'"), vr);
        }

	    inline void enumInits(Formatter &out, api::Enum const &e)
	    {
            int idx = 0;
            out << "def initEnum(self, x,e):";
            BOOST_FOREACH(api::Enum::Members::const_reference em, e.members)
            {
                out("IDX",idx)("LABEL", em.second.label) 
                    << "   if e == %IDX%: x.set_%LABEL%()";
                ++idx;
            }
            out << "";
	    }

        inline void enumChoicesEx(Formatter &out, api::Enum const &e)
        {
            out << "def getFields(self, e):";
            int idx = 0;
            BOOST_FOREACH(api::Enum::Members::const_reference em, e.members)
            {
                out("IDX",idx)("LABEL", em.second.label) 
                    << +(seq, 
                    "if e == %IDX%: return [", +(seq, 
                        foreach_x(em.second.params, enumFields)),
                    "]","");
                ++idx;
            }
        }

		struct include_enums : var_visitor<include_enums>
		{
			include_enums(Formatter & out) : out(out){}

         template <class T>
			void operator () (Numeric<T> const & i) 
			{
			}

			void operator () (std::string const & i)  
			{
			}

			void operator () (std::vector<double> const & i) 
			{
			}
			
			void operator() (EnumValue const & e) 
			{
				out("ENUM_NAME", e.type->label) 
				   << "from kspremia.enum.%ENUM_NAME% import %ENUM_NAME%";
			}
			
			
		private:
			Formatter & out;
		};

      inline void includeEnums(Formatter &out, NamedVar const &vr)
		{
		   include_enums(out).apply(vr.value);
		}

      inline void commonHeader(Formatter &out)
      {
         out << (seq,
            "from kspremia.scalar import Scalar",
            "from kspremia.vector import Vector",
            "from kspremia.vector_compact import VectorCompact"
            );            
      }

		template <class Scalar>
			std::string const & tostr(api::Range<Scalar> const * c)
		{
			// since number of different ranges is small, let's memoize their representations
			typedef std::map<Range<Scalar> const *, std::string> Cache;

			static Cache	cache;

			// looking in the cache
			typename Cache::const_iterator it = cache.find(c);

			if (it != cache.end())
				return it->second;
			else	// if not in the cache
			{	
				std::string res;

				if (c != 0) 
				{
					res = "&#8745;";

					res += c->low_inclusive ? "[" : "(";

					res += c->has_low() ? str(c->low) : "-&infin;";

					res += ";";

					res += c->has_hi() ? str(c->hi) : "+&infin;";

					res += c->hi_inclusive ? "]" : ")";
				}

				cache[c] = res;

				return cache[c];
			}
		}

		template <class T> const char * symbol() { return "&#8484;"; }
		template <> const char * symbol<double>(){ return "&#8477;"; }

		template <class T> const char * converter() { return "int"; }
		template <> const char * converter<double>(){ return "float"; }
		template <> const char * converter<long>  (){ return "long"; }

		struct field_printer : var_visitor<field_printer>
		{
			field_printer(Formatter & out, VAR const *src, bool iterable) : out(out), src(src), iterable(iterable) {}
			
			template <class Scalar>
				void operator () (Numeric<Scalar> const & i) 
			{
				out("CONSTR",tostr(i.constraint))
				   ("SYMB", symbol<Scalar>())
				   ("CONV", converter<Scalar>())
				   ("ITERABLE", iterable ? "True" : "False")	
				   << "Scalar('%VAR_NAME%', '%FRIENDLY_NAME%', %PREFIX%+'_%VAR_NAME%', '%SYMB%%CONSTR%', %ONCHANGE%, %ITERABLE%, %CONV%),"
;
			}

			void operator () (std::string const & i)  
			{
			   out << "Scalar('%VAR_NAME%', '%FRIENDLY_NAME%', %PREFIX%+'_%VAR_NAME%', '', '', False, str),";
			}

			void operator () (std::vector<double> const & i) 
			{
			   if (src && src->Vtype==PNLVECTCOMPACT)
			      out 
			      << "VectorCompact('%VAR_NAME%', '%FRIENDLY_NAME%', %PREFIX%+'_%VAR_NAME%'),";
			   else
			      out 
			      << "Vector('%VAR_NAME%', '%FRIENDLY_NAME%', %PREFIX%+'_%VAR_NAME%'),";
			}
			
			void operator() (EnumValue const & e) 
			{
			    out("ENUM_TYPE", e.type->label) 
				   << "%ENUM_TYPE%('%VAR_NAME%', '%FRIENDLY_NAME%', %PREFIX%+'_%VAR_NAME%'),";
			}
		private:
			Formatter & out;
         VAR const * src;
         bool iterable;
		};
		
		inline void FieldVal(Formatter & out, NamedVar const & vr)
		{
		   if (vr.has_setter) std::cout << vr.name << std::endl;
		   std::string submit = vr.has_setter ? ", onchange='submit();'": "";
		   std::string onchange = vr.has_setter ? "'submit();'": "''";
			field_printer(out("SUBMIT", submit)("ONCHANGE",onchange), vr.src, vr.iterable).apply(vr.value);
		}

		inline void FieldEx(Formatter & out, NamedVar const & vr)
		{
		    std::string star;
			out("VAR_NAME", vr.name)
			   ("VAR_ID", symbol_utils::replaceNonAlnumCharacters(out.lookupVar("OBJ").c_str()) + "_" + vr.name)
               ("FRIENDLY_NAME", vr.src->Vname)
			   << (seq, call(boost::bind(FieldVal, _1, vr)));
		}

		inline void Field(Formatter & out, NamedVar const & vr)
		{
		   FieldEx(out("PREFIX", std::string("'")+symbol_utils::replaceNonAlnumCharacters(out.lookupVar("OBJ").c_str())+"'"), vr);
		}
	}


}}}
