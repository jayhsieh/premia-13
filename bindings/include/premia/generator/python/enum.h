#pragma once 

namespace premia {
namespace pygen  {
namespace python {

	namespace print {

		/// binds with the formatter label for the enumeration member and its ordinal number
		template <class Pair>
		inline Formatter & operator && (Formatter & out, Pair const & p)
		{
			return out("LABEL", p.second.label)("ID", p.first);
		}

		/// prints initializer part for python map that implements the enum type
		void Ini(Formatter & out, Enum::Members::const_reference p)
		{
			out("LABEL", p.second.quoted_original_label)
			   ("ID", p.first) 
			   << "%ID% : '%LABEL%',";
		}

		void Meta(Formatter & out, Enum::Members::const_reference p)
		{
			out("LABEL", p.second.quoted_original_label)
			   ("CLASS", p.second.label) 
			   << "('%LABEL%' , %CLASS%),";
		}

		void Assign(Formatter & out, Enum::Members::const_reference p)
		{
			out("LABEL", p.second.quoted_original_label)
			   ("CLASS", p.second.label) 
			   << "if label=='%LABEL%': self.set_%CLASS%()";
		}

		/// prints accessors to enum members
		void Prop(Formatter & out, Enum::Members::const_reference p)
		{
			(out && p) << (seq,
					//"@staticmethod", 
					//"def %LABEL%(): return %CLASS%(%ID%)", 
					"def set_%LABEL%(self): self._value = %LABEL%()", 
					"def get_%LABEL%(self):",
					"   assert type(self._value) == %LABEL%",
					"   return self._value", 
					"");
		}
		
		void Label(Formatter & out, Enum::Members::const_reference p)
		{
			(out && p) << "%LABEL%,";
		};
		
		void Choice(Formatter & out, Enum::Members::const_reference p)
		{
			(out && p)("MEMBERS_LEN", p.second.params.size()) << (seq,
			    "", 
			    "class %LABEL%(object):", +(seq,
			        "def __init__(self):", 
			        +call(boost::bind(print::Initializers, _1, p.second.params)), "",
						"def meta(self): ", +(seq, 
							"return [", +foreach_x(p.second.params, print::meta), "]"),
			        "def key(self): return %ID%", "",
			        foreach_x(p.second.params, boost::bind(print::PropertyEx, _1, _2, boost::cref(p.second.params))),
		 				"def assign(self, args, iterables):", +(seq,
		 					"assert(len(args) == %MEMBERS_LEN%)",
							"it = args.__iter__()",
							foreach_x(p.second.params, print::assign_param)
		 					), "",
			        "def export(self):", +(seq,
			            "from premia import interop",
			            "interop.write_enum(%ID%)",
			            foreach_x(p.second.params, print::copy_param)), "",
			        "def export_ignoring(self):", +(seq,
			            "from premia import interop",
			            "interop.ignore_enum(%ID%)",
			            foreach_x(p.second.params, print::ignore_param)), "",
			        "def __repr__(self): return '%LABEL%'"));
        }		
	}

	/// generates an enumeration wrapper 
	PyCtx const & operator << (PyCtx const & ctx, Enum const & e)  
	{
		Formatter f(ctx.filename(e));
		f("CLASS", e.label) << (seq,
		         "from premia.common import *", "",
				 "class %CLASS%:", +(seq,
					"_labels = {", +(
						foreach_x(e.members, print::Ini)),
					"}",
					"def assign(self, label, data, iterables):", +(seq, 
						foreach_x(e.members, print::Assign),
						"self._value.assign(data, iterables)"), "",
					"@staticmethod",
					"def meta():", +(seq, "choices = [", 
						+foreach_x(e.members, print::Meta),
					"]",
					"return [(k, v().meta()) for (k,v) in choices]"),
					"def __init__(self, newval): self._value = newval", "",
					"",
					"def __repr__(self): return self._value.__repr__()",
					"",
					"def export(self): self._value.export()",
					"",
					foreach_x(e.members, print::Prop)), "",
				//"class Base(object):", +(seq,
				//    "@staticmethod",
    			//	"def choices():", +(seq, 
    			//	    "return [", +foreach_x(e.members, print::Label), "]")),
				foreach_x(e.members, print::Choice)
			);
			
			
		/* ctx.out(2) << "   " << e.label << std::endl; */

		return ctx;
	}

	/// generates wrappers for all Premia enumerations
	PyCtx const & operator << (PyCtx const & ctx, Enums const & e) 
	{
		ctx.out(1) << "Generating enums:...";

		ctx.create(ctx.enumDir());
		
		BOOST_FOREACH(Enums::EnumsToLabels::const_reference p, e.enums())
		{
		    ctx << p.second;
		}
		
		//boost::for_each(e.enums() | map_values, ctx << lm::_1);

		ctx.out(1) << "ok!" << std::endl;
		return ctx;
	}

}}}
