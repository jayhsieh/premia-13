#pragma once

namespace premia {
namespace pygen {
namespace python {

	/// generates wrapper for an option
	PyCtx const & operator << (PyCtx const & ctx, Option const & opt)
	{
		Formatter f(ctx.filename(opt));
		f	("NAME",	  opt.name)
			("OPTION_ID", id(opt))
			("ASSET_ID",  id(opt.family.asset))
			("FAMILY_ID", id(opt.family))
			("FAMILY_NAME", opt.family.name)
			("MEMBERS_LEN", opt.vars.size())
			<< (seq, 
				"from ...common import *", "",
				"class %NAME%(object):", +(seq, "",
					"def __init__(self):", 
						+call(boost::bind(print::Initializers, _1, boost::ref(opt.vars))), "", 
					foreach_x(opt.vars, boost::bind(print::PropertyEx, _1, _2, boost::cref(opt.vars))), "", 
					"def __repr__(self): return getRepr(self, 'Option')", 
					"@staticmethod",
					"def familyName(): return '%FAMILY_NAME%'",
 					"def makeCurrent(self):", +(seq,
						"from premia import interop",
 						"interop.setCurrentAsset(%ASSET_ID%)",
 						"interop.setCurrentOption(%FAMILY_ID%, %OPTION_ID%)",
						foreach_x(opt.vars, print::copy_param),
						"interop.stopWriteParameters()"
 					), "",
					"@staticmethod",
	 				"def create(args, iterables):", +(seq,
	 					"self = %NAME%()",
 						"assert(len(args) == %MEMBERS_LEN%)",
						"it = args.__iter__()",
						foreach_x(opt.vars, print::assign_param),
						"return self"
 						), "",
					"@staticmethod",
					"def parameters(): ", +(seq, 
						"return [", +foreach_x(opt.vars, print::member), "]"),
					"def meta(self): ", +(seq, 
						"return [", +foreach_x(opt.vars, print::meta), "]")
					))
			;

		return ctx;
	}

	namespace print {

		/// prints import option clause in a family file
		void OptionUsing(Formatter & out, Option const & opt)
		{
			out("NAME", opt.name) << "from .%NAME% import %NAME%";
		}

		/// mentions an option in a family file
		void OptType(Formatter & out, Option const & opt)
		{
			out("NAME", opt.name) << "%NAME%,";
		}
	}

	/// generates wrappers for all options in the family and a file listing options in the family
	PyCtx const & operator << (PyCtx const & ctx, Family const & f)
	{
		ctx.create(ctx.dir(f));
		// generating option wrappers
		std::for_each(f.options.begin(), f.options.end(), boost::ref(ctx) << lm::_1);

		// creating family's option listing
		Formatter ff(ctx.filename(f));
		ff << (seq, 
			foreach_x(f.options, print::OptionUsing), "", 
			"def all(): return [", +foreach_x(f.options, print::OptType), "]");

		/* ctx.out(2) << "   " << f.name << std::endl; */

		return ctx;
	}

	namespace print 
	{
		/// prints import family clause for all families file
		void importFamily(Formatter & out, std::string const & family_name)
		{
			out("NAME", family_name) << "from .opt.%NAME% import %NAME%";
		}

		/// mentions the family in all families file
		void family(Formatter & out, std::string const & family_name)
		{
			out("NAME", family_name) << "%NAME%,";
		}
	}

	/// generates all options of all families and creates all families file
	PyCtx const& operator << (PyCtx const& ctx, Families const & fs)
	{
		ctx.out(1) << "Generating option families:...";

		ctx.create(ctx.opt());

		std::for_each(fs.families.begin(), fs.families.end(), ctx << lm::_1);
		
		std::list<std::string> names; // ugly copy but we cannot use a new boost =(
		
		BOOST_FOREACH(Families::Names::const_reference n, fs.names)
		{
		    names.push_back(n.first);
		}

		Formatter f(ctx.optionsPy());
		f << (seq, 
			foreach_x(names, print::importFamily), "", 
			"def all(): return [", +foreach_x(names, print::family), "]");

		/* ctx.out(1) << "ok!" << std::endl; */

		return ctx;
	}
}}}
