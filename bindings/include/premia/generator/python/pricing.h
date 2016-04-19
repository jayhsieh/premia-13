#pragma once

namespace premia {
namespace pygen {
namespace python {

	namespace print {

		/// prints an option name compatible with the method 
		void Compatible(Formatter & out, Option const * opt)
		{
			out("OPT_NAME", opt->name) << "%OPT_NAME%,";
		}

		/// returns code for getting result from a method according to the type of result parameter
		const char* getResultFmt(int var_type)
		{
			switch (var_type)
			{
			case BOOL: return "('%RES_NAME%', interop.getResultBool(%RES_IDX%)),";
			case DOUBLE: return "('%RES_NAME%', interop.getResultDouble(%RES_IDX%)),";
			case PNLVECT: return "('%RES_NAME%', getResultArray(%RES_IDX%)),";
			default:
				throw std::logic_error(("Unexpected result type " + str(var_type)).c_str());
			}
		}
		
		/// prints code for getting result parameter
		/// \param e tuple (variable, its index in result parameters)
		void ResElementEx(Formatter & out, boost::tuple<NamedVar, int> const & e)
		{
			NamedVar const & v = boost::get<0>(e);

			out("RES_NAME", v.name)("RES_IDX", boost::get<1>(e)) << getResultFmt(v.src->Vtype);
		}
		
		void OutParamName(PyCtx const &ctx, Formatter & out, NamedVar const & v)
		{
			karrigell::ResultKinds::Kind k = ctx.resultKinds().lookup(v.name);
			out("RES_NAME", v.src->Vname)
				("KIND", k.second)
				("VISIBLE", k.first ? "True" : "False") 
				<< 
					"('%RES_NAME%', %KIND%, %VISIBLE%),";
		}
		
	}

	inline void printResults(Formatter & out, VarList const & results)
	{
		int i = 0;
		BOOST_FOREACH(NamedVar const vr, results)
		{
			print::ResElementEx(out, boost::make_tuple(vr, i++));	
		}
	}        

	/// generates wrapper for a pricing method
	PyCtx const & operator << (PyCtx const & ctx, PricingMethod const & method)
	{
		try {
			Formatter f(ctx.filename(method));
			f  ("CLASS", method.name)
			   ("MODEL_NAME", method.pricing.model->name)
			   ("FAMILY_ID", method.pricing.family->name)
			   ("ASSET_ID", id(method.pricing.asset))
			   ("PRICING_ID", id(method.pricing))
			   ("METHOD_ID", id(method))
				("MEMBERS_LEN", method.members.size())
				<< (seq, 
				"from ..model import %MODEL_NAME%", 
				"from ....opt.%FAMILY_ID%.%FAMILY_ID% import *",
				"from ....common import *", ""
				"class %CLASS%(object):", +(seq, 
					"def __init__(self):", 
						+call(boost::bind(print::Initializers, _1, boost::ref(method.members))), "", 
					foreach_x(method.members, boost::bind(print::PropertyEx, _1, _2, boost::cref(method.members))), "", 
					"def __repr__(self): return getRepr(self, 'Pricing Method')", "" 
	 				"def makeCurrent(self):", +(seq,
						"from premia import interop",
	 					"interop.setCurrentAsset(%ASSET_ID%)",
	 					"interop.setCurrentMethod(%PRICING_ID%, %METHOD_ID%)",
						foreach_x(method.members, print::copy_param),
	 					"interop.stopWriteParameters()"
	 					), "",
	 				"def compute(self, opt, mod):", +(seq,
						"from premia import interop",
	 					"mod.makeCurrent()", 
	 					"opt.makeCurrent()",
	 					"self.makeCurrent()",
						"interop.compute_3()",
						"return [", 
							+call(boost::bind(printResults, _1, boost::cref(method.results))),
						"]"), "",

					"@staticmethod",
					"def results(): ", +(seq, 
						"return [", +foreach_x(method.results, boost::bind(print::OutParamName, boost::cref(ctx), _1, _2)), "]"), "",
					"@staticmethod",
					"def parameters(): ", +(seq, 
						"return [", +foreach_x(method.members, print::member), "]"), "",
					"@staticmethod",
	 				"def create(args, iterables):", +(seq,
	 					"self = %CLASS%()",
	 					"assert(len(args) == %MEMBERS_LEN%)",
						"it = args.__iter__()",
						foreach_x(method.members, print::assign_param),
						"return self"
	 					), "",
					"def meta(self): ", +(seq, 
						"return [", +foreach_x(method.members, print::meta), "]"),
					"@staticmethod",
					"def model(): return %MODEL_NAME%", ""
					"@staticmethod",
					"def options():", +(seq, 
						"return [", +foreach_x(method.compatible_options, print::Compatible), "]"), "",
					"def __call__(self, option, model): return self.compute(option, model)"
				));
		}
		catch (...) {
			ctx.out(1) << "When processing " 
				<< method.pricing.model->name << "." 
				<< method.pricing.family->name << "." 
				<< method.name;
			throw;
		}

		return ctx;
	}

	namespace print {
		
		/// prints import method clause in a pricing description file
		void importMethod(Formatter & out, PricingMethod const & m)
		{
			out("NAME", m.name) << "from .%NAME% import %NAME%";
		}

		/// mentions method in a pricing description file
		void Method(Formatter & out, PricingMethod const & m)
		{
			out("NAME", m.name) << "%NAME%,";
		}

		/// mentions method in a pricing description file
		void pMethod(Formatter & out, PricingMethod const * m)
		{
			out("NAME", m->name) << "%NAME%,";
		}

		void MethodsForOptions(Formatter & out, std::pair<Option const *, std::list<PricingMethod*> > const & p)
		{
			out("OPT_NAME", p.first->name) << (seq, "%OPT_NAME% : [", +foreach_x(p.second, pMethod), "],");
		}
	}

	/// generates wrappers for all methods in a pricing and creates a file with the pricing description
	PyCtx const & operator << (PyCtx const & ctx, Pricing const & p)
	{
		ctx.create(ctx.dir(p));

		// generating wrappers for methods
		std::for_each(p.methods.begin(), p.methods.end(), boost::ref(ctx) << lm::_1);

		// creating a file with the pricing description
		Formatter f(ctx.filename(p));
		f	("FAMILY_NAME",p.family->name) << (seq, 

			foreach_x(p.methods, print::importMethod), "",
			"from ....opt.%FAMILY_NAME%.%FAMILY_NAME% import *", "",
			"def all(): return [", +foreach_x(p.methods, print::Method), "]", "",
			"def methods_for_options(): return {", +foreach_x(p.methods_for_options, print::MethodsForOptions), "}");

		/* ctx.out(2) << "   " << (*p.source)->ID << std::endl; */

		return ctx;
	}

	/// generates all methods of all pricings 
	PyCtx const & operator << (PyCtx const & ctx, Pricings const & p)
	{
		ctx.out(1) << "Generating pricing methods:...";
		std::for_each(p.pricings.begin(), p.pricings.end(), boost::ref(ctx) << lm::_1);
		/* ctx.out(1) << "ok!" << std::endl; */
		return ctx;
	}
}}}
