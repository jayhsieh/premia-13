#pragma once

#include <premia/generator/api/pricing.h>
#include <premia/generator/fsharp/ctx.h>
#include <premia/generator/fsharp/var.h>

namespace premia {
namespace pygen {
namespace fsharp {

	namespace print {

		/// \return a code retrieving from Premia value of a result parameter
		const char* getResultFmt(int var_type)
		{
			switch (var_type)
			{
			case BOOL: return "%RES_NAME% = get_result_bool(%RES_IDX%)";
			case DOUBLE: return "%RES_NAME% = get_result_double(%RES_IDX%)";
			case PNLVECT: return "%RES_NAME% = getResultArray(%RES_IDX%)";
			default:
				throw std::logic_error("Unexpected result type " + str(var_type));
			}
		}

		/// \return string representation for F# type of a result VAR
		const char* getResultType(int var_type)
		{
			switch (var_type)
			{
			case BOOL: return "bool";
			case DOUBLE: return "float";
			case PNLVECT: return "float list";
			default:
				throw std::logic_error("Unexpected result type " + str(var_type));
			}
		}

		/// prints a code for retrieving i-th result parameter from a Premia method
		void ResElementEx(Formatter & out, boost::tuple<NamedVar, int> const & e)
		{
			NamedVar const & v = boost::get<0>(e);

			out("RES_NAME", v.name)("RES_IDX", boost::get<1>(e)) << getResultFmt(v.src->Vtype);
		}

		/// prints a declaration for the result VAR 
		void ResElementDecl(Formatter & out, boost::tuple<NamedVar, int> const & e)
		{
			NamedVar const & v = boost::get<0>(e);

			out("RES_NAME", v.name)("RES_TYPE", getResultType(v.src->Vtype)) << "%RES_NAME% : %RES_TYPE%";
		}

		/// prints compute method for a compatible option
		void compute(Formatter & out, Option const * opt)
		{
			out("OPT_NAME", opt->name) << "member x.compute(opt : %OPT_NAME%, m : Model) = x.compute'(opt,m)";
		}

		/// prints makeCurrent method 
		void makeCurrent(Formatter & out, PricingMethod const & method)
		{
			out << (seq,
				"with", "interface IPremiaObj with", +(seq, 
					"member x.makeCurrent() = ", +(seq,
						"setCurrentAsset(%ASSET_ID%)",
						"setCurrentMethod(%PRICING_ID%, %METHOD_ID%)",
						foreach_x(method.members, print::copy_param_fs),
						"stopWriteParameters()"
				)), "");
		}
		
	    inline void printResults(Formatter & out, VarList const & results)
	    {
		    int i = 0;
		    BOOST_FOREACH(NamedVar const vr, results)
		    {
			    print::ResElementEx(out, boost::make_tuple(vr, i++));	
		    }
	    }        
	    
		/// prints generic compute method
		void computeGeneric(Formatter & out, PricingMethod const & method)
		{
			out << (seq, 
				"member private x.compute'(opt : #IPremiaObj, model : #IPremiaObj) = ", +(seq,
					"Util.premiaCompute(model, opt, x)", 
					block("{}",
						call(boost::bind(printResults, _1, boost::cref(method.results)))
					)));
		}

		/// prints common part for methods with and without parameters
		void common(Formatter & out, PricingMethod const & method)
		{
			makeCurrent(out, method);
			computeGeneric(out, method);
			out << foreach_x(method.compatible_options, print::compute);
		}

		/// prints import clauses 
		void imports(Formatter & out)
		{
			out << (seq, 
				"open FsPremia", "open FsPremia.Interop", "open FsPremia.Util", "",
				"open FsPremia.Types.mods.%MODEL_NAME%", "open FsPremia.Types.opt.%FAMILY_NAME%", ""
				);
		}

		/// prints a class containing result parameters for the method
		void result(Formatter & out, PricingMethod const & method)
		{
			out << (seq, 
				"and %CLASS%_Result = ", 
					+block("{}", foreach_x(method.results, print::ResElementDecl))
				);
		}
	}

	/// generates an F# wrapper for a pricing method
	FsCtx const & generateMethod (FsCtx const & ctx, PricingMethod const & method)
	{
		if (method.members.empty())
		{
			Formatter f(ctx.filename(method));
			f	("CLASS", method.name)
				("MODEL_NAME", method.pricing.model->name)
				("FAMILY_NAME", method.pricing.family->name)
				("ASSET_ID", id(method.pricing.asset))
				("PRICING_ID", id(method.pricing))
				("METHOD_ID", id(method)) << (seq,
					"namespace FsPremia.Types.mods.%MODEL_NAME%.%FAMILY_NAME%", "", +(seq, 
						call(print::imports),
						"type %CLASS%() = ", +(seq, 
							"class", "end",  
							call(boost::bind(print::common, _1, boost::ref(method))),
							"static member Create() = %CLASS%()"),
						call(boost::bind(print::result, _1, boost::ref(method)))));

		}
		else {

			Formatter f(ctx.filename(method));
			f	("CLASS", method.name)
				("MODEL_NAME", method.pricing.model->name)
				("FAMILY_NAME", method.pricing.family->name)
				("ASSET_ID", id(method.pricing.asset))
				("PRICING_ID", id(method.pricing))
				("METHOD_ID", id(method)) << (seq,
					"namespace FsPremia.Types.mods.%MODEL_NAME%.%FAMILY_NAME%", "", +(seq, 
						call(print::imports),
						"type %CLASS% = ", +(seq, 
							block("{}", foreach_x(method.members, print::memberDecl)),  
							call(boost::bind(print::common, _1, boost::ref(method))),
							"static member Create() = ",  
								+block("{}", foreach_x(method.members, print::memberIni))),
							call(boost::bind(print::result, _1, boost::ref(method)))));
		}
		return ctx;
	}

	namespace print {

		/// prints a shortcut to the method constructor
		void methodUsing(Formatter & out, PricingMethod const & m)
		{
			out("NAME", m.name) << "let %NAME% = FsPremia.Types.mods.%MODEL_NAME%.%FAMILY_NAME%.%NAME%.Create()";
		}
	}

	/// generates wrappers for all methods in the pricing
	FsCtx const & generatePricing (FsCtx const & ctx, Pricing const & p)
	{
		if (!p.methods.empty())
		{
			std::for_each(p.methods.begin(), p.methods.end(), boost::bind(generateMethod, boost::ref(ctx), _1));

			Formatter f(ctx.filename(p));
			f	("MODEL_NAME", p.model->name)
				("FAMILY_NAME", p.family->name) << (seq, 
					"namespace FsPremia.methods.%MODEL_NAME%", +(seq, "",
						"module %FAMILY_NAME% = ", 
							+foreach_x(p.methods, print::methodUsing)));

			ctx.out(2) << "   " << p.model->name << "_" << p.family->name << std::endl;
		}

		return ctx;
	}

	/// generates wrappers for all methods in all pricings
	FsCtx const & operator << (FsCtx const & ctx, Pricings const & p)
	{
		ctx.out(1) << "Generating pricing methods:...";
		std::for_each(p.pricings.begin(), p.pricings.end(), boost::bind(generatePricing, boost::ref(ctx), _1));
		ctx.out(1) << "ok!" << std::endl;
		return ctx;
	}

}}}
