#pragma once

#include <premia/generator/api/option.h>
#include <premia/generator/fsharp/ctx.h>
#include <premia/generator/fsharp/var.h>

namespace premia {
namespace pygen {
namespace fsharp {

	/// generates an F# wrapper for the option
	FsCtx const & generateOption (FsCtx const & ctx, Option const & opt)
	{
		Formatter f(ctx.filename(opt));
		f	("NAME",	  opt.name)
			("OPTION_ID", id(opt))
			("ASSET_ID",  id(opt.family.asset))
			("FAMILY_ID", id(opt.family))
			("FAMILY_NAME", opt.family.name)
			<< (seq, 
				"namespace FsPremia.Types.opt.%FAMILY_NAME%", "", +(seq, 
					"open FsPremia", "open FsPremia.Interop", "open FsPremia.Util", "",
					"type %NAME% = ", +(seq, 
						block("{}", foreach_x(opt.vars, print::memberDecl)),  
						"with", "interface IPremiaObj with", +(seq, 
							"member x.makeCurrent() = ", +(seq,
 								"setCurrentAsset(%ASSET_ID%)",
 								"setCurrentOption(%FAMILY_ID%, %OPTION_ID%)",
								foreach_x(opt.vars, print::copy_param_fs),
								"stopWriteParameters()"
 						)), "",
						"static member Create() = ",
							+block("{}", foreach_x(opt.vars, print::memberIni)))));

		return ctx;
	}

	namespace print {

		/// prints shortcut to the option constructor
		void optionUsing(Formatter & out, Option const & opt)
		{
			out("NAME", opt.name) << "let %NAME% = FsPremia.Types.opt.%FAMILY%.%NAME%.Create()";
		}
	}

	/// generates wrappers for all options in the family and a file with shortcuts to their constructors
	FsCtx const & generateFamily (FsCtx const & ctx, Family const & f)
	{
		std::for_each(f.options.begin(), f.options.end(), boost::bind(generateOption, boost::ref(ctx), _1));

		Formatter ff(ctx.filename(f));
		ff	("FAMILY", f.name) << (seq, 
				"namespace FsPremia.opt", +(seq, "",
					"module %FAMILY% = ", 
						+foreach_x(f.options, print::optionUsing)));

		ctx.out(2) << "   " << f.name << std::endl;

		return ctx;
	}

	/// generates wrappers for all options for all families
	FsCtx const& operator << (FsCtx const& ctx, Families const & fs)
	{
		ctx.out(1) << "Generating option families:...";

		std::for_each(fs.families.begin(), fs.families.end(), boost::bind(generateFamily, boost::ref(ctx), _1));

		ctx.out(1) << "ok." << std::endl;

		return ctx;
	}

}}}
