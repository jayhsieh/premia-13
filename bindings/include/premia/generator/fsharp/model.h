#pragma once

#include <premia/generator/api/model.h>
#include <premia/generator/fsharp/ctx.h>
#include <premia/generator/fsharp/var.h>

namespace premia { 
namespace pygen {
namespace fsharp {

	/// generates F# wrapper for a model
	FsCtx const & generateModel (FsCtx const & ctx, Model const & m)
	{
		Formatter f(ctx.filename(m));
		f	("NAME", m.name)
			("ASSET_ID", id(m.asset))
			("MODEL_ID", id(m)) << (seq, 
				"namespace FsPremia.Types.mods.%NAME%", "", +(seq, 
					"open FsPremia", "open FsPremia.Interop", "open FsPremia.Util", "",
					"type Model = ", +(seq, 
						block("{}", foreach_x(m.members, print::memberDecl)),  
						"with", "interface IPremiaObj with", +(seq, 
							"member x.makeCurrent() = ", +(seq,
								"setCurrentAsset(%ASSET_ID%)",
								"setCurrentModel(%MODEL_ID%)",
								foreach_x(m.members, print::copy_param_fs),
								"stopWriteParameters()"
						)), "",
						"static member Create() = ",  
							+block("{}", foreach_x(m.members, print::memberIni)))));

		ctx.out(2) << "   " << m.name << std::endl;

		return ctx;
	}

	namespace print 
	{
		/// prints a shortcut to a model constructor
		void modelUsing(Formatter & out, Model const & m)
		{
			out("NAME", m.name) << "let %NAME% = FsPremia.Types.mods.%NAME%.Model.Create()";
		}
	}

	/// generates wrappers for all Premia models and creates shortcuts to their constructors
	FsCtx const & operator <<(FsCtx const & ctx, Models const & m)
	{
		ctx.out(1) << "Generating models:...";

		std::for_each(m.models.begin(), m.models.end(), boost::bind(generateModel, boost::ref(ctx), _1));

		ctx.out(1) << "ok!" << std::endl;

		Formatter f(ctx.modelsFilename());
		f << (seq, 
			"namespace FsPremia", +(seq, "",
				"module models = ", 
					+foreach_x(m.models, print::modelUsing)
				)
			);

		return ctx;
	}

}}}
