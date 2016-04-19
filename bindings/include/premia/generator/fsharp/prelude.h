#pragma once

#include <boost/filesystem/fstream.hpp>
#include <premia/generator/fsharp/ctx.h>

namespace premia {
namespace pygen {
namespace fsharp {

	struct Common {

		/// generates premia_interop.fs and common.fs files
		friend FsCtx & operator << (FsCtx & ctx, Common const & )
		{
			Formatter f(ctx.dllimportFs()) ;
			f("DLLPATH", ctx.premiaDLL())
			    .process_file(ctx.templateDir() / "premia_interop.fs.template");

			Formatter ff(ctx.commonFs());
			ff ("DATA_PATH", ctx.dataDir())
			    .process_file(ctx.templateDir() / "common.fs.template");
			
			return ctx;
		}
		
	};

	namespace print {
		/// prints an include clause into F# project file
		void FsProjFile(Formatter & out, fs::path const & p)
		{
			out("FILE", p.string()) << "    <Compile Include=\"%FILE%\" />";
		}
	}

	struct FsProj {

		/// creates F# project file for Premia wrappers
		friend FsCtx const & operator << (FsCtx const & ctx, FsProj const &)
		{
			Formatter f(ctx.fsprojName()) ;
			(f	  .process_file(ctx.templateDir() / "FsPremia.fsproj.template")
				<< foreach_x(ctx.generatedFilesList(), print::FsProjFile))
				  .process_file(ctx.templateDir() / "FsPremia_end.fsproj.template");

			return ctx;
		}
	};

}}}
