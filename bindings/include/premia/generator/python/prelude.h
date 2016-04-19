#pragma once

#include <boost/filesystem/fstream.hpp>
#include <premia/generator/python/ctx.h>

namespace premia {
namespace pygen {
namespace python {

	/// \brief represents fixed files to be generated first
	struct CommonPy
	{
		static std::string portable(fs::path const & t)
		{
			std::string s = t.string();
			boost::algorithm::replace_all(s, "\\", "/");
			return s;
		}

		friend PyCtx const & operator << (PyCtx const & ctx, CommonPy const & x)
		{
			Formatter cf(ctx.commonPy());
			cf.process_file(ctx.templateDir() / "common.py.template");

			Formatter ff(ctx.interopPy());
			ff	("PREMIA_DLL", portable(ctx.premiaDll()))
				("PREMIA_DATA", portable(ctx.premiaData()))
				.process_file(ctx.templateDir() / "interop.py.template");

			return ctx;
		}
	};

}}}
