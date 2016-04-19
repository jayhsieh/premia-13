#pragma once

#include <premia/generator/api/enum.h>
#include <premia/generator/utils/formatter_dsl.h>

namespace premia {
namespace pygen	{
namespace fsharp{

	namespace print {

		/// binds the formatter with the label of the enum member and its ordinal value
		inline Formatter & operator && (Formatter & out, Enum::Members::const_reference p)
		{
			return out("LABEL", p.second.label)("ID", p.first);
		}

		/// prints initialization clause for the enum member
		void Ini(Formatter & out, Enum::Members::const_reference p)
		{
			(out && p) << "| %LABEL% = %ID%";
		}
	}

	/// generates a wrapper for an enumeration type
	FsCtx const & generateEnum (FsCtx const & ctx, Enum const & e)  
	{
		Formatter(ctx.filename(e))
			("CLASS", e.label) << (seq,
				"namespace FsPremia", "", +(seq,
				"type %CLASS% = ", 
					+foreach_x(e.members, print::Ini)));

		ctx.out(2) << "   " << e.label << std::endl;

		return ctx;
	}
	/// generates wrappers for all enumerations
	FsCtx const & operator << (FsCtx const & ctx, Enums const & e) 
	{
		ctx.out(1) << "Generating enums:...";

        BOOST_FOREACH(Enums::EnumsToLabels::const_reference r, e.enums())
        {
            generateEnum(ctx, r.second);
        } 

		ctx.out(1) << "ok!" << std::endl;
		return ctx;
	}

}}}
