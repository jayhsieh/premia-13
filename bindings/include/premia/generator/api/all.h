#pragma once

#include <premia/generator/api/model.h>
#include <premia/generator/api/option.h>
#include <premia/generator/api/pricing.h>

namespace premia {
namespace pygen  {
namespace api    {

	struct PremiaCtx : boost::noncopyable
	{
		Assets		assets;
		Enums		enums;
		Models		models;
		Families	families;
		Pricings	pricings;

		PremiaCtx()
			:	models(*this)
			,   families(*this)
			,	pricings(*this)
		{
		}
	};
}

using api::PremiaCtx;

}}