#pragma once

#include <premia/generator/utils/formatter_dsl.h>
#include <premia/generator/api/model.h>
#include <premia/generator/python/ctx.h>

namespace premia {
namespace pygen {
namespace python {
	
	namespace print 
	{
		/// prints import clause in the pricings file for a model
		void importModelPricings(Formatter & out, Pricing const *p)
		{
			out("FAMILY_ID", p->family->name) << "from .%MODEL_ID%_%FAMILY_ID% import %MODEL_ID%_%FAMILY_ID%";
		}

		/// prints pricing name in the pricings file for a model
		void PricingType_2(Formatter & out, Pricing const *p)
		{
			out("FAMILY_ID", p->family->name) << "%MODEL_ID%_%FAMILY_ID%, ";
		}

		void ModelFamilies(Formatter & out, Pricing const *p)
		{
			out("FAMILY_ID", p->family->name) << "'%FAMILY_ID%', ";
		}
	}

	/// generates a wrapper for a model and a file with enumeration of all pricings for the model
	PyCtx const & operator << (PyCtx const & ctx, Model const & m)
	{
		ctx.create(ctx.dir(m));

		Formatter f(ctx.filename(m));
		
		f	("NAME", m.name)
			("ID",   m.ID())
			("ASSET_ID", id(m.asset))
			("ASSET_TYPE", m.asset.source->name)
			("MEMBERS_LEN", m.members.size())
			("MODEL_ID", id(m)) << (seq, 
				"from ...common import *",
				"",
				"class %NAME%(object):", +(seq, 
					"def __init__(self):", 
						+call(boost::bind(print::Initializers, _1, m.members)),
					"",
					"@staticmethod",
					"def ID(): return '%ID%'", 
					"",
					"@staticmethod",
					"def assetType(): ",
					"   from ...assets import %ASSET_TYPE%", 
					"   return %ASSET_TYPE%",
					"",
					"@staticmethod",
					"def families():", +(seq,
					    "from premia import options",
						"return [", +foreach_x(m.pricings, print::ModelFamilies), "]"
					), "",
					foreach_x(m.members, boost::bind(print::PropertyEx, _1, _2, boost::cref(m.members))), 
					"def __repr__(self): return getRepr(self, 'Model')", "", 
 					"def makeCurrent(self):", +(seq,
						"from premia import interop",
 						"interop.setCurrentAsset(%ASSET_ID%)",
 						"interop.setCurrentModel(%MODEL_ID%)",
						foreach_x(m.members, print::copy_param),
						"interop.stopWriteParameters()"
 						), "",
 					"def sync(self):", +(seq,
						"from premia import interop",
						"self.makeCurrent()",
						"interop.readCurrentModel()",
						foreach_x(m.members, print::load_param),
						"interop.stopReadParameters()"
 						), "",
					"@staticmethod",
	 				"def create(args, iterables):", +(seq,
	 					"self = %NAME%()",
 						"assert(len(args) == %MEMBERS_LEN%)",
						"it = args.__iter__()",
						foreach_x(m.members, print::assign_param),
						"return self"
 						), "",
					"@staticmethod",
					"def parameters(): ", +(seq, 
						"return [", +foreach_x(m.members, print::member), "]"),
					"def meta(self): ", +(seq, 
						"return [", +foreach_x(m.members, print::meta), "]")
				));

		Formatter ff(ctx.modelPricingsName(m));
		ff	("MODEL_ID", m.ID()) << (seq, 
			foreach_x(m.pricings, print::importModelPricings), "", 
			"def all(): return [", +foreach_x(m.pricings, print::PricingType_2), "]");

		/* ctx.out(2) << "   " << m.ID() << std::endl; */

		return ctx;
	}

	namespace print
	{
		/// prints import clause in all models file
		void ModelUsing(Formatter & out, Model const & m)
		{	
			if (!m.pricings.empty())
				out("ID", m.ID())("NAME", m.name) << "from .mod.%ID%.model import %NAME%";
		}

		/// mentions model name in all models file
		void ModelType(Formatter & out, Model const & m)
		{
			if (!m.pricings.empty())
				out("NAME", m.name) << "%NAME%,";
		}

		/// prints import clause in all pricings file
		void PricingUsing(Formatter & out, Model const & m)
		{
			if (!m.pricings.empty())
				out("ID", m.ID()) << "from .mod.%ID% import %ID%";
		}

		/// mentions pricings for the model in all pricings file
		void PricingType(Formatter & out, Model const & m)
		{
			if (!m.pricings.empty())
				out("ID", m.ID()) << "%ID%,";
		}
	}

	/// generates wrappers for all models, creates listing of all models and pricings
	PyCtx const & operator <<(PyCtx const & ctx, Models const & m)
	{
		ctx.out(1) << "Generating models:...";

		ctx.create(ctx.mod());

		// generating wrappers for all models
		std::for_each(m.models.begin(), m.models.end(), ctx << lm::_1);

		// creating list of all models
		Formatter f(ctx.modelsPy());
		f << (seq,
			foreach_x(m.models, print::ModelUsing), "",
			"def all(): return [", +foreach_x(m.models, print::ModelType), "]");

		// creating list of all pricings
		Formatter ff(ctx.pricingsPy()) ;
		ff << (seq,
			foreach_x(m.models, print::PricingUsing), "",
			"def all(): return [", +foreach_x(m.models, print::PricingType), "]");

		ctx.out(1) << "ok!" << std::endl;

		return ctx;
	}
}}}
