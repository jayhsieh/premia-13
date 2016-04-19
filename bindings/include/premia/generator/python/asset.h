#pragma once 

namespace premia {
namespace pygen  {
namespace python {

	inline void printAssetModel(Formatter & out, Model const * m)
	{
		if (!m->pricings.empty())
			out("MODEL_NAME", m->name) << "%MODEL_NAME%,";
	}

	inline void printAsset (Formatter & out, Asset const & a)
	{
		typedef std::map<int, api::Model const*> Sorted;
		Sorted sorted;
		BOOST_FOREACH(api::Model const * m, a.models)
			sorted[id(*m)] = m;
		std::vector<api::Model const*> models;
		BOOST_FOREACH(Sorted::const_reference p, sorted)
			models.push_back(p.second);

		out("ASSET_NAME", a.source->name) << (seq, 
			"class %ASSET_NAME%(object):", +(seq, 
				"@staticmethod",
				"def models(): ", +(seq, 
					"return [", +foreach_x(models, printAssetModel), "]"), ""
			));
	}

	inline void printAssetName(Formatter & out, Asset const & a)
	{
		out("ASSET_NAME", a.source->name) << "%ASSET_NAME%, ";
	}

	/// generates correspondence between Premia assets and models
	PyCtx const & operator << (PyCtx const & ctx, Assets const & a) 
	{
		Formatter f(ctx.assetsPy());
		f << (seq, 
			"from models import *", "",
			foreach_x(a.assets, printAsset),
			"def all(): return [", +foreach_x(a.assets, printAssetName), "]" 
		);

		return ctx;
	}

}}}
