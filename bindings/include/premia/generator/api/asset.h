#pragma once

#include <map>
#include <set>
#include <boost/ptr_container/ptr_list.hpp>
#include <boost/filesystem.hpp>
#include <premia/import.h>

namespace premia {

    namespace fs = boost::filesystem;

namespace pygen  {
namespace api {

	struct Model;

	/// \brief a wrapper over a native Premia asset
	struct Asset : boost::noncopyable
	{
		/// \brief creates an asset by a pointer to a native asset
		Asset(::PremiaAsset * source) : source(source) {}

		std::set<Model const *> models;

		/// \brief a pointer to the native asset
		::PremiaAsset * source;
	};

	/// \brief represents all Premia assets
	struct Assets 
	{
		/// \brief a list of Premia assets
		boost::ptr_list<Asset>	assets;

		/// \brief retrieves a wrapper for the given native asset
		Asset & lookup(::PremiaAsset * src)
		{
			// checking if we created a wrapper for it ones
			std::map< ::PremiaAsset*, Asset*>::iterator it = source_to_asset.find(src);

			// if no,
			if (it == source_to_asset.end())
			{
				// create a wrapper and put it into the map
				Asset * a = new Asset(src);
				assets.push_back(a);
				source_to_asset[src] = a;
				return *a;
			}

			return *it->second;
		}

	private:
		/// \brief maps from native assets to wrappers over them
		std::map<PremiaAsset*, Asset*>	source_to_asset;
	};

	/// \brief gets the asset id (0,1,2...)
	inline int id(Asset const & a)
	{
		return  a.source - premia_assets;
	}

}

using api::Asset;
using api::Assets;
using api::id;

}}
