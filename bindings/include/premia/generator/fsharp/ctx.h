#pragma once

#include <set>
#include <list>
#include <string>
#include <boost/filesystem/path.hpp>
#include <premia/generator/api/all.h>
#include <premia/generator/utils/null_stream.h>

namespace premia {
namespace pygen {
namespace fsharp {

	namespace fs = boost::filesystem;

	inline std::string operator - (std::string const & lhs, std::string const & rhs){ return lhs + "_" + rhs; }

	/// \brief context responsible for generated F# library structure on disk
	struct FsCtx : boost::noncopyable
	{
		/// constructs the context
		/// \param output_dir path to the directory where Premia is installed (so it contains directories data, src etc.)
		FsCtx(fs::path const & output_dir, fs::path const &data_dir, fs::path const &dll_path, fs::path const &template_dir, int verbosity) 
			:	output_dir_(output_dir)
			,	dll_path_(dll_path)
			, data_dir_(data_dir)
			, template_dir_(template_dir)
			, verbosity_(verbosity)
		{}
		
		std::ostream& out(int v) const 
		{
	   	   static null_stream nullstr;

		   if (verbosity_< v) return nullstr;
		   else               return std::cout;
		}
		
		fs::path templateDir() const { return template_dir_; }

		/// returns path to Premia runtime DLL
		fs::path premiaDLL() const { return dll_path_; }

		/// returns path to 'data' directory of premia installation
		fs::path dataDir() const { return data_dir_; }

		/// returns path to common.fs file
		fs::path commonFs() const {	return fs("common"); }

		/// returns path to FsPremia.fs file
		fs::path dllimportFs() const { return fs("FsPremia"); }

		/// returns path to an enum wrapper
		fs::path filename(Enum const & e) const { return fs("enum" - e.label); }

		/// returns path to a model wrapper
		fs::path filename(Model const & m) const { return fs("mod_" - m.ID()); }

		/// returns path to all models file
		fs::path modelsFilename() const { return fs("models");	}

		/// returns path to an option wrapper
		fs::path filename(Option const & opt) const { return fs("opt" - opt.family.name - opt.name);	}

		/// returns path to a file listing options of the family
		fs::path filename(Family const & f) const {	return fs("family" - f.name);	}

		/// returns a label for the pricing
		static std::string name(Pricing const & p) { return p.model->ID() - p.family->name; }

		/// returns path to a pricing method wrapper
		fs::path filename(PricingMethod const & m) const {	return fs("mod" - name(m.pricing) - m.name);	}

		/// returns path to a pricing wrapper
		fs::path filename(Pricing const & p) const { return fs("pricing" - name(p)); }

		/// returns list of file that were generated
		std::list<fs::path>  const & generatedFilesList() const  {	return generated_files;	}

		/// returns path to generated F# project file
		fs::path fsprojName() const {	return output_dir_ / "FsPremia.fsproj";	}

	private:
		/// directory where F# library will be generated
		//fs::path const base_path_;
		/// path to Premia installation
		fs::path const output_dir_;
		fs::path const dll_path_;
		fs::path const data_dir_;
		fs::path const template_dir_;
		int const      verbosity_;

		/// list of generated F# source files
		mutable std::list<fs::path>	generated_files;

		/// list of generated F# source files as a set in order to speed up lookup
		mutable std::set<fs::path>  generated_files_set;

		/// adds F# extension to the filename and makes it from relative path to absolute path
		/// also relative filename is saved to generated_files
		fs::path fs(std::string const & s) const
		{
			fs::path p = s + ".fs";

			if (generated_files_set.find(p) == generated_files_set.end())
			{
				generated_files_set.insert(p);
				generated_files.push_back(p);
			}
			return output_dir_ / p;
		}
	};
}}}

