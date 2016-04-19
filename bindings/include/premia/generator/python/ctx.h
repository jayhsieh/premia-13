#pragma once

#include <set>
#include <list>
#include <string>
#include <boost/filesystem.hpp>
#include <premia/generator/api/all.h>
#include <premia/generator/utils/null_stream.h>
#include <premia/generator/karrigell/result_kinds.h>

namespace premia {
namespace pygen  {
namespace python {

	using karrigell::ResultKindsInitialized;

	namespace fs = boost::filesystem;

	/// creates a python module in a given directory
	/// \param p path to the directory, its parent directory should exist
	inline void createDir(fs::path const & p)
	{
		if (!fs::exists(p))  
			fs::create_directory(p);
		fs::ofstream(p / "__init__.py");
	}

	struct CommonPy;
	
	
	/// \brief A context responsible for generated python wrappers directory structure
	struct PyCtx
	{
		/// constructs a context
		/// \param base_path path to a directory where wrappers will be generated
		PyCtx(fs::path const & data_dir, 
		      fs::path const & dll_dir,
		      fs::path const & template_dir,
		      fs::path const & output_dir,
		      int verbosity) 
			: data_dir     (data_dir)
			, output_dir   (output_dir)
			, dll_dir      (dll_dir)
			, template_dir (template_dir)
			, verbosity    (verbosity)
		{}
		
		std::ostream& out(int v) const 
		{
	   	   static null_stream nullstr;

		   if (verbosity < v) return nullstr;
		   else               return std::cout;
		}

		/// returns path to Premia Runtime DLL
		fs::path premiaDll() const { return dll_dir; }		
		
		fs::path templateDir() const { return template_dir; }
		
		int getVerbosity() const { return verbosity; }

		/// returns path to Premia /data directory
		fs::path premiaData() const { return data_dir; } 

		/// returns path to common.py
		fs::path commonPy() const {	return output_dir / "common.py"; }

		/// returns path to interop.py
		fs::path interopPy() const { return output_dir / "interop.py"; } 

		/// returns path to the directory where enumerations will be generated
		fs::path enumDir() const { return output_dir / "enum"; } 

		/// creates a python directory
		/// \param p path to directory
		void create(fs::path const & p) const 	{	createDir(p);	}

		/// adds .py extension to a string
		static std::string py(std::string const & p) { return p + ".py"; }

		/// \param e a reference to a enum type
		/// returns path to a file where the enumeration should be generated
		fs::path filename(Enum const & e) const { 	return enumDir() / py(e.label);	}

		fs::path assetsPy() const { return output_dir /  "assets.py"; }

		/// returns a path to the directory where models and pricings will be generated
		fs::path mod() const {	return output_dir / "mod"; 	}

		/// returns a path to the directory where the model and its pricing methods will be generated
		fs::path dir(Model const & m) const {	return mod() / m.ID();		}

		/// returns path to a file where the model should be generated
		fs::path filename(Model const & m) const { 	return dir(m) / "model.py"; }

		/// returns path to file with all models
		fs::path modelsPy() const { return output_dir / "models.py"; }

		/// returns path to directory where options will be generated
		fs::path opt() const { return output_dir / "opt"; }

		/// returns path to directory where the family options will be generated
		fs::path dir(Family const & f) const {	return opt() / f.name;	}

		/// returns path to file where the the option will be generated
		fs::path filename(Option const & opt) const	{	return dir(opt.family) / py(opt.name);	}

		/// returns path to a file where a description of the family will be generated
		fs::path filename(Family const & f) const	{	return dir(f) / py(f.name);	}

		/// returns path to a file where a description of all options will generated
		fs::path optionsPy() const { return output_dir / "options.py"; }

		/// returns a label for the pricing
		static std::string name(Pricing const & p) { return p.model->ID() + "_" + p.family->name; }

		/// returns path to a directory where methods of the pricing will be generated
		fs::path dir(Pricing const & p) const {	return dir(*p.model) / name(p);	}

		/// returns path to a file where a description of the family will be generated
		fs::path filename(Pricing const & p) const	{ return dir(p) / py(name(p)); }

		/// returns path to a file where the pricing method will be generated
		fs::path filename(PricingMethod const & m) const { return dir(m.pricing) / py(m.name); }

		/// returns path to a file where description of all pricings for the given method will generated
		fs::path modelPricingsName(Model const & m) const { return dir(m) / py(m.ID()); }

		/// returns path to a file where a description of all pricings will be generated
		fs::path pricingsPy() const { return output_dir / "pricings.py"; }

		ResultKindsInitialized const & resultKinds() const { return result_kinds_; }

	private:
		fs::path const data_dir; //!< path to Premia data directory 
		fs::path const output_dir;	//!< root directory for the library being generated
		fs::path const dll_dir; // !< directory containing premia binding library
		fs::path const template_dir;
		int const verbosity;
		ResultKindsInitialized const result_kinds_;
	};
}}}
