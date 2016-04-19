#ifdef _MSC_VER
#  define BOOST_FILESYSTEM_VERSION 3
#endif

#include <premia/import.h>

#include <istream>
#include <boost/filesystem.hpp>
#include <boost/lambda/lambda.hpp>


namespace fs = boost::filesystem;
typedef fs::path  path_t;

#include <premia/generator/utils/fspath.h>
#include <premia/generator/utils/symbols.h>
#include <premia/generator/utils/formatter.h>
#include <premia/generator/utils/formatter_dsl.h>
#include <premia/generator/python/ctx.h>

typedef boost::format fmt;
namespace lm = boost::lambda;

const std::string premia_lib_name = "premia";

typedef std::logic_error premia_exception;
 
#include <premia/generator/api/all.h>
#include <premia/generator/python/var.h>

#include <premia/generator/python/prelude.h>
#include <premia/generator/python/asset.h>
#include <premia/generator/python/enum.h>
#include <premia/generator/python/model.h>
#include <premia/generator/python/option.h>
#include <premia/generator/python/pricing.h>

#include <boost/program_options.hpp>
#include <boost/foreach.hpp> 

/* Get shared library extension in SHEXT */
//#include "../../Src/config.h"


using premia::pygen::formatter_dsl::Formatter;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    using namespace premia::pygen;

	try {
    
		FsPath current_path = fs::current_path();
		FsPath root         = current_path->parent_path().parent_path();
		FsPath dll_dir;
		FsPath template_dir;
		FsPath data_dir;
		FsPath output_path;
    
		int verbosity;

		po::options_description desc("Allowed options");
		desc.add_options()
			("help", "produce help message")
			("output-dir", po::value(&output_path)->default_value(*root / "bindings" / "python" / "premia"), "output directory")
			("dll-dir",    po::value(&dll_dir)->default_value(*current_path), "directory containing python binding library")
			("data-dir",   po::value(&data_dir)->default_value(*root / "data"), "directory containing premia data files")
			("template-dir", po::value(&template_dir)->default_value(*current_path / "templates"), "directory containing templates for pygen")
			("verbosity,v", po::value(&verbosity)->default_value(1), "verbosity level (0 - no output, 1 - basic, 2 - detailed)")
		;
    
		po::variables_map vm;
		po::store(po::parse_command_line(argc, argv, desc), vm);
		po::notify(vm); 
    
		// making them absolute
		data_dir = !data_dir->root_directory().empty() ? *data_dir : fs::current_path() / *data_dir;

		python::PyCtx ctx(*data_dir, *dll_dir, *template_dir, *output_path, verbosity);
    
		ctx.out(1)
				<< "premia-root: "  << root         << std::endl
				<< "output-dir: "   << output_path  << std::endl
				<< "dll-dir: "      << dll_dir      << std::endl
				<< "data-dir: "     << data_dir     << std::endl
				<< "template-dir: " << template_dir << std::endl;

		if (vm.count("help")) {
			std::cout << desc << "\n";
			return 1;
		}    

		if (!fs::exists(*output_path))
			fs::create_directories(*output_path);
    
		python::createDir(*output_path);

		InitVar();

		strcpy(premia_data_dir, data_dir->string().c_str());
	
		ctx.out(1) << "Initializing...";

		PremiaCtx p;

		ctx.out(1) << "ok!" << std::endl;
    
		ctx << p.models << p.families << p.pricings << p.enums << p.assets;
	} catch (std::exception const &e) {
		std::cout << "exception caught: " << e.what() << std::endl;
	}   
	return 0;
}
