#include <istream>
#include <boost/filesystem.hpp>
#include <boost/lambda/lambda.hpp>


namespace fs = boost::filesystem;
typedef fs::path  path_t;

#include <premia/import.h>
#include <premia/generator/utils/symbols.h>
#include <premia/generator/utils/formatter.h>
#include <premia/generator/utils/formatter_dsl.h>
#include <premia/generator/fsharp/ctx.h>

typedef boost::format fmt;
namespace lm = boost::lambda;

const std::string premia_lib_name = "premia";

typedef std::logic_error premia_exception;
 
#include <premia/generator/api/all.h>
#include <premia/generator/fsharp/var.h>

#include <premia/generator/fsharp/prelude.h>
#include <premia/generator/fsharp/enum.h>
#include <premia/generator/fsharp/model.h>
#include <premia/generator/fsharp/option.h>
#include <premia/generator/fsharp/pricing.h>

#include <boost/program_options.hpp>
#include <boost/foreach.hpp>

using premia::pygen::formatter_dsl::Formatter;

namespace po = boost::program_options;

int main(int argc, char *argv[])
{
    using namespace premia::pygen;
    
    path_t current_path = fs::current_path();
    path_t root         = current_path.parent_path().parent_path();
    path_t dll_dir;
    path_t template_dir;
    path_t data_dir;
    path_t output_path;
    
    int verbosity;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("output-dir", po::value(&output_path)->default_value(root / "bindings" / "fsharp" / "premia"), "output directory")
        ("dll-dir",    po::value(&dll_dir)->default_value(current_path), "directory containing F# binding library")
        ("data-dir",   po::value(&data_dir)->default_value(root / "data"), "directory containing premia data files")
        ("template-dir", po::value(&template_dir)->default_value(current_path / "templates" / "fsharp"), "directory containing templates for fsgen")
        ("verbosity,v", po::value(&verbosity)->default_value(1), "verbosity level (0 - no output, 1 - basic, 2 - detailed)")
    ;
    
    
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm);   
    
	fsharp::FsCtx ctx(output_path, data_dir, dll_dir / "fspremia.so", template_dir,verbosity);
    
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

    fs::create_directories(output_path);

	InitVar();

	strcpy(premia_data_dir, data_dir.string().c_str());
	
    ctx.out(1) << "Initializing...";

	PremiaCtx p;

    ctx.out(1) << "ok!" << std::endl;
    
	ctx << fsharp::Common() 
	    << p.models << p.families << p.pricings << p.enums
	    << fsharp::FsProj();
   
    return 0;
}
