#ifdef _MSC_VER
#  define BOOST_FILESYSTEM_VERSION 3
#endif

#include <fstream>
#include <iostream>

#include <boost/filesystem.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options.hpp>

#include <premia/generator/api/all.h>
#include <premia/generator/utils/fspath.h>

namespace fs = boost::filesystem;

#include <set>

namespace premia {

    template <fs::path (*Func)(fs::path const&)>
        struct FullPath
        {
            FullPath(fs::path const &base) : base(base) {}

            template <class T> fs::path operator () (T const &x) const 
            {
                return base / Func(relativeDocPath(x));
            }

            bool to_check() const 
            {
                return !base.empty();
            }

        private:
            fs::path const base;
        };

    template <class T, class PathFunc>
        void check(T const &x, PathFunc const &fp, std::ostream &out)
        {
            fs::path f = fp(x);

            if (!fs::exists(f))
            {
                out << f << " doesn't exist" << std::endl;
            }
        }

    struct PdfTeXCheck
    {
        PdfTeXCheck(fs::path const &pdf_base, fs::path const &tex_base, std::ostream &out)
            :   pdf(pdf_base), tex(tex_base), out(out)
        {}

        template <class T>
            void operator () (T const &x) const 
            {
                if (pdf.to_check())
                    check(x, pdf, out);

                if (tex.to_check())
                    check(x, tex, out);
            }

    private:
        FullPath<pygen::api::pdf> const pdf;
        FullPath<pygen::api::tex> const tex;
        std::ostream &out;
    };

    void checkHelpFileNames(fs::path const &pdf_base, fs::path const &tex_base, std::ostream &out)
    {
        using pygen::api::pdf;

        InitVar();
        pygen::PremiaCtx ctx;

        PdfTeXCheck check(pdf_base, tex_base, out);

        BOOST_FOREACH(pygen::Model const &m, ctx.models.models)
        {
            check(m);
        }

        BOOST_FOREACH(pygen::Family const &f, ctx.families.families) {
            std::for_each(f.options.begin(), f.options.end(), check);
        }

        BOOST_FOREACH(pygen::Pricing const &p, ctx.pricings.pricings) {
            std::for_each(p.methods.begin(), p.methods.end(), check);
        }
    }

}

int main(int argc, char* argv[])
{
    namespace po = boost::program_options;

    FsPath pdf_base;
    FsPath tex_base;

    po::options_description desc("Allowed options");
    desc.add_options()
        ("help", "produce help message")
        ("pdf-dir", po::value(&pdf_base), "directory for PDF documentation")
        ("tex-dir", po::value(&tex_base), "directory for TeX documentation")
        ;

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, desc), vm);
    po::notify(vm); 

    if (vm.count("help") || vm.size() == 0) 
    {
        std::cout << "Checks that for every Premia model, option and pricing method there is a corresponding documentation file.\n"
            "It may look for PDF or/and TeX documentation." << std::endl;
        std::cout << desc << "\n";
        return 1;
    }    

    if (pdf_base->empty() && tex_base->empty())
    {
        std::cerr << "Please specify either --pdf-dir or --tex-dir" << std::endl;
    }
    else 
    {
        premia::checkHelpFileNames(*pdf_base, *tex_base, std::cout);
    }

	return 0;
}

