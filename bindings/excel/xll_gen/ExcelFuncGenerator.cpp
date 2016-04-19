// ExcelFuncGenerator.cpp : définit le point d'entrée pour l'application console.
//

#include "stdafx.h"
#include "..\addin\premia.h"

namespace prxl {

    /// generates a header for file generated.cpp
void generateHeader(std::ostream & out)
{
    out <<
"#include <xlw/pragmas.h>\n"
"#include <xlw/MyContainers.h>\n"
"#include <xlw/CellMatrix.h>\n"
"#include <xlw/xlw.h>\n"
"#include <xlw/XlFunctionRegistration.h>\n"
"#include <stdexcept>\n"
"#include <xlw/XlOpenClose.h>\n"
"#include <ctime>\n"
"#include <xlw/ArgList.h>\n"

"#include <xlw/ArgListFactory.h>\n"

"#include <xlw/xlarray.h>\n"

"        namespace {\n"
"            const char* LibraryName = \"MyTestLibrary\";\n"
"    };\n"

"    namespace {\n"
"        void DummyFunction()\n"
"        {\n"
"            xlAutoOpen();\n"
"            xlAutoClose();\n"
"        }\n"
"    }\n";
}

/// to be replaced by boost::lexical_cast
std::string toString(int i)
{
    std::ostringstream s;
    s << i;
    return s.str();
}

/// Generates a source for an entity function definition
/// \param out - output stream
/// \param name - name of the entity
/// \param vars - an array of input parameters of the entity
/// \param nvar - number of input parameters of the entity
/// \param role - a string describing role of the entity
void generateFunc(std::ostream & out, const char * name, VAR const * vars, int nvar, const char * role)
{
    out << "namespace \n{\n  XLRegistration::Arg\n\t\t";
    out << name << "Args[]=\n{\n";

    // printing entity parameters
    for (int i = 0; i != nvar; ++i)
    {
        VAR const & v = vars[i];

        int vtype = true_typeV[v.Vtype];

        out << "      { \"" << v.Vname << "\", \" default value = " << (vtype == 2 ? v.Val.V_DOUBLE : v.Val.V_INT) << " \"},\n";
    }

    out << "    };\n    XLRegistration::XLFunctionRegistrationHelper register" << name << 
        "(\"xl" << name << "\", \n\t\"" << name << "\", \"  \", \n\t\"Premia" << role << "\", \n\t" 
        << name << "Args, \""
        << std::string(nvar, 'R') << "\");\n}\n";

    out << "extern \"C\"\n{        LPXLOPER EXCEL_EXPORT xl" << name << "(\n";

    for (int i = 0; i != nvar; ++i)
    {
        VAR const & v = vars[i];

        out << "      LPXLOPER p_" << i << "a" << (i == nvar - 1 ? ")" : ",") << "\n";
    }

    out << "\t{\n\t\tEXCEL_BEGIN;\n";

    out << "\t\tstd::string result = \"" << name << "(\";\n";

    for (int i = 0; i != nvar; ++i)
    {
        VAR const & v = vars[i];

        out << "\t\tresult += XlfOper(p_" << i << "a).AsString(\"\")";

        if (i != nvar - 1)
            out << "+std::string(\";\")";

        out << ";\n";
    }

    out << "\t\tresult += \")\";\n";

    out << "\t\treturn XlfOper(CellMatrix(result));\n";

    out << "\t\tEXCEL_END;\n\t}\n}\n";

}

/// generates 
void generateFunc(std::ostream & out, const char * name, std::vector<VAR> pars, const char * role)
{
    lib::removeCalculableFields(pars);

    if (pars.size() == 0)
    {
        VAR v;
        v.Vname = "dummy";
        v.Vtype = INT;
        v.Val.V_INT = 0;
        pars.push_back(v);
    }

    generateFunc(out, name, &pars[0], (int)pars.size(), role);
}

void generateFunc(std::ostream & out, const char * name, std::set<VAR, lib::EqualNames> const & pars, const char * role)
{
    std::string n = name;

    std::replace(n.begin(),n.end(), ' ', '_');

    generateFunc(out, n.c_str(), std::vector<VAR>(pars.begin(), pars.end()), role);    
}

}
int _tmain(int argc, _TCHAR* argv[])
{
    using namespace prxl;
    InitVar();

    std::ofstream out(argc == 1 ? L"generated.cpp" : argv[1]);

    generateHeader(out);

    for (Model ** m_it = models_e; *m_it; ++m_it)
    {
        Model * mod = *m_it;
        (mod->Init)(mod);
        generateFunc(out, mod->Name, lib::getParameters(mod), "Model");
    }

    for (Family ** f_it = families_e; *f_it; ++f_it)
    {
        for (Option ** o_it = **f_it; *o_it; ++o_it)
        {
            Option * opt = *o_it;
            (opt->Init)(opt, &BSND_model);
            generateFunc(out, opt->Name, lib::getParameters(opt), (std::string("Options_") + opt->ID).c_str());
        }
    }

    std::set<std::string>   methods_processed;

    for (Pricing ** pr_it = pricings_e; *pr_it; ++pr_it)
    {
        for (PricingMethod ** m_it = (*pr_it)->Methods; *m_it; ++m_it)
        {
            PricingMethod * met = *m_it;

            if (methods_processed.find(met->Name) == methods_processed.end())
            {
                methods_processed.insert(met->Name);

                generateFunc(out, met->Name, lib::methodParameters()[met->Name], (std::string("Methods_") + (*pr_it)->ID).c_str());
            }
        }
    }

	return 0;
}

