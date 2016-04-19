#pragma once

namespace prxl {
namespace fmt {
/// Handling zero coupon bond prices region
namespace zcb {

/// Tests is the string a key for ZCB Prices region parameter.
/// If so, trims the string
/// \param pKey - a pointer to a string being processed
/// \return true iff the string is a key for ZCB Prices region parameter.
inline bool isZcbParameter(std::string *pKey)
{
    return false;

    std::string::size_type zcb = pKey->find("ZCB Prices in data/initialyield.dat");

    if (zcb != std::string::npos)
    {
        *pKey = pKey->substr(0, zcb - 2);
        return true;
    }

    return false;
}

/// Create a region with ZCB prices.
/// The region will be filled with data from file initialyield.dat which should be in the current directory.
/// \param reference_to_zcb_region - a cell where the reference to the zcb region will be put to
/// \param topleft - the top-left corner for zcb region
/// \return a reference to the cell down to the region 
inline Excel::RangePtr createZcbRegion(Excel::RangePtr reference_to_zcb_region, Excel::RangePtr topleft)
{
    Excel::RangePtr C = topleft;

    printInBoldCell("Term structure", C, 34);
    printInBoldCell("(NO_FLAT case)", Right(C), 34);

    C = Down(C);

    printInBoldCell("Maturities", C, 46);
    printInBoldCell("Market prices", Right(C), 46);

    C = Down(C);

    putRef(reference_to_zcb_region, C);

    std::ifstream in("initialyield.dat");

    while (in)
    {
        double val;
        std::string stime;

        in >> val >> stime;

        if (stime[0] == 't' && stime[1] == '=')
        {
            double time = atof(stime.c_str() + 2);

            C->Value2 = time;
            Right(C)->Value2 = val;

            C = Down(C);
        }
    }

    return C;
}

/// Read data from a zcb pricies region and puts them to initialyield.dat file.
/// \param C - the top-left corner of zcb prices region
inline void readZcbRegion(Excel::RangePtr C)
{
    std::ofstream out("initialyield.dat");

    while (C->Value2 != _variant_t())
    {
        out << (double)Right(C)->Value2 << " t=" << (double)C->Value2 << std::endl;

        C = Down(C);
    }
}

/// Detects whether the string is a key for zcb prices parameter.
/// \param key - a string to be tested
/// \return true iff the string is key for zcb prices parameter
inline bool isZcbParameterInCell(_bstr_t const & key)
{
    return key == _bstr_t("ZCB Prices Region");
}

}}}