#pragma once

namespace prxl {

/// \param c - an excel cell
/// \return a cell right next to c
inline Excel::RangePtr Right(Excel::RangePtr c) { return c->GetOffset(0,1); }

/// \param c - an excel cell
/// \return a cell down next to c
inline Excel::RangePtr Down(Excel::RangePtr c) { return c->GetOffset(1,0); }

/// \param c - an excel cell
/// \return a string in A1 format containing address of the cell
inline _bstr_t address_of(Excel::RangePtr C)
{
    return C->GetAddress(VARIANT_TRUE, VARIANT_TRUE, Excel::xlA1);
}

inline _bstr_t full_address_of(Excel::RangePtr C)
{
    return C->Worksheet->Name + "!" + address_of(C);
}

/// global pointer to Excel application
extern CComQIPtr<Excel::_Application>  g_app;

/// \return a pointer to Excel application
inline CComQIPtr<Excel::_Application> & excelApp()
{
    if (g_app == 0)
        throw std::exception("g_app is not initialized");
    return g_app;
}

/// A wrapper over COM objects created directly
/// \param T - COM object class
template <class T>
    struct com_object
    {
        /// Constructor. 
        /// Creates an instance of the COM class
        com_object()
        {
            CComObject<T>::CreateInstance(&p_);
        }

        /// \return a pointer to the underlying COM object
        T * operator -> () 
        {
            return p_;
        }

    private:
        /// the underlying com object
        CComObject <T>   *p_;
    };

}