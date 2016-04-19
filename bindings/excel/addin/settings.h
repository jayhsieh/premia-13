#pragma once

namespace prxl
{



/// entity role
enum Role { 
    /// an entity is a model
    stModel, 
    /// an entity is an option
    stOption, 
    /// an entity is a pricing method
    stMethod, 
    /// an entity is a problem result 
    stResult, 
    /// the mark of the end of the enumeration
    stEnd 
};

/// \fn role_of
/// overloaded function returning role of the entity

inline Role role_of(Model *) { return stModel; }
inline Role role_of(Option *) { return stOption; }
inline Role role_of(PricingMethod *){ return stMethod; }

template <class T>
Role role_of()
{
    return role_of((T*)0);
}

namespace env
{
    /// return string to be put as label to top-left corner of an entity region
    /// \param r - the entity role
    /// \return a string with role of the entity
    inline std::string getLabel(Role r)
    {
        return 
            r == stModel ? "Model" :
            r == stOption ? "Option" : 
            r == stMethod ? "Method" : 
            r == stResult ? "Result" :
            "-=Unknown=-";
    }

    /// returns db sheet name for the role 'r'
    /// \param r - sheet role
    /// \return  db sheet name
    inline std::string getSheetName(Role r)
    {
        return std::string("premia_11_") + getLabel(r);
    }

    /// returns Excel color for header cell background for role 'r'
    /// \param r - entity role
    /// \return Excel color for header cell background 
    inline long getHeaderBkColor(Role r)
    {
        return 
            r == stModel ? 4 :
            r == stOption ? 6 :
            r == stMethod ? 37 :
            r == stResult ? 3 :
            -1;
    }

    /// returns Excel color for header font for role 'r'
    /// \param r - entity role
    /// \return Excel color for header font
    inline long getHeaderFontColor(Role r)
    {
        return 
            r == stResult ? 2 : 1;
    }

    /// returns Excel color for cell background for role 'r'
    /// \param r - entity role
    /// \return Excel color for cell background
    inline long getCellBkColor(Role r)
    {
        return 
            r == stResult ? 35 : 15;
    }

    /// returns Excel color for cell font for role 'r'
    /// \param r - entity role
    /// \return Excel color for cell font 
    inline long getCellFontColor(Role r)
    {
        return 
            r == stResult ? 55 : 1;
    }
}}