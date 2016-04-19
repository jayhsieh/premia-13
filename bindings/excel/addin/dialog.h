#pragma once

#include "premia.h"

#undef BOOL
#define BOOL int

#include "excel.h"
#include "db.h"
#include "readingcells.h"

#include "history.h"
#include "mode_panel.h"

namespace prxl {
namespace wzd {

/// Clears content of WTL listbox or combobox
/// \param ctrl - a control to be cleared
template <class Control>
    void ClearContent(Control & ctrl)
{
    while (ctrl.GetCount() > 0)
        ctrl.DeleteString(0);
}

/// Add a string to the (list,combo)-box and associates with it some pointer.
/// This pointer can be retrieved later by means of GetItemDataPtr function.
/// \param ctrl - a control where to add the string
/// \param str - a string to be added
/// \param ptr - data to be associated
template <class Control>
    void AddStringWithPtr(Control & ctrl, const char * str, void * ptr)
{
    ctrl.SetItemDataPtr(ctrl.AddString(CComBSTR(str)), ptr);
}
}}

#include "existing_models_ctrl.h"
#include "create_new_model_ctrl.h"
#include "existing_create_new_panel.h"
#include "create_new_method_ctrl.h"
#include "existing_methods_ctrl.h"
#include "create_new_option_ctrl.h"
#include "existing_options_ctrl.h"

namespace prxl {
/// Namespace for Premia Excel Wizard
namespace wzd {
/// Premia wizard dialog.
class CWizard : public CSimpleDialog<IDD_DIALOG1>
{
private:
    /// panel for choosing mode
    ModePanel   mode_;

    template <class T> struct Panel;
    //---------------------------------------------------------- Models

    template <>
        struct Panel<Model>
        {
            typedef 
                ExistingCreateNewPanel<
                    ExistingModelsControl<CWizard>, 
                    CreateNewModelControl<CWizard>, 
                    IDC_RB_EXISTING_MODEL, 
                    IDC_RB_NEW_MODEL 
                > 
                type;
        };

    /// panel for choosing model
    Panel<Model>::type   model_;

    //---------------------------------------------------------- Options
    template <>
        struct Panel<Option>
        {
            typedef 
                ExistingCreateNewPanel<
                    ExistingOptionsControl<CWizard>, 
                    CreateNewOptionControl<CWizard>, 
                    IDC_RB_OPTIONS_EXISTING, 
                    IDC_RB_OPTION_NEW 
                > 
                type;
        };

    /// panel for choosing option
    Panel<Option>::type   option_;

    //---------------------------------------------------------- Methods
    template <>
        struct Panel<PricingMethod>
        {
            typedef 
                ExistingCreateNewPanel<
                    ExistingMethodsControl<CWizard>,
                    CreateNewMethodControl<CWizard>,
                    IDC_RB_METHOD_EXISTING,
                    IDC_RB_METHOD_NEW
                >
                type;
        };

    /// panel for choosing method
    Panel<PricingMethod>::type      methods_;

public:

    /// \return current mode
    lib::Mode    getMode() const { return mode_.getValue(); }

    // \return current family
    Family* getFamily()
    {
        return wzd::getFamily<Family*>(option_);
    }

    template <class T>
        T * get();

    // \return current option 
        template <> Option * get<Option>() { return getValue<Option*>(option_); }
        template <> Model  * get<Model>() { return getValue<Model*>(model_); }
    // \return current pricing method
        template <> PricingMethod * get<PricingMethod>() { return getValue<PricingMethod*>(methods_); }

    template <class T>
        typename Panel<T>::type  & getPanel();

        template <> Panel<Model>::type & getPanel<Model>() { return model_; }
        template <> Panel<Option>::type & getPanel<Option>() { return option_; }

        // \return reference to the panel for choosing a pricing method
        template <> Panel<PricingMethod>::type & getPanel<PricingMethod>() { return methods_; }

    // \return current pricing
    Pricing * getPricing() const 
    { 
        return wzd::getPricing<Pricing*>(option_);
    }

 public:
    enum { IDD = IDD_DIALOG1 };

    BEGIN_MSG_MAP(CWizard)
        MESSAGE_HANDLER(WM_INITDIALOG, OnInitDialog)
        CHAIN_MSG_MAP_MEMBER(mode_)
        CHAIN_MSG_MAP_MEMBER(model_)
        CHAIN_MSG_MAP_MEMBER(option_)
        CHAIN_MSG_MAP_MEMBER(methods_)
        CHAIN_MSG_MAP(CSimpleDialog<IDD_DIALOG1>)
    END_MSG_MAP()

    /// A constructor.
    CWizard()
        :   model_(*this)
        ,   option_(*this)
        ,   methods_(*this)
    {
    }

    /// Calls Init function of the first control and ties the first control with update function of the second control
    /// \param what - a control to be initialized and tied with update member of 'with'
    /// \param with - a control 'what' should be tied to
    template <class What, class With>
        void init_and_tie(What & what, With & with)
        {
            what.Init(*this);
            what.Changed = boost::bind(&With::update, boost::ref(with));
        }

    /// Initializes dialog controls.
    /// Called on dialog creation.
    LRESULT OnInitDialog(UINT, WPARAM, LPARAM, BOOL&)
    {
        init_and_tie(mode_, model_);
        init_and_tie(model_, option_);
        init_and_tie(option_, methods_);

        methods_.Init(*this);

        mode_.update();

        return 0;
    }



    LRESULT OnClose(WORD, WORD, HWND, BOOL&)
    {
        return 0;
    }
};
}}