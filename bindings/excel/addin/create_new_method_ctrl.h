#pragma once

#include "premia.h"
#include "history.h"

namespace prxl {
namespace wzd {

/// A wrapper over listbox control that adds a logic for choosing type of a new pricing method instance.
template <class Parent>
    struct CreateNewMethodControl
{
    /// Constructor.
    /// \param p - parent control which can be queried for current option, model and so on
    CreateNewMethodControl(Parent & p) : parent_(p) 
    {
    }

    /// Binds the control to its resource ID
    /// \param wnd - a window that has a list-box with ID = IDC_METHODS 
    void Init(CDialogImplBase & wnd)
    {
        m_methods = wnd.GetDlgItem(IDC_METHODS);
    }

    /// \return a pointer to currently selected pricing method
    PricingMethod * getValue() const
    {
        return method_.find(parent_.getPricing())->second.top();;
    }

    /// updates the control content respecting currently selected model and option. 
    void update()
    {
        ClearContent(m_methods);

        // currently selected pricing
        Pricing * p = parent_.getPricing();

        HistoryChooser<PricingMethod*> chooser(method_[p]);

        // added - how many items have been added to the list-box
        // idx - index of the list-box element to be selected
        int added = 0, idx = 0;

        // current model
        Model * m = parent_.get<Model>();
        // current option
        Option * opt = parent_.get<Option>();

        // for each method in the pricing
        for (PricingMethod ** pm_it = p->Methods; *pm_it; ++pm_it)
        {
            (*pm_it)->Init(*pm_it, opt);

            // if it is compatible with option and model
            if ((*pm_it)->CheckOpt(opt, m) == OK)
            {
                // ... add it to the list-box
                AddStringWithPtr(m_methods, (*pm_it)->Name, *pm_it);

                // and if it was selected before 
                if (chooser.contains(*pm_it))
                    // select it again
                    idx = added;

                ++added;
            }
        }

        m_methods.SetCurSel(idx);
    }

    /// Selection change event handler.
    /// Updates selection history.
    LRESULT OnClicked(WORD d1, WORD d2, HWND d3, BOOL& d4)
    {
        // put selected pricing method into the history
        method_[parent_.getPricing()].push(
            reinterpret_cast<PricingMethod*>(m_methods.GetItemDataPtr(m_methods.GetCurSel())));

        Changed();

        return 0;
    }

    /// A slot for selection change event.
    boost::function<void ()>    Changed;

    BEGIN_MSG_MAP(CreateNewMethodControl)
        COMMAND_HANDLER(IDC_METHODS, LBN_SELCHANGE, OnClicked)
    END_MSG_MAP()

private:
    /// a reference to a parent control
    Parent      &   parent_;
    /// ATL list-box
    CListBox        m_methods;
    /// selection history: Pricing --> History of selected pricing methods
    std::map<Pricing*, History<PricingMethod *> > method_;
};
}}