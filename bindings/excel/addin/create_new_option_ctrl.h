#pragma once

#include "premia.h"
#include "history.h"

namespace prxl {
namespace wzd {
/// A wrapper over a listbox control and a combo-box control 
/// that adds a logic for choosing type of a new option instance.
/// The combobox with ID = IDC_FAMILY is used to choose between option families
/// The listbox with ID = IDC_OPTION is used choose an option type
template <class Parent>
struct CreateNewOptionControl
{
    /// Constructor.
    /// \param p - parent control which can be queried for current model
    CreateNewOptionControl(Parent & p) : parent_(p) 
    {}

    /// Binds the control to its resource ID
    /// \param wnd - a window that has a list-box with ID = IDC_OPTION and a combobox with ID = IDC_FAMILY
    void Init(CDialogImplBase & wnd)
    {
        m_families = parent_.GetDlgItem(IDC_FAMILY);
        m_options = parent_.GetDlgItem(IDC_OPTION);

		PremiaAsset * pa;
		for (pa = premia_assets; pa->name != NULL; ++pa);
		family_.resize(pa - premia_assets);
    }

    /// \return a pointer to currently selected option
    Option * getValue() const
    {
        return getOption();
    }

    /// updates the control content respecting current model. 
    void update()
    {
        ClearContent(m_families);

        // index of combo-box element to be selected
        int idx = 0;

        // current mode
        lib::Mode mode = parent_.getMode();

        HistoryChooser<Family*>  chooser(family_[mode]);

        // number of families have been added to the combo-box
        int added = 0;

        // ID of the currently selected model
        CComBSTR model_id = parent_.get<Model>()->ID;

        // for each families of the mode
        for (Family ** f_it = lib::getFamilies(mode); *f_it; ++f_it)
        {
            // family name is ID of any option of the family
            char* family_name = (**f_it)[0]->ID;

            // if the family is compatible with the current model
            if (lib::getPricing(mode, parent_.get<Model>(), (**f_it)[0]))
            {
                // add family to the combo-box
                AddStringWithPtr(m_families, family_name, *f_it);

                // and if it was been selected earlier, select it again
                if (chooser.contains(*f_it))
                    idx = added;

                ++added;
            }
        }

        // triggers option list update
        m_families.SetCurSel(idx);
    }

    /// Updates options list-box respecting currently selected family and model
    void fillOptionsList()
    {
        ClearContent(m_options);

        // current model
        Model * m = parent_.get<Model>();
        // current family
        Family * f = getFamily();

        m->Init(m);

        // pricing for the model and the family
        // it exists due to the way we fill the families combo-box
        Pricing * p = lib::getPricing(parent_.getMode(), m, (*f)[0]);

        HistoryChooser<Option*>  chooser(option_[f]);

        // added - number of elements added to the list-box
        // idx - element to be selected
        int added = 0, idx = 0;

        // for option in the family
        for (Option** o_it = *f; *o_it; ++o_it)
        {
            Option * opt = *o_it;

            opt->Init(opt, m);

            // look through all pricing methods in the corresponding pricing
            for (PricingMethod ** pm_it = p->Methods; *pm_it; ++pm_it)
            {
                // if some method is compatible with the option
                if ((*pm_it)->CheckOpt(opt, m) == OK)
                {
                    // add the option to combo-box
                    AddStringWithPtr(m_options, opt->Name, opt);

                    // check whether it was selected before
                    if (chooser.contains(opt))
                        idx = added;

                    ++added;

                    break;
                }
            }
        }

        m_options.SetCurSel(idx);
    }


    /// Family selection change event handler.
    /// Updates family selection history and the options list.
    LRESULT OnClicked(WORD d1, WORD d2, HWND d3, BOOL& d4)
    {
        return OnFamily(d1,d2,d3,d4);
    }

    /// Family selection change event handler.
    /// Updates family selection history and the options list.
    LRESULT OnFamily(WORD d1, WORD d2, HWND d3, BOOL& d4)
    {
        // selected family
        Family * f = reinterpret_cast<Family*>(m_families.GetItemDataPtr(m_families.GetCurSel()));

        // current mode
        lib::Mode mode = parent_.getMode();

        // add the family into the history
        family_[mode].push(f);

        // update current pricing
        pricing_ = lib::getPricing(mode, parent_.get<Model>(), (*f)[0]);

        fillOptionsList();

        return OnOption(d1,d2,d3,d4);
    }

    /// Option selection change event handler.
    /// Updates option selection history
    LRESULT OnOption(WORD d1, WORD d2, HWND d3, BOOL& d4)
    {
        option_[getFamily()].push(reinterpret_cast<Option*>(m_options.GetItemDataPtr(m_options.GetCurSel())));

        Changed();

        return 0;
    }

    /// A slot for option selection change event.
    boost::function<void ()>    Changed;

    BEGIN_MSG_MAP(CreateNewOptionControl)
        COMMAND_HANDLER(IDC_FAMILY, CBN_SELCHANGE, OnFamily)
        COMMAND_HANDLER(IDC_FAMILY, CBN_SETFOCUS, OnFamily)
        COMMAND_HANDLER(IDC_OPTION, LBN_SELCHANGE, OnOption)
    END_MSG_MAP()

    /// \return current family
    Family* getFamily() const
    {
        return family_[parent_.getMode()].top(); 
    }

    /// \return current option
    Option* getOption() const
    { 
        return option_.find(getFamily())->second.top(); 
    }

    /// \return current pricing
    Pricing*  getPricing() const { return pricing_; }

private:
    /// a reference to a parent control
    Parent      &   parent_;
    /// ATL combo-box for families
    CComboBox       m_families;
    /// ATL list-box for options
    CListBox        m_options;

    // family selection history: Mode --> Family
	std::vector<History<Family *> >         family_;
    // option selection history: Family --> history of options
    std::map<Family*, History<Option *> >   option_;
    // current pricing
    Pricing *                               pricing_;
};
}}