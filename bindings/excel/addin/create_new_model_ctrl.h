#pragma once

#include "premia.h"

namespace prxl {
namespace wzd {
/// A wrapper over listbox control that adds a logic for choosing type of a new model instance.
template <class Parent>
struct CreateNewModelControl
{
    /// Constructor.
    /// \param p - parent control which can be queried for current mode
    CreateNewModelControl(Parent & p) : parent_(p) 
    {
		for (PremiaAsset* a = premia_assets; a->name != NULL; ++a)
		{
			history_.push_back(0);
		}
    }

    /// Binds the control to its resource ID
    /// \param wnd - a window that has a list-box with ID = IDC_MODELS 
    void Init(CDialogImplBase & wnd)
    {
        m_models = wnd.GetDlgItem(IDC_MODELS);
    }

    /// \return a pointer to currently selected model
    Model * getValue() const
    {
        return history_[parent_.getMode()];
    }

    /// updates the control content respecting current mode. 
    void update()
    {
        ClearContent(m_models);

        // index of a list-box element to be selected
        int idx = 0;

        // for each model in current mode
        for (Model ** m_it = lib::getModels(parent_.getMode()); *m_it; ++m_it)
        {
            // add it into the list
            AddStringWithPtr(m_models, (*m_it)->Name, *m_it);

            // if it was selected earlier restore selection on it
            if (*m_it == getValue())
                idx = int(m_it - lib::getModels(parent_.getMode()));
        }

        m_models.SetCurSel(idx);
    }

    /// Selection change event handler.
    /// Updates selection history.
    LRESULT OnClicked(WORD d1, WORD d2, HWND d3, BOOL& d4)
    {
        history_[parent_.getMode()] = reinterpret_cast<Model*>(m_models.GetItemDataPtr(m_models.GetCurSel()));

        Changed();

        return 0;
    }

    /// A slot for selection change event.
    boost::function<void ()>    Changed;

    BEGIN_MSG_MAP(CreateNewModelControl)
        COMMAND_HANDLER(IDC_MODELS, LBN_SELCHANGE, OnClicked)
    END_MSG_MAP()

private:
    /// a reference to a parent control
    Parent      &   parent_;
    /// ATL list-box
    CListBox        m_models;
    /// selection history: Mode --> the last selected model
	std::vector<Model*>  history_;
};
}}