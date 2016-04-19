#pragma once

namespace prxl {
namespace wzd  {

/// A wrapper over combo-box adding a logic for choosing an existing entity instance.
/// \param ID - resource ID of the underlying combobox
template <int ID>
struct ComboboxOfExisting
{
    BEGIN_MSG_MAP(ComboboxOfExisting)
        COMMAND_HANDLER(ID, CBN_SELCHANGE, OnClicked)
        COMMAND_HANDLER(ID, CBN_SETFOCUS, OnClicked)
    END_MSG_MAP()

    /// Binds the control to its resource ID
    /// \param wnd - a window that has a combo-box with ID = this->ID
    void Init(CDialogImplBase & wnd)
    {
        m_existing = wnd.GetDlgItem(ID);
    }

    /// Updates combobox contents
    void update()
    {
        ClearContent(m_existing);

        // ask a derived class to fill the combobox
        fillList();

        int added = m_existing.GetCount();

        // if nothing has been added
        if (added == 0)
        {
            // disable the window
            m_existing.EnableWindow(false);
        }
        else    // otherwise
        {
            // enable it
            m_existing.EnableWindow(true);
            m_existing.SetCurSel(0);
        }

        // trigger the event only if there are some elements
        Enabled(added > 0);
    }

    /// Selection change event handler.
    LRESULT OnClicked(WORD d1, WORD d2, HWND d3, BOOL& d4)
    {
        OnChanged(m_existing.GetItemData(m_existing.GetCurSel()));

        Changed();

        return 0;
    }

    /// A slot for selection change event
    boost::function<void ()>        Changed;
    /// A slot for an event about the control enabling/disabling 
    boost::function<void (bool)>    Enabled;

protected:

    /// A hook for derived classes letting them to fill the list
    virtual void fillList() = 0;

    /// Adds an element to the list and associates some data with it
    /// \param label - element label
    /// \param data - data to be associated
    void AddString(_bstr_t label, int data)
    {
        m_existing.SetItemData(m_existing.AddString(label), data);
    }

    /// A hook for derived classes.
    /// Called when combo-box selection has changed
    /// \param data - data saved in selected element (id of existing entity)
    virtual void OnChanged(int data) = 0;

private:
    /// the underlying combo-box
    CComboBox   m_existing;
};
}}