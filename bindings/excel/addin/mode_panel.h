#pragma once

namespace prxl {
namespace wzd {

/// A panel with 3 radio buttons for choosing a premia mode
struct ModePanel
{
    /// A constructor
	ModePanel() : mode_(0) {}

    /// Binds controls with their resource ids
    /// \param parent - window containing radio buttons with ids = {IDC_EQUITY,IDC_INTEREST_RATE,IDC_ENERGY}
    void Init(CDialogImplBase & parent)
    {
        m_combo = parent.GetDlgItem(IDC_MODE);

        int id = 0;

		for (int i = 0; premia_assets[i].name != 0; ++i)
        {
            m_combo.SetItemData(m_combo.AddString(_bstr_t(premia_assets[i].name)), i);
            
            if (i == mode_)
                id = i;
        }

        m_combo.SetCurSel(id);
    }

    BEGIN_MSG_MAP(ModePanel)
        COMMAND_HANDLER(IDC_MODE, CBN_SELCHANGE, OnClicked)
    END_MSG_MAP()

    /// Click event handler
    LRESULT OnClicked(WORD d1, WORD d2, HWND d3, BOOL& bHandled)
    {
        update();

        bHandled = TRUE;

        return 0;
    }

    /// Updates current mode in respect to state of the radio buttons
    void update()
    {
        mode_ = lib::Mode(m_combo.GetItemData(m_combo.GetCurSel()));

        Changed();
    }

    /// \return current mode
    lib::Mode getValue() const { return mode_; }

    /// A slot for mode change event
    boost::function<void ()>    Changed;

private:
    WTL::CComboBox  m_combo;
    /// current mode
    lib::Mode   mode_;
};
}}