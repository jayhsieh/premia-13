#pragma once

namespace prxl {
namespace wzd {
/// Mode how to use an entity for a problem
enum CreationMode { 
    /// use existing entity instance
    UseExisting, 
    /// create a new entity instance
    CreateNew 
};

/// retrieves the principal value (e.g. pointer to current model for ModelPanel) from a composite control
/// \param T - type of the value
/// \param ExistingNew - type of the composite control
/// \param x - the composite control
/// \return the principal value
template <class T, class ExistingNew>
T  getValue(ExistingNew const & x)
{
    return 
        x.getCreationMode() == UseExisting
        ?   x.existing().getValue() 
        :   x.create_new().getValue();
}

/// retrieves current family from OptionsPanel
/// \param x -- OptionPanel
template <class T, class ExistingNew>
T  getFamily(ExistingNew const & x)
{
    return 
        x.getCreationMode() == UseExisting
        ?   x.existing().getFamily() 
        :   x.create_new().getFamily();
}

/// retrieves current pricing from OptionsPanel
/// \param x - OptionPanel
template <class T, class ExistingNew>
T  getPricing(ExistingNew const & x)
{
    return 
        x.getCreationMode() == UseExisting
        ?   x.existing().getPricing() 
        :   x.create_new().getPricing();
}

/// Composite control consisting of two controls: one for choosing existing entity instance and
/// another one - for choosing a type of an entity to be created
/// \param ExistingCtrl - type of a control for choosing existing entity instance
/// \param CreateNewCtrl - type of a control for choosing a type of an entity to be created
/// \param ID_EXISTING - resource id of a radio button for switching to Existing creation mode 
/// \param ID_CREATE_NEW - resource id of a radio button for switching to CreateNew creation mode
template <class ExistingCtrl, class CreateNewCtrl, int ID_EXISTING, int ID_CREATE_NEW>
struct ExistingCreateNewPanel
{
    /// A constructor
    template <class Parent>
    ExistingCreateNewPanel(Parent & p) 
        :   existing_(p)
        ,   create_new_(p)
        ,   creation_mode_(CreateNew)
    {
        existing_.Changed = boost::bind(&ExistingCreateNewPanel::OnChanged, this, UseExisting);
        create_new_.Changed = boost::bind(&ExistingCreateNewPanel::OnChanged, this, CreateNew);

        existing_.Enabled = boost::bind(&ExistingCreateNewPanel::OnEnableExisting, this, _1);
    }

    BEGIN_MSG_MAP(ExistingCreateNewPanel)
        COMMAND_HANDLER(ID_EXISTING,   BN_CLICKED, OnMode)
        COMMAND_HANDLER(ID_CREATE_NEW, BN_CLICKED, OnMode)
        CHAIN_MSG_MAP_MEMBER(create_new_)
        CHAIN_MSG_MAP_MEMBER(existing_)
    END_MSG_MAP()

    /// Binds controls with its resource ids
    /// \param wnd - a window that can queried for resource handles
    void Init(CDialogImplBase & wnd)
    {
        existing_.Init(wnd);
        create_new_.Init(wnd);
        m_rb_existing = wnd.GetDlgItem(ID_EXISTING);
        m_rb_create_new = wnd.GetDlgItem(ID_CREATE_NEW);
    }

    /// \return current creation mode
    CreationMode  getCreationMode() const 
    {
        return creation_mode_;
    }

    /// \return reference to a panel for choosing existing entity instances
    ExistingCtrl const & existing() const { return existing_; }

    /// \return reference to a panel for choosing a type for an entity to be created
    CreateNewCtrl const & create_new() const { return create_new_; }

    /// Handles radio buttons clicks.
    LRESULT OnMode(WORD d1, WORD d2, HWND d3, BOOL& d4)
    {
        creation_mode_ = m_rb_create_new.GetCheck() == BST_CHECKED ? CreateNew : UseExisting;

        return 
            creation_mode_ == UseExisting
            ?   existing_.OnClicked(d1,d2,d3,d4)
            :   create_new_.OnClicked(d1,d2,d3,d4);
    }

    /// Converts a boolean value into radio button state
    static int toBST(bool b) { return b ? BST_CHECKED : BST_UNCHECKED; }

    /// Updates child windows and sets CreateNew mode
    void update()
    {
        create_new_.update();
        existing_.update();

        setMode(CreateNew);

        BOOL b;
        OnMode(0,0,0,b);
    }

    /// A slot for handling an event about current entity has changed
    boost::function<void ()>    Changed;

private:

    /// Enables or disable window for choosing an existing entity instance
    /// \param bEnable - if true then enable otherwise disable the window
    void OnEnableExisting(bool bEnable)
    {
        m_rb_existing.EnableWindow(bEnable);
    }

    /// handles events about selection changes in child windows
    /// \param m -- identifies the source of the event
    void OnChanged(CreationMode m)
    {
        setMode(m);

        if (Changed)
            Changed();
    }

    /// Sets the current creation mode
    /// \param m - creation mode 
    void setMode( CreationMode m )
    {
        creation_mode_ = m;

        m_rb_existing.SetCheck(toBST(m == UseExisting));
        m_rb_create_new.SetCheck(toBST(m == CreateNew));
    }
private:
    /// radio button for switching into Existing creation mode
    CButton         m_rb_existing;
    /// radio button for switching into CreateNew creation mode
    CButton         m_rb_create_new;
    /// current creation mode
    CreationMode    creation_mode_;

    /// control for choosing existing entity instance
    ExistingCtrl    existing_;
    /// control for choosing a type of an entity to be created
    CreateNewCtrl   create_new_;
};
}}