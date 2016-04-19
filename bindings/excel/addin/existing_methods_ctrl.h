#pragma once

#include "result_entry.h"
#include "existing_combo.h"

namespace prxl {
namespace wzd {
/// Defines whether a pricing contains a pricing method
/// \param pricing - a pricing
/// \param method - a pricing method
/// \return true iff the pricing contains the method
inline bool contains(Pricing * pricing, PricingMethod * method)
{
    for (PricingMethod ** pm_it = pricing->Methods; *pm_it; ++pm_it)
    {
        if (*pm_it == method)
            return true;
    }

    return false;
}

/// A wrapper over a list-box adding a logic for selecting an existing pricing method instance.
/// \param Parent - a class that can be queried for a current model, option etc.
template <class Parent>
struct ExistingMethodsControl 
    :   ComboboxOfExisting<IDC_METHODS_EXISTING>
{
    /// A constructor.
    /// \param p - reference to Parent
    ExistingMethodsControl(Parent & p) : parent_(p), current_(0), id_of_selected_(db::id<PricingMethod>::notExist()) {}

    /// \return currently selected existing pricing method instance
    PricingMethod * getValue() const 
    {
        return current_;
    }

    /// \return db sheet ID of currently selected existing pricing method instance
    db::id<PricingMethod> getId() const
    {
        return id_of_selected_;
    }

protected:

    /// Fills the combo-box with labels of existing pricing method instances.
    void fillList()
    {
        // current mode
        lib::Mode m = parent_.getMode();

        // current model
        Model * aModel = parent_.get<Model>();
        // current option
        Option * anOption = parent_.get<Option>();
        // current family
        Family * aFamily = parent_.getFamily();
        // current pricing
        Pricing * aPricing = parent_.getPricing();

        // for each existing pricing method
        for (db::Storage<PricingMethod>::iterator i(db::theMethods()); i; ++i)
        {
            // get i-th pricing method entry
            db::Entry<PricingMethod> e = db::theMethods().at(*i);

            // .. and a pointer to the underlying premia structure
            PricingMethod * method = e.getData();

            // if a current pricing contains the method
            if (contains(aPricing, method))
            {
                // ... check it for compatibility with current option and model
                if (method->CheckOpt(anOption, aModel) == OK)
                {
                    // ... and if it is compatible add it to the list remembering its db sheet id
                    this->AddString(e.getLabel(), (*i).getValue());
                }
            }
        }
    }

    /// Selection change event handler.
    /// \param data - db sheet id of a currently selected existing pricing method instance
    void OnChanged(int data)
    {
        current_ = db::theMethods().at(id_of_selected_ = db::id<PricingMethod>(data)).getData();
    }

private:
    /// reference to parent control
    Parent        & parent_;
    /// current pricing method
    PricingMethod * current_;
    /// id of a current pricing method
    db::id<PricingMethod> id_of_selected_;
};
}}