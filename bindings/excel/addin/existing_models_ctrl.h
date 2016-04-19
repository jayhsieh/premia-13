#pragma once

#include "existing_combo.h"
#include "storage.h"

namespace prxl {
namespace wzd {
/// A wrapper over a list-box adding a logic for selecting an existing model instance.
/// \param Parent - a class that can be queried for a current mode
template <class Parent>
struct ExistingModelsControl : ComboboxOfExisting<IDC_MODELS_CREATED>
{
    /// A constructor.
    /// \param p - reference to Parent
    ExistingModelsControl(Parent & p) : parent_(p), current_(0), id_of_selected_(db::id<Model>::notExist()) {}

    /// \return db sheet ID of currently selected existing model instance
    db::id<Model> getId() const
    {
        return id_of_selected_;
    }

    /// \return currently selected existing model instance
    Model * getValue() const 
    { 
        return current_;
    }   

protected:

    /// Fills the combo-box with labels of existing pricing method instances.
    void fillList()
    {
        // current mode
        lib::Mode m = parent_.getMode();

        // for each existing model
        for (db::Storage<Model>::iterator i(db::theModels()); i; ++i)
        {
            // get its entry
            db::Entry<Model> const & model = db::theModels().at(*i);

            // and if it of suitable mode
            if (model.getMode() == m)
            {
                // ... add it the list
                this->AddString(model.getLabel(), (*i).getValue());
            }
        }
    }

    /// Selection change event handler.
    /// \param data - db sheet id of a currently selected existing model instance
    void OnChanged(int data)
    {
        current_ =  db::theModels().at(id_of_selected_ = db::id<Model>(data)).getData();
    }

private:
    /// reference to parent control
    Parent    & parent_;
    /// id of a current model
    db::id<Model> id_of_selected_;   
    /// current model
    Model   *   current_;
};
}}