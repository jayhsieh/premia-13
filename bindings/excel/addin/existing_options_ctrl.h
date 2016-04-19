#pragma  once

namespace prxl {
namespace wzd {
/// A wrapper over a list-box adding a logic for selecting an existing option instance.
/// \param Parent - a class that can be queried for a current model
template <class Parent>
struct ExistingOptionsControl 
    :   ComboboxOfExisting<IDC_OPTIONS_EXISTING>
{
    /// A constructor.
    /// \param p - reference to Parent
    ExistingOptionsControl(Parent & p) 
        : parent_(p), current_(0), family_(0), pricing_(0), id_of_selected_(db::id<Option>::notExist()) {}


    /// \return currently selected existing option instance
    Option * getValue() const 
    {
        return current_;
    }

    /// \return a family of currently selected existing option instance
    Family * getFamily() const
    {
        return family_;
    }

    /// \return a pricing for currently selected existing model option and a model type
    Pricing * getPricing() const 
    {
        return pricing_;
    }

    /// \return db sheet ID of currently selected existing option instance
    db::id<Option> getId() const
    {
        return id_of_selected_;
    }

protected:

    /// Fills the combo-box with labels of existing pricing method instances.
    void fillList()
    {
        // current model
        Model * m = parent_.get<Model>();
        // current mode
        lib::Mode mode = parent_.getMode();

        // for each existing option instance
        for (db::Storage<Option>::iterator i(db::theOptions()); i; ++i)
        {
            // get its entry
            db::Entry<Option> const & e = db::theOptions().at(*i);

            // if it is in proper mode
            if (e.getMode() == mode)
            {
                Option * opt = e.getData();
                Family * f = lib::getFamily(mode, opt);

                // ... check whether it is compatible with current model
                if (Pricing * p = lib::getPricing(mode, m, (*f)[0]))
                {
                    opt->Init(opt, m);

                    // ... and if so then check are there any pricing methods 
                    for (PricingMethod ** pm_it = p->Methods; *pm_it; ++pm_it)
                    {
                        // ... compatible with the option and a current model
                        if ((*pm_it)->CheckOpt(opt, m) == OK)
                        {
                            // if so then add the option to the list
                            this->AddString(e.getLabel(), (*i).getValue());

                            break;
                        }
                    }
                }
            }
        }
    }

    /// Selection change event handler.
    /// \param data - db sheet id of a currently selected existing option instance
    void OnChanged(int data)
    {
        current_ = db::theOptions().at(id_of_selected_ = db::id<Option>(data)).getData();
        family_ = lib::getFamily(parent_.getMode(), current_);
        pricing_ = lib::getPricing(parent_.getMode(), parent_.get<Model>(), current_);
    }


private:
    /// reference to parent control
    Parent        & parent_;
    /// current option
    Option        * current_;
    /// family of current option
    Family        * family_;
    /// pricing for current option and model
    Pricing       * pricing_;
    /// db sheet id of current option
    db::id<Option> id_of_selected_;
};

}}