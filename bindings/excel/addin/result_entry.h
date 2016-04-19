#pragma once

#include "premia.h"
#include "storage.h"

namespace prxl {
namespace db
{
    /// Tag class for result regions.
    struct Result {};

    /// An entry of the db containing problem results.
    /// Besides fields inherited from EntryBase 
    /// it contains IDs of the problem model, option and method
    template <> 
    struct Entry<Result> : EntryBase
    {
        /// A constructor for just created problems.
        /// \param r - top-left corner of the result region
        /// \param h - number of rows in the result region
        /// \param m - the problem mode
        /// \param model_id - the problem model ID
        /// \param option_id - the problem option ID
        /// \param method_id - the problem method ID
        Entry(Excel::RangePtr r, int h, lib::Mode m, 
            id<Model> model_id, id<Option> option_id, id<PricingMethod> method_id)
            :   EntryBase("",r,h,m)
            ,   model_id_(model_id)
            ,   option_id_(option_id)
            ,   method_id_(method_id)
        {}

        /// A constructor restoring a result entry from a db sheet row.
        /// \param O - a row
        Entry(Excel::RangePtr O)
            :   EntryBase(O)
            ,   model_id_(At(O, fld::ModelId)->Value2)
            ,   option_id_(At(O, fld::OptionId)->Value2)
            ,   method_id_(At(O, fld::MethodId)->Value2)
        {
            // we have read all parameters
            // we know which model and option corresponds to the given method
            // and now we may initialize data_ field in its entry

            // the method entry
            Entry<PricingMethod>  &method_entry = theMethods().at(method_id_); 

            // corresponding model
            Option * opt = theOptions().at(option_id_).getData();
            // and option
            Model * mod = theModels().at(model_id_).getData();

            // get pricing for the model and the option
            Pricing * pr = lib::getPricing(getMode(), mod, opt);

            for (PricingMethod ** pm_it = pr->Methods; *pm_it; ++pm_it)
            {
                // if the pricing contains a method with same name
                if (_bstr_t((*pm_it)->Name) == method_entry.getName())
                {
                    // check whether it is compatible with the model and the option
                    if ((*pm_it)->CheckOpt(opt, mod) == OK)
                    {
                        method_entry.setData(*pm_it);
                        break;
                    }
                }
            }
        }

        /// serializes the entry into the db sheet row.
        /// \param row - a row where the problem description should be serialized
        void serialize(Excel::RangePtr row) const
        {
            EntryBase::serialize(row);

            At(row, fld::ModelId)->Value2 = model_id_.getValue();
            At(row, fld::OptionId)->Value2 = option_id_.getValue();
            At(row, fld::MethodId)->Value2 = method_id_.getValue();
        }

        /// returns problem model ID
        id<Model> getModelId() const { return model_id_; }

        /// return problem option ID
        id<Option> getOptionId() const { return option_id_; }

        /// returns problem method ID
        id<PricingMethod> getMethodId() const { return method_id_; }

    private:
        /// the problem model ID
        id<Model>     model_id_;
        /// the problem option ID
        id<Option>     option_id_;
        /// the problem method ID
        id<PricingMethod>     method_id_;
    };

    /// the storage of problem descriptions
    inline Storage<Result>&  theResults()
    {
        return theStorage<Result>();
    }

    /// returns Compute button position for idx-th problem
    /// \param idx - a problem index
    /// \return Compute button position
    inline Excel::RangePtr getButtonPosition(id<Result> idx)
    {
        Entry<Result> const & e = theResults().at(idx);
        return e.getTopLeft()->GetOffset(e.getHeight(), 1);
    }

}}