#pragma once

#include "premia.h"
#include "entry_base.h"

namespace prxl {
namespace db
{
    /// Searches for a premia entity instance by its name.
    /// Specialized for concrete object types
    template <class T> T * createEntity(lib::Mode m, _bstr_t const & name);

    /// Searches for a premia model instance by its name
    /// \param m - mode
    /// \param name - name of the model to be returned
    /// \returns a pointer to Premia model instance
    template <> inline Model * createEntity<Model>(lib::Mode m, _bstr_t const & name)
    {
        return lib::getModel(m, name);
    }

    /// Searches for a premia option instance by its name
    /// \param m - mode
    /// \param name - name of the option to be returned
    /// \returns a pointer to Premia option instance
    template <> inline Option * createEntity<Option>(lib::Mode m, _bstr_t const & name)
    {
        return lib::getOption_2(m, name);
    }

    /// searches a premia pricing method instance by its name.
    /// \param m - mode
    /// \param name - name of a pricing method to be returned
    /// \returns null pointer since (m,name) is not sufficient for identifying a pricing method
    template <> inline PricingMethod * createEntity<PricingMethod>(lib::Mode m, _bstr_t const & name)
    {
        return 0;
    }

    /// A wrapper over premia objects (model, option, method).
    /// Stores: 
    /// - Pointer to the underlying premia object
    /// - fields inherited from EntryBase
    template <class T>
    struct Entry : EntryBase
    {
        /// A constructor for just created entities.
        /// \param n - name of the entity
        /// \param r - pointer to the top-left corner of the entity region
        /// \param h - number of rows in the entity region
        /// \param data - a pointer to the underlying premia object 
        Entry(_bstr_t const & n, Excel::RangePtr r, int h, lib::Mode m, T * data)
            :   EntryBase(n,r,h,m), data_(data)
        {}

        /// A constructor restoring entities from a db sheet.
        /// Fills its fields reading a db sheet row
        /// \param O - a db sheet row to be read
        Entry(Excel::RangePtr O)
            :   EntryBase(O)
        {
            data_ = createEntity<T>(getMode(), getName());
        }

        /// Stores the pointer to the premia object.
        /// This function is used only for pricing methods 
        /// \param x - a pointer to the underlying premia object
        void setData(T * x)
        {
            data_ = x;
        }

        /// Returns the pointer to the underlying premia object.
        /// This function is used only for pricing methods 
        /// \returns the pointer to the underlying premia object.
        T * getData() const
        {
            return data_;
        }


    private:
        /// Pointer to the underlying premia object
        T            *  data_;
    };
}}