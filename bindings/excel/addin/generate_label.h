#pragma once

#include "storage.h"

namespace prxl {
namespace db
{
    /// Searches for an entity of type T with given 'label'.
    /// \param label - a label of en entity to be found
    /// \returns ID of an entity if found and -1 otherwise
    template <class T>
        id<T> lookupLabel(_bstr_t label)
    {
        for (Storage<T>::iterator i(theStorage<T>()); i; ++i)
        {
            if (getLabel(theStorage<T>().at(*i)) == label)
                return *i;
        }

        return id<T>::notExist();
    }


    /// Generates a unique label for an entity instance.
    /// \param prefix - an entity type string, e.g. BlackScholes1dim
    /// \returns a label for the entity
    template <class T>
    _bstr_t generateNewLabel(_bstr_t prefix)
    {
        for (int i = 0;; ++i)
        {
            _bstr_t label = prefix + ":" + toString(i);

            if (lookupLabel<T>(label) == id<T>::notExist())
                return label;
        }

        return "...strange...";
    }
}}
