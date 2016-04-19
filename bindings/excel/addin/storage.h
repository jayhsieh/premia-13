#pragma once

#include "entry.h"

namespace prxl {
namespace db
{
    /// Class representing a db sheet.
    /// Contains an array of Entry<T> representing a single db row
    template <class T>
    struct Storage
    {
        /// A constructor reading the storage content from the db sheet
        Storage() 
        {
            Excel::RangePtr O = db::getRow(role_of<T>(), 0);

            // number of rows
            int N = O->Value2;

            // for each row 
            for (int i = 0; i != N; ++i)
            {
                // read it and push its representation into the array
                storage_.push_back(Entry<T>(O->GetOffset(i,0)));        
            }
        }

        struct iterator
        {
            iterator(Storage<T> & s) : idx_(0), stopper_(s.size().getValue()) {}

            id<T> operator * () const { return id<T>(idx_); }

            operator bool () const { return idx_ != stopper_; }

            iterator & operator ++ () { ++idx_; return *this; }

        private:
            int         idx_;
            const int   stopper_;
        };

        /// adds an entry both to the array and to the db sheet
        /// \param e - an entry to be added
        id<T> push(Entry<T> const & e)
        {
            Excel::RangePtr O = getRow(role_of<T>(), 0);
            // number of rows in db sheet
            int N = O->Value2;
            // add the entry to the array
            storage_.push_back(e);
            // serialize it into the db sheet
            storage_.back().serialize(O->GetOffset(N,0));
            // increment db sheet rows count
            O->Value2 = N + 1;
            return id<T>(N);
        }

        Excel::RangePtr locate(id<T> idx) const
        {
            return getRow(role_of<T>(), idx.getValue());            
        }

        /// clears the array (but not db sheet).
        /// NB! to be used only at the add-in shutdown to release resources acquired
        void clear()
        {
            storage_.clear();
        }

        typedef std::vector<Entry<T> >  storage_t;

        /// looks up an entry with label 'name'
        /// \param name - a name of entry being looked for
        /// \return a pointer to it if found, otherwise 0
        Entry<T> const * lookup(_bstr_t const & name) const
        {
            for (storage_t::const_iterator it = storage_.begin(); it != storage_.end(); ++it)
            {
                if (name == it->getLabel())
                    return &*it;
            }

            return 0;
        }

        /// returns an entry at idx-th position in the array
        /// \param idx - entry index
        /// \return the entry
        Entry<T> const & at(id<T> idx) const
        {
            return storage_.at(idx.getValue());
        }

        /// returns an entry at idx-th position in the array
        /// \param idx - entry index
        /// \return the entry
        Entry<T>       & at(id<T> idx)      
        {
            return storage_.at(idx.getValue());
        }

        /// \return number of rows in the db sheet
        id<T> size() const
        {
            return id<T>((int)storage_.size());
        }

    private:
        /// an array of entries
        storage_t   storage_;
    };

    template <class T>
        std::map<IDispatch*, Storage<T> > & storageMap()
        {
            static std::map<IDispatch*, Storage<T> >  storage;
            return storage;
        }

    /// Meyer's singleton returning the instance of the storage for entities of type T
    template <class T>
        Storage<T> & theStorage()
    {
        return storageMap<T>()[g_app->ActiveWorkbook];
    }

    template <class T> 
        void cleanUpStorage()
        {
            std::map<IDispatch*, Storage<T> >   & m = storageMap<T>();
            std::map<IDispatch*, Storage<T> >::iterator it = m.find(g_app->ActiveWorkbook);
            if (it != m.end())
                m.erase(it);
        }

    /// Meyer's singleton returning the instance of the storage for created model instances
    inline Storage<Model>&  theModels()
    {
        return theStorage<Model>();
    }

    /// Meyer's singleton returning the instance of the storage for created option instances
    inline Storage<Option>&  theOptions()
    {
        return theStorage<Option>();
    }

    struct Result;
    Storage<Result> & theResults();

    /// Meyer's singleton returning the instance of the storage for created pricing method instances
    inline Storage<PricingMethod>&  theMethods()
    {
		//  theResults();   // to make sure that all methods data is initialized (????)
        return theStorage<PricingMethod>();
    }
}}