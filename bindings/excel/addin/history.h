#pragma once

namespace prxl {
namespace wzd {

template <class T> struct HistoryChooser;

/// A selection history for entities of type T
/// \param T - entity type
template <class T>
struct History
{
    typedef std::list<T>  history_t;

    friend struct HistoryChooser<T>;

    /// adds an element to the history
    /// \param x - an element to be added
    void push(T const & x) { history_.push_front(x); }

    /// \return most recently added element
    T    top() const { return history_.empty() ? 0 : history_.front(); }

private:
    
    history_t   history_;
};

/// A class for choosing the most recent element from a selection history
/// \param T -- element type
template <class T>
struct HistoryChooser
{
    typedef typename History<T>::history_t  history_t;

    /// A constructor.
    /// \param h - selection history
    HistoryChooser(History<T> const & h) 
        : history_(h.history_), stopper_(history_.end())
    {}

    /// Looks in the history for an element
    /// and if finds one better then the best updates the best one and returns true
    /// \param x - an element to be tested
    /// \return true iff x is better than the best
    bool contains(T const & x)
    {
        for (typename history_t::const_iterator it = history_.begin(); it != stopper_; ++it)
        {
            if (*it == x)
            {
                stopper_ = it;
                return true;
            }
        }

        return false;
    }

private:
    history_t const &  history_;
    /// the best element
    typename history_t::const_iterator stopper_;
};
}}