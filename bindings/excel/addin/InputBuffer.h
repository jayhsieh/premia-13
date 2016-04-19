#pragma once

#include <string>
#include <map>
#include <list>

namespace prxl {
/// Makes from a key and a value a parameter string in Premia parameter format
/// \param key - a string containing a parameter key
/// \param value - a string containing a parameter value
/// \return a parameter string in Premia parameter format
inline std::string concatMKeyValue(std::string const & key, std::string const & value)
{
    std::string res = "m ";
    res += key;
    res += " ";
    res += value;

    return res;
}

/// Accumulates parameters in a string array for feeding them to FGet method
struct InputBuffer
{
#ifndef MAX_LINE
    static const int MAX_LINE = 200;
#endif

    /// A constructor
    InputBuffer()
    {
        for (int i = 0; i != MAX_LINE; ++i)
            buf_[i] = "";
    }

    /// Adds a parameter to the buffer.
    /// \param key - a parameter key
    /// \param value - a parameter value
    void PutLine(std::string const & key, std::string value)
    {
        // convert doubles from Excel format to Premia format
        std::replace(value.begin(), value.end(), ',', '.');
        data_[key] = value;
    }


    /// \return a string array with accumulated parameters in Premia format
    char ** getBuffer()
    {
        generated_.clear();

        int line_no = 0;

        for (std::map<std::string, std::string>::const_iterator it = data_.begin(); it != data_.end(); ++it)
        {
            generated_.push_back(concatMKeyValue(it->first, it->second));
            buf_[line_no++] = const_cast<char*>(generated_.back().c_str());
        }

        return buf_;
    }

private:
    char *buf_[MAX_LINE];
    std::map<std::string, std::string> data_;
    std::list<std::string>  generated_;
};
}