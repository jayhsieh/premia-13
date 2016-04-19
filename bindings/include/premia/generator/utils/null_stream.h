#ifndef _premia_generator_utils_null_stream_h_included_
#define _premia_generator_utils_null_stream_h_included_

#include <iostream>

namespace premia {
namespace pygen  {

    struct null_stream : std::basic_ostream<char>
    {
        struct null_buf : std::basic_streambuf<char>
        {
            std::streamsize xsputn(const char *, std::streamsize c){
                return c;
            } 
        } buf;
        
#ifdef _MSC_VER
        null_stream() : std::basic_ostream<char>(&buf) {}
#else
        null_stream()
        {
            rdbuf(&buf);
        }
#endif
    };    
    
}}

#endif
