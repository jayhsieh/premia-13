#ifndef _fspath_included
#define _fspath_included

#include <istream>
#include <ostream>
#include <boost/filesystem.hpp>
namespace fs = boost::filesystem;

struct FsPath 
{
    FsPath(fs::path const & s = fs::path())
        : p(s)
    {}

    FsPath(const char * s)
        : p(s)
    {}

    friend std::istream & operator >> (std::istream &in, FsPath &x)
    {
        std::string ss;
        std::getline(in, ss);
        x.p = fs::path(ss);
        return in;
    }

    friend std::ostream& operator << (std::ostream &out, FsPath const& x)
    {
        std::string s = x.p.string();
        out << s;
        return out;
    }

    fs::path const & operator *() const 
    {
        return p;
    };

    fs::path const * operator ->() const 
    {
        return &p;
    };

private:
    fs::path p;
};

#endif