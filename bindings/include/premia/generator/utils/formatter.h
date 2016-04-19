#pragma once

#include <map>
#include <string>
#include <list>
#include <istream>
#include <ostream>
#include <exception>

#include <boost/format.hpp>
#include <boost/array.hpp>
#include <boost/foreach.hpp>
#include <boost/filesystem.hpp>
#include <boost/filesystem/fstream.hpp>

namespace premia {
namespace pygen {
namespace formatter_dsl {

	namespace fs = boost::filesystem;

	/// \brief pretty printer for a source code to a file
	/// Tracks current indentation and substitutes variables
	struct Formatter 
	{
		/// mapping from variable names to their values
		typedef std::map<std::string, std::string> Variables;

		/// constructs a formatter living on a given file
		/// \param p path to the file
		/// \param indent_size size for indentation (default value: 3)
		explicit Formatter(fs::path const & p, int indent_size = 3)
			: indent_size_(indent_size), indent_(0), buf_(), out_(p)
		{
		    if (!out_)
		        throw std::logic_error("Unable to open " + p.string() + " for writing");
		
			for (int i = 0; i < 200; ++i) buf_[i] = ' ';
			setIndent(0);
			push_scope();
		}

		/// introduces a new variables scope 
		/// (variables added afterwards will be discarded by calling to pop_scope)
		void push_scope()
		{
			variables_.push_front(Variables());
		}

		/// discards the current variables scope
		void pop_scope()
		{
			variables_.pop_front();
		}
		
		/// increments indentation level
		void incindent() 
		{
			setIndent(indent_ + indent_size_);
		}

		/// decrements indentation level
		void decindent() 
		{
			setIndent(indent_ - indent_size_);
		}

		/// prints a single line given by const char *
		friend Formatter& operator << (Formatter & out, const char * str)
		{
			return out.print(str);
		}

		/// prints a single line given by std::string
		friend Formatter& operator << (Formatter & out, std::string const & str)
		{
			return out.print(str.c_str());
		}

		/// prints a single line given by boost::format
		template<class Ch, class Tr, class Alloc>  
		friend Formatter& operator << (Formatter & out, const boost::basic_format<Ch, Tr, Alloc>& f) 
		{
			return out.print(f.str().c_str());
		}

		/// adds a variable to current scope
		/// \param T a type convertible to std::string by boost::lexical_cast<std::string>
		/// \param name name of the variable, it is possible to refer on it in a printable string as %name%
		/// \param value a value string representation of which will be bound to name
		template <class T>
			Formatter& operator () (const char * name, T const & value)
			{
				variables_.front()[name] = make_substitutions(str(value));
				return *this;
			}
			
		Formatter& process_file(boost::filesystem::path const & p)
		{
		    if (!exists(p))
		         throw std::logic_error(p.string() + " doesn't exist");
		         
		    boost::filesystem::ifstream input(p);
		    
		    return *this << input;
		}

		/// loads strings from std::istream and formats them
		friend Formatter& operator << (Formatter & out, std::istream & input)
		{
			std::string buf;
			while (getline(input, buf))
			{
				out << buf;
			}
			return out;
		}

		void putch(char ch)
		{
			out_ << &buf_[0] << ch << "\n";
		}

		/// looks up a variable
		std::string const & lookupVar(std::string const & var) const
		{
			BOOST_FOREACH(Variables const & scope, variables_)
			{
				Variables::const_iterator it = scope.find(var);

				if (it != scope.end())
					return it->second;		
			}

			throw std::logic_error(std::string("unknown variable '") + var  + "'");
		}

		void changeVar(std::string const & var, std::string const & newVal) 
		{
			BOOST_FOREACH(Variables & scope, variables_)
			{
				Variables::iterator it = scope.find(var);

				if (it != scope.end())
			    {
			        it->second = newVal;
			        return;		
			    }
			}

			throw std::logic_error(std::string("unknown variable '") + var  + "'");
		}

	private:

		/// makes variable substitutions in a string
		/// \param s string to transform
		/// \return transformed string
		std::string make_substitutions(std::string const & s) const
		{
			std::stringstream ss;

			print_replaced(ss, s.c_str());

			return ss.str();
		}

		/// makes variable substitutions in a string and output it to a std::ostream
		void print_replaced(std::ostream & out, const char * str) const
		{
			while (*str)
			{
				const char * p = str;

				for (; *p && *p != '%'; ++p);

				out.write(str, p - str);

				if (*p == '%')
				{
					const char *q = p + 1;
					
					if (*q == '%')
					    out << "%";
					else
					{
					    for (; *q && *q != '%'; ++q);

					    if (*q == 0)
						    throw std::logic_error(std::string("unmatched % in string '") + str + "'");

					    std::string  v(p + 1, q);

					    out << lookupVar(v);
                    }
					str = q + 1;
				}
				else 
					str = p;
			}

		}

		/// prints a single line
		Formatter& print(const char * str)
		{
			// makes an indentation
			out_ << &buf_[0];

			// prints the string replacing variables in it
			print_replaced(out_, str);
			
			out_ << "\n";
			return *this;
		}

		/// changes indentation level
		void setIndent(int newIndent)
		{
			if (newIndent > 199) throw std::logic_error("too many indents...");
			if (newIndent   < 0) throw std::logic_error("indent level is below 0...");
			buf_[indent_] = ' ';
			buf_[indent_ = newIndent] = 0;
		}


	private:
		int const				indent_size_;	//!< indentation size
		int						indent_;		//!< current indentation level
		boost::array<char, 200>	buf_;			//!< holds a string of spaces representing current indentation

		fs::ofstream			out_;			//!< file stream for output
		std::list<Variables>	variables_;		//!< stack of variable scopes. the global scope is at the end of the list
	};
}

using formatter_dsl::Formatter;

}}
