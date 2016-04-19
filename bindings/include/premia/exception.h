#ifndef _PREMIA_API_EXCEPTION_H_INCLUDED_
#define _PREMIA_API_EXCEPTION_H_INCLUDED_

/*! \file exception.h
	\brief Defines premia::api::exception
*/

#include <exception>
#include <string>

namespace premia {
namespace api	{

	//! \brief Base class for any Premia API exceptions
	//  TODO: use Boost.Exception
	struct exception : std::exception 
	{
		//! \brief constructs an exception given error message
		exception(std::string const & s) : s(s) {}

		//! \brief return a string representation of the exception
		const char * what() const throw () 
		{
			return s.c_str();
		}
                ~exception() throw () {}

	private:
		std::string const s;
	};

}}

#endif
