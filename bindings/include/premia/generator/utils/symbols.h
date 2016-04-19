#pragma once

#include <set>
#include <map>
#include <string>

#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/assign/list_of.hpp>


namespace premia {
namespace pygen {

	/// Contains utilities for assuring identifier correctness according to target language rules
	namespace symbol_utils
	{
		namespace ba = boost::algorithm;

		/// assures that the first character of an identifier would be valid
		/// \param ch symbol to correct
		/// \return corrected symbol
		inline char correctFirstSymbol(char ch)
		{
			return ba::is_alpha()(ch) ? ch : '_';
		}

		/// assures that a non-first character of an identifier would be valid
		/// \param ch symbol to correct
		/// \return corrected symbol
		inline char correctSymbol(char ch)
		{
			return ba::is_alnum()(ch) ? ch : '_';
		}

		/// makes an identifier valid
		/// \param src an original identifier
		/// \return valid identifier
		inline std::string replaceNonAlnumCharacters(const char * src)
		{
			if (*src == 0)
				return "_";

            char last;
			std::string dst;
			dst.push_back(last = correctFirstSymbol(*src++));
			
			for (; *src; ++src)
			{
			    char ch = correctSymbol(*src);
				if (last == '_' && ch == '_')
				    ; // skip it
				else
    				dst.push_back(ch);
    			last = ch;
			}

			return dst;
		}
		/// converts a value of type T to std::string using internally std::stringstream
		template <class T>
		inline std::string str(T const & x)
		{
			return boost::lexical_cast<std::string>(x);
		}

		/// maps used id to its repetition number
		typedef std::map<std::string, int> UsedIds;

		using boost::assign::list_of;

		/// list of forbidden identifiers. to be extended and parameterized by target language
		const std::set<std::string> python_keywords = list_of("lambda");

		// TODO: parametrize by naming conventions for the target language

		/// makes an identifier valid
		/// \param src original identifier
		/// \param used_names set of already used names
		/// \return corrected identifier
		inline std::string correctIdentifierName(std::string const &src, UsedIds & used_names)
		{
			// making the identifier valid according to the target language rules
			std::string alnum = symbol_utils::replaceNonAlnumCharacters(src.c_str());

			// checking for collision with reserved names for the target language
			if (python_keywords.find(alnum) != python_keywords.end())
				used_names[alnum] = 0;

			// checking whether the name was already used
			UsedIds::iterator it = used_names.find(alnum);

			// if used
			if (it != used_names.end())
			{
				std::string s;

				do {
					// generate new numerical suffix
					s = alnum + "_" + str(++it->second);
					// and check uniqueness of generated name
				} while (used_names.find(s) != used_names.end());

				alnum = s;
			}
			else
				// remembering that the name is used
				used_names[alnum] = 0;

			return alnum;
		}
	}

	using symbol_utils::str;
	
	namespace api 
	{
		using symbol_utils::UsedIds;
		using symbol_utils::correctIdentifierName;	
	}
}}
