#ifndef _PREMIA_API_INPUT_CHECK_H_INCLUDED_
#define _PREMIA_API_INPUT_CHECK_H_INCLUDED_

/*! \file input_check.h
	\brief	Utilities to check input data
*/
#include <stdio.h>
#include <fstream>
#include <premia/import.h>
#include <premia/exception.h>

namespace premia {
namespace api {
namespace details {

	/// \brief checks that the option and the model are initialized mutually correctly for methods of the given pricing
	inline bool checkMixing(Option * option, Model * model, Pricing * pricing)
	{
		// if there are any errors when CheckMixing is called, they will be reported to 'out_stream' 
		g_dup_printf = true;
		FILE* h = fopen("checking_premia_input", "w");
                out_stream = h;

		int res = (*pricing->CheckMixing)(option, model);

		fclose(out_stream);
		g_dup_printf = false;

		if (res != OK)
		{
			std::ifstream infile("checking_premia_input");

			std::string s;
			char buf[256];

			// collecting lines from 'out_stream' file in order to provide readable error message
			while (infile)
			{
				infile.getline(buf, 256);
				s += buf;
				s += "\n";
			}

			throw premia::api::exception("CheckMixed failed: " + s);
		}

		return true;
	}

	/// \brief Launches premia's mechanism for checking an entity's parameters
	/// \param Target - an entity type
	/// \param target - an entity 
	/// \return true iff there are no errors
	template <class Target>
	bool are_there_any_errors_in_params(Target * target)
	{
       	
		// if there any errors when Check method is called, they will reported to 'g_dup_file' 
		g_dup_printf = true;
		FILE *h = fopen("checking_premia_input", "w");
                g_dup_file = h;

		Planning plan;
		plan.Action = 'p';
		plan.VarNumber = 0;

		int result = target->Check(TOFILE, &plan, target);

		fclose(g_dup_file);
		g_dup_printf = false;

		if (result != OK)
		{
			std::ifstream infile("checking_premia_input");
			char buf[256];

			std::string err_msg;

			while (infile)
			{
				infile.getline(buf, 256);  
				std::string s = buf;

				if (s == "" || s == "All Right (ok: Return, no: n, h for Help) ? 	");
				else
					if (s == "Bad value:")
					{
						infile.getline(buf, 256);
						err_msg += buf;

						infile.getline(buf, 256);
						err_msg += ". ";
						err_msg += buf;
						err_msg += "\n";
					}
					else
					{
						err_msg += buf;
						err_msg += "\n";
					}
			}


			if (err_msg.size() > 0)
			{
				throw premia::api::exception("Some params have errors: " + err_msg);
			}
		}

		return  true;
	}

}}}


#endif
