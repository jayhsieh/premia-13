#include <premia/runtime/core.h>
#include <boost/bind.hpp>

std::string g_error_msg;

template <class F>
	int wrap(F f)
	{
		try {
			f();
			return 0;
		} catch (std::exception & ex)
		{
			g_error_msg = ex.what();
			return -1;
		} catch (...) {
			g_error_msg = "Unknown error in PremiaDLL";
			return -1;
		}
	}

extern "C"
{
	PREMIADLL_API const char *  get_error_msg()
	{
		return g_error_msg.c_str();
	}

	int PREMIADLL_API  init_premia(const char * path)
	{
		return wrap(boost::bind(premia::api::init_premia, path));
	}

	int PREMIADLL_API setCurrentAsset(int asset_type)
	{
		return wrap(boost::bind(premia::api::setCurrentAsset, asset_type));
	}

	int PREMIADLL_API reset()
	{
		return wrap(premia::api::reset);
	}

	int PREMIADLL_API stopWriteParameters()
	{
		return wrap(premia::api::stopWriteParameters);
	}

	int PREMIADLL_API setCurrentModel(int model_id)
	{
		return wrap(boost::bind(premia::api::setCurrentModel, model_id));
	}

	int PREMIADLL_API setCurrentOption(int family_id, int option_id)
	{
		return wrap(boost::bind(premia::api::setCurrentOption, family_id, option_id));
	}

	int PREMIADLL_API setCurrentMethod(int pricing_id, int method_id)
	{
		return wrap(boost::bind(premia::api::setCurrentMethod, pricing_id, method_id));
	}

	int PREMIADLL_API write_double(double x)
	{
		return wrap(boost::bind(premia::api::write_double, x));
	}

	int PREMIADLL_API write_long(long x)
	{
		return wrap(boost::bind(premia::api::write_long, x));
	}

	int PREMIADLL_API write_int(int x)
	{
		return wrap(boost::bind(premia::api::write_int, x));
	}

	int PREMIADLL_API write_filename(const char * x)
	{
		return wrap(boost::bind(premia::api::write_filename, x));
	}

	int PREMIADLL_API write_enum(int x)
	{
		return wrap(boost::bind(premia::api::write_enum, x));
	}

	int PREMIADLL_API write_array(int sz)
	{
		return wrap(boost::bind(premia::api::write_array, sz));
	}

	// TODO: range check
	double  PREMIADLL_API get_result_double(int i)
	{
		return premia::api::get_result_double(i);
	}

	bool  PREMIADLL_API get_result_bool(int i)
	{
		return premia::api::get_result_bool(i);
	}

	int  PREMIADLL_API get_result_array_size(int i)
	{
		return premia::api::get_result_array_size(i);
	}

	double  PREMIADLL_API get_result_array_item(int i, int j)
	{
		return premia::api::get_result_array_item(i,j);
	}

	int PREMIADLL_API compute_3()
	{
		return wrap(premia::api::compute);
	}

};

