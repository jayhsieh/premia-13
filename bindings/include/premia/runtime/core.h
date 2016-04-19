#ifndef _PREMIA_API_CORE_H_INCLUDED_
#define _PREMIA_API_CORE_H_INCLUDED_

/*!	\file	core.h
	\brief	Defines the elementary API to Premia

	\warning To be included into only one translation unit of a project since contains some initialization code
*/

#ifdef _MSC_VER
#include <direct.h>
#endif


#include <boost/scoped_ptr.hpp>
#include <premia/import.h>
#include <premia/runtime/VarIterator.h>
#include <premia/runtime/input_check.h>

/// Contains Premia-related C++ stuff
namespace premia	{

/// Contains Premia API
namespace api		{

	/// Implementation details
	namespace details	{

		/// Selected asset Id
		int g_asset_id = 0;

		/// Selected model
		::Model *	g_model = 0;
		/// Selected option family
		::Family*   g_family = 0;
		/// Selected option
		::Option*   g_option = 0;
		/// Selected pricing
		::Pricing*  g_pricing = 0;
		/// Selected pricing method
		::PricingMethod* g_method = 0;
		
		VAR*  g_params = 0;

		/// If in parameter writing mode, points to the current writer
		boost::scoped_ptr<VarIterator>		g_var_iterator;

		/// \brief Utility function to check that the API is in writing mode
		void check_iterating()
		{
			if (!g_var_iterator)
				throw exception("Unexpected write operation");
		}

		/// Initialization code to be called automatically before any work with Premia
		struct StartUp
		{
			StartUp()
			{
				InitVar();
			}
		} startup;

	}

	using namespace details;


	/// \brief Initializes Premia
	/// \param base_path path to Premia 'data' directory 
	void init_premia(const char *base_path)
	{
		chdir(base_path);
		strcpy(premia_data_dir, base_path);
		path_sep = "/";
	}

	/// \brief sets up current asset
	void setCurrentAsset(int asset_type)
	{
		g_asset_id = asset_type;
	}

	/// \brief sets up current model and switches API to parameter writing mode
	/// \param model_id model index
	void  setCurrentModel(int model_id)
	{
		if (g_var_iterator)
			throw exception("unexpected setCurrentModel");
		g_model = premia_assets[g_asset_id].models[model_id];
		(*g_model->Init)(g_model);
		g_params = reinterpret_cast<VAR*>(g_model->TypeModel);
		g_var_iterator.reset(new VarIterator(g_params, g_model->nvar));
	}

	void  readCurrentModel()
	{
		if (g_var_iterator)
			throw exception("unexpected readCurrentModel");

		(*g_model->Init)(g_model);
		g_var_iterator.reset(new VarIterator(reinterpret_cast<VAR*>(g_model->TypeModel), g_model->nvar));
	}

	/// \brief sets up current option and switches API to parameter writing mode
	/// \param family_id  option family index
	/// \param option_id  option index in the family
	/// \warning should be called after having set up a model
	void  setCurrentOption(int family_id, int option_id)
	{
		if (g_var_iterator)
			throw exception("unexpected setCurrentOption");
		g_family = premia_assets[g_asset_id].families[family_id];
		g_option = (*g_family)[option_id];

		(*g_option->Init)(g_option, g_model);
		g_params = reinterpret_cast<VAR*>(g_option->TypeOpt);
		g_var_iterator.reset(new VarIterator(g_params, g_option->nvar));
	}
	
	void readCurrentOption()
	{
		if (g_var_iterator)
			throw exception("unexpected setCurrentOption");

		(*g_option->Init)(g_option, g_model);
		g_var_iterator.reset(new VarIterator(reinterpret_cast<VAR*>(g_option->TypeOpt), g_option->nvar));
	}

	/// \brief sets up current method and switches API to parameter writing mode
	/// \param pricing_id  pricing index
	/// \param method_id  method index in the pricing
	/// \warning should be called after having set up an option
	void  setCurrentMethod(int pricing_id, int method_id)
	{
		if (g_var_iterator)
			throw exception("unexpected setCurrentMethod");
		g_pricing = premia_assets[g_asset_id].pricings[pricing_id];
		g_method = g_pricing->Methods[method_id];
		(*g_method->Init)(g_method, g_option);
		g_params = g_method->Par;
		g_var_iterator.reset(new VarIterator(g_method->Par, MAX_PAR));
	}
	
	void  readCurrentMethod()
	{
		if (g_var_iterator)
			throw exception("unexpected readCurrentMethod");

		(*g_method->Init)(g_method, g_option);
		g_var_iterator.reset(new VarIterator(g_method->Par, MAX_PAR));
	}
	
    template <class T>
    struct Assign {         
        Assign(T x) :x(x) {}
        void operator() (T *p, VAR * pv) 
        {
            *p = x; 
            if (pv->setter)
                (*pv->setter)(g_params);
        }
        T x;
    };

    template <class T> inline Assign<T> assign(T x) { return Assign<T>(x); }
    
    struct Nope 
    {
        template <class T> void operator () (T const &) {}        
        template <class T> void operator () (T const &, VAR * pv) {}        
        template <class T, class U> void operator () (T const &, U const &) {}        
    } nope;
	
    template <class T>
        struct Read {
            Read(T *value) : value(value) {}
            void operator() (T *p, VAR * pv) { *value = *p; }
        private:
            T   *value;
        };
        
    template <class T>
        Read<T> read(T *p) 
        {
            return Read<T>(p);
        }

	/// \brief writes a double value to the current parameter and goes to the next one
	void  write_double(double x)
	{
		check_iterating();
		g_var_iterator->process_double(assign(x));
	}
	
	void ignore_double(double x)
	{
		check_iterating();
		g_var_iterator->process_double(nope);
	}

	double read_double()
	{
	    double x;
		check_iterating();
		g_var_iterator->process_double(read(&x));
		return x;
	}

	/// \brief writes a long int value to the current parameter and goes to the next one
	void  write_long(long x)
	{
		check_iterating();
		g_var_iterator->process_long(assign(x));
	}
	
	void  ignore_long(long x)
	{
		check_iterating();
		g_var_iterator->process_long(nope);
	}

	long read_long()
	{
	    long x;
		check_iterating();
		g_var_iterator->process_long(read(&x));
		return x;
	}
	

	/// \brief writes an int value to the current parameter and goes to the next one
	void  write_int(int x)
	{
		check_iterating();
		g_var_iterator->process_int(assign(x));
	}
	
	void  ignore_int(int x)
	{
		check_iterating();
		g_var_iterator->process_int(nope);
	}

	int read_int()
	{
	    int x;
		check_iterating();
		g_var_iterator->process_int(read(&x));
		return x;
	}
	
	inline void assign_filename(char ** p, VAR * pv, const char *x)
	{
		if (*p) free(*p);
		*p = strdup(x); 
	}

	/// \brief writes a string/filename value to the current parameter and goes to the next one
	void  write_filename(const char * x)
	{
		check_iterating();
		g_var_iterator->process_filename(boost::bind(assign_filename, _1, _2, x));
	}
	
	void  ignore_filename(const char *)
	{
		check_iterating();
		g_var_iterator->process_filename(nope);
	}

	const char* read_filename()
	{
	    char* x;
		check_iterating();
		g_var_iterator->process_filename(read(&x));
		return x;
	}
	

	/// \brief writes a value of enumeration type to the current parameter and goes to the next one
	void  write_enum(int x)
	{
		check_iterating();
		g_var_iterator->process_enum(assign(x));
	}
	
	void  ignore_enum(int x)
	{
		check_iterating();
		g_var_iterator->process_enum(nope);
	}
	
	
	int read_enum()
	{
	    int x;
		check_iterating();
		g_var_iterator->process_int(read(&x));
		return x;
	}
	
    static int constant(int x, VAR *v) { return x; } 
    
    std::vector<double>     array_members_;
    
	//! \brief starts writing an array into the current VAR
	void write_array(int size)
	{
		check_iterating(); 
		g_var_iterator->process_array_begin(boost::bind(constant, size, _1));
	}
	
	void ignore_array(int size)
	{
		check_iterating(); 
		g_var_iterator->process_array_begin(boost::bind(constant, size, _1));
	}

	struct ReadArraySize {
	    ReadArraySize(int * x) : x(x) {}
	    int operator () (VAR *v) {
	        switch (v->Vtype) {
	            case PNLVECT: *x = v->Val.V_PNLVECT->size; break;
	            case PNLVECTCOMPACT: *x = v->Val.V_PNLVECTCOMPACT->size; break;
	            default: throw exception("VAR is expected to be either PNLVECT or PNLVECTCOMPACT");
	        }
	        return *x;
	    }
	    private:
	        int *x;
	};
	
	int read_array_size()
	{
		check_iterating();
		int x;
		g_var_iterator->process_array_begin(ReadArraySize(&x));
		return x;
	}

	static void assign_array_double(std::vector<double> *dst, double x, VAR*, int idx)
	{
	    dst->push_back(x);
	}
	
	void write_array_double(double x)
	{
		check_iterating();
	    g_var_iterator->process_array_double(boost::bind(assign_array_double, &array_members_, x, _1, _2));
	}
	
	void ignore_array_double(double x)
	{
		check_iterating();
	    g_var_iterator->process_array_double(nope);
	}

	struct ReadArrayDouble {
	    ReadArrayDouble(double *x) : value(x) {}	    
	    void operator () (VAR *v, int idx) {
	        switch (v->Vtype) {
	            case PNLVECT: *value = pnl_vect_get(v->Val.V_PNLVECT, idx); break;
	            case PNLVECTCOMPACT: *value = pnl_vect_compact_get(v->Val.V_PNLVECTCOMPACT, idx); break;
	            default: throw exception("VAR is expected to be either PNLVECT or PNLVECTCOMPACT");
	        }
	    }
	private:
	    double *value;
	};
	
	double read_array_double()
	{
	    check_iterating();
	    double x;
	    g_var_iterator->process_array_double(ReadArrayDouble(&x));
	    return x;
	}
	
	static void flush_array(std::vector<double> *array_members, VAR* v)
	{
		if (v->Vtype == PNLVECT)
		{
			pnl_vect_free(&v->Val.V_PNLVECT);
			v->Val.V_PNLVECT = pnl_vect_create_from_ptr(array_members->size(), &(*array_members)[0]);
		}

		if (v->Vtype == PNLVECTCOMPACT)
		{
			pnl_vect_compact_free(&v->Val.V_PNLVECTCOMPACT);
			v->Val.V_PNLVECTCOMPACT = pnl_vect_compact_create_from_ptr(array_members->size(), &(*array_members)[0]);
		}

		array_members->clear();
	}
	
	//! \brief all values for the array are accumulated so we can create pnl_vect and modify current VAR
	void write_array_end()
	{
		check_iterating();
	    g_var_iterator->process_array_end(boost::bind(flush_array, &array_members_, _1));
	}
	
	void ignore_array_end()
	{
		check_iterating();
	    g_var_iterator->process_array_end(nope);
	}

	void read_array_end()
	{
		check_iterating();
	    g_var_iterator->process_array_end(nope);
	}

	/// \brief utility function that interrupts parameter writing mode
	/// \warning to be used only in REPL mode
	void reset()
	{
		g_var_iterator.reset(0);
	}

	/// \brief checks that all parameters for the current entity have been written and switches from parameter writing mode
	void  stopWriteParameters()
	{
		if (!g_var_iterator)
			throw exception("unexpected stopWriteParameters");
		if (!g_var_iterator->finished())
			throw exception("not enough parameters for VAR");
		g_var_iterator.reset(0);
	}

	/// \brief Runs the current pricing method on the current model and the current option
	/// \warning should be called after having set up a model, an option and a method
	void  compute()
	{
		are_there_any_errors_in_params(g_model);
		are_there_any_errors_in_params(g_option);
		are_there_any_errors_in_params(g_method);
		checkMixing(g_option, g_model, g_pricing);

		if ((*g_method->CheckOpt)(g_option, g_model) != OK)
			throw exception("the method cannot work on an option of this type");

		(*g_method->Compute)(g_option->TypeOpt, g_model->TypeModel, g_method);
	}

	/// \brief retrieves a pricing method result of type double
	/// \param i result parameter index
	/// \warning should be called after compute
	double  get_result_double(int i)
	{
		return g_method->Res[i].Val.V_DOUBLE;
	}

	/// \brief retrieves a pricing method result of type bool
	/// \param i result parameter index
	/// \warning should be called after compute
	bool  get_result_bool(int i)
	{
		return g_method->Res[i].Val.V_BOOL != 0;
	}

	/// \brief retrieves a pricing method result of array of double type 
	/// \param i result parameter index
	/// \return the array's size
	/// \warning should be called after compute
	int  get_result_array_size(int i)
	{
		return g_method->Res[i].Val.V_PNLVECT->size;
	}

	/// \brief retrieves a pricing method result of array of double type 
	/// \param i result parameter index
	/// \param j index in the result array
	/// \return j-th item of the array
	/// \warning should be called after compute
	double  get_result_array_item(int i, int j)
	{
		return g_method->Res[i].Val.V_PNLVECT->array[j];
	}

}}

#endif
