#ifndef _PREMIA_API_VARWRITER_H_INCLUDED_
#define _PREMIA_API_VARWRITER_H_INCLUDED_

/*!	\file VarWriter.h
	\brief Defines an output iterator into an array of VARs
*/

#include <stack>
#include <vector>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/scoped_ptr.hpp>
#include <premia/exception.h>

namespace premia	{
namespace api		{
namespace details	{

//! \brief Writer into an array of doubles
struct ArrayIterator 
{
	template <class Func>
	    ArrayIterator(VAR *v, Func f)
	    :   stopper_(f(v))
	    ,   idx_(0)
	    ,   v_(v)
	{
		if (stopper_ < 1) throw exception("Array size should be positive");
	}

	template <class Func>
	    void process_double(Func f)
	{
	    assert(idx_ < stopper_);
		f(v_, idx_);
		++idx_;
	}

	//! \brief checks if all data have been written into the array
	bool finished() { return stopper_ >= idx_; }

private:
	/// \brief how many elements are to be read
	int const stopper_;
	int       idx_;
	VAR * const v_;
};

//! \brief An output iterator into an array of VARs
struct VarIterator 
{
	/// Invariant:
	///    either in finished state or ready to write data
	///    finished = layers.empty()
	///    after every increment we validate the iterator

	//! \brief constructs a writer
	//! \param vars a pointer to the first VAR
	//! \param stopper maximal number of variables to be written
	VarIterator(VAR * vars, int stopper)
	{
		layers.push(Layer(vars, stopper));
		validate();
	}
	
	template <class Func>
    	void process_double(Func f)
    {
	    assert(!in_array_mode());
		check_type(DOUBLE, "DOUBLE");
		f(&current()->Val.V_DOUBLE, current());
		increment();
    }
    
	template <class Func>
    	void process_int(Func f)
    {
	    assert(!in_array_mode());
		check_type(INT, "INT");
		f(&current()->Val.V_INT, current());
		increment();
    }
    
	template <class Func>
    	void process_long(Func f)
    {
	    assert(!in_array_mode());
		check_type(LONG, "LONG");
		f(&current()->Val.V_LONG, current());
		increment();
    }

    template <class Func>    
    	void process_enum(Func f)
	{
	    assert(!in_array_mode());
		check_type(ENUM, "ENUM");
		f(&current()->Val.V_ENUM.value, current());
		
	    PremiaEnumMember * em = lookup_premia_enum(current(), current()->Val.V_ENUM.value);
	    if (em->nvar > 0)
	        push(em->Par, em->nvar);
		else
		    increment();
	}

	template <class Func>
    	void process_filename(Func f)
    {
	    assert(!in_array_mode());
		check_type(FILENAME, "FILENAME");
		f(&current()->Val.V_FILENAME, current());
		increment();
    }
    
    template <class Func>
        void process_array_begin(Func f)
    {
	    assert(!in_array_mode());
      assert_in_valid_position();
      if (current()->Vtype == PNLVECT || current()->Vtype == PNLVECTCOMPACT)
        {
          // start collecting values for the PNLVECTOR
          array_writer_.reset(new ArrayIterator(current(), f));
        }
      else
        {
          throw exception("current VAR is expected to be PNLVECT or PNLVECTCOMPACT");
        }				
    }

	
	bool in_array_mode() const { return array_writer_; }
	
	template <class Func>
	    void process_array_double(Func f)
	{
	    assert(in_array_mode());
	    array_writer_->process_double(f);
	}   
	
	template <class Func>
	    void process_array_end(Func f)
	{
	    assert(in_array_mode());
	    assert(array_writer_->finished());

        f(current());
        	    
		// exiting from array writing mode
		array_writer_.reset(0);
		increment();
	}
	
	//! \brief checks if all values are written into the VARs
	bool finished()  
	{
		return layers.empty();
	}

private:

	/// \brief assures that we are in position where we can write a VAR
	void assert_in_valid_position()
	{
		// we are supposed to be in valid position
		// otherwise an exception is thrown
		if (finished())
			throw premia::api::exception("VAR range_error");
	}

	//! \brief checks that current VAR accepts value of type 'type'
	//! \param label if the type is wrong an error with this label is issued
	void check_type(int type, std::string const & label)
	{
		assert_in_valid_position();
		if (true_typeV[current()->Vtype] != type)
			throw premia::api::exception("Expected data of " + label + " type");
	}

	//! \brief recurses into a compound VAR (i.e. NUM_FUNC_1)
	void push(VAR * vars, int stop = MAX_PAR)
	{
		--stopper();
		++current();

		// parameters of a NUM_FUNC are PREMIA_NULLTYPE terminated
		layers.push(Layer(vars, stop));

		validate();
	}

	//! \brief makes iterator valid
	void validate()
	{
		// ?
		if (layers.empty())
			return;

		// if current level is full, pop it and go to the previous one
		if (stopper() == 0)
		{
			pop(); validate();
			return;
		}

		VAR * c = current();

		// if current VAR is compound, recurse into its parameters introducing a layer
		switch (c->Vtype)
		{	
		case NUMFUNC_1:  push(c->Val.V_NUMFUNC_1->Par); return;
		case NUMFUNC_2:  push(c->Val.V_NUMFUNC_2->Par); return;
		case NUMFUNC_ND: push(c->Val.V_NUMFUNC_ND->Par); return;
		
		// PREMIA_NULLTYPE marks the end of a VAR sequence
		case PREMIA_NULLTYPE: pop(); validate(); return;
		}

		// if the current variable is unsetable bypass it 
		// decrementing the stopper since stopper takes them into account
		if (c->Vsetable == UNSETABLE) 
		{ 
			++current(); 
			--stopper();
			validate(); 
		}
	}


private:
	/// \brief represents VARs to be filled of one level (NUMFUNC_x introduce new levels)
	struct Layer {
		// pointer to variable to be written in
		VAR	* current;
		// how many variables need to be written in this layer
		// so, the range to write is [current, current + stopper)
		int   stopper;

		Layer(VAR * current, int stopper) 
			: current(current), stopper(stopper) 
		{}
	};

	/// \brief if in array writing mode, this variable points to the writer into the array
	boost::scoped_ptr<ArrayIterator>	array_writer_;

	/// \brief a VAR we are writing to
	VAR *& current() { return layers.top().current; }

	/// \brief how many elements left to be written in the current layer 
	int	 & stopper() { return layers.top().stopper; }

	/// \brief pops up the last layer
	void   pop() { layers.pop(); }

	/// \brief moves forward the iterator and validates it
	void   increment() 
	{
		++current(); --stopper();
		validate();
	}

	/// \brief the stack of layers
	std::stack<Layer>	layers;
};


}}}

#endif
