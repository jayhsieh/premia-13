#pragma once

#include <boost/bind.hpp>

namespace premia {
namespace pygen {

	namespace formatter_dsl
	{

		/// \brief represents a list of lines to be formatted individually
		/// \param Prev type of the tail of the list
		/// \param T current line type
		template <class Prev, class T>
			struct Seq
		{
			/// constructs a list from its tail and head
			/// \param prev a reference to tail 
			/// \param cur a reference to head
			/// NB! We use references for performance reasons
			/// It implies that a sequence construction and formatting must be in a single expression
			Seq(Prev const & prev, T const & cur)
				:	prev_(prev), cur_(cur) 
			{}

			/// outputs the list of lines to a stream or formatter
			template <class Stream>
				friend Stream& operator << (Stream & out, Seq<Prev, T> const & x) 
			{
				return out << x.prev_ << x.cur_;
			}

		private:
			Prev const & prev_;		//!< list tail
			T    const & cur_;		//!< list head
		};

		/// \brief empty sequence of lines to be formatted
		template <> struct Seq<void,void>
		{
			/// does nothing
			template <class Stream>
				friend Stream& operator << (Stream & out, Seq<void, void> const & x) 
			{
				return out;
			}
		};

		/// a symbol to start constructing a list of lines to format
		Seq<void, void>  seq;

		/// constructs a new list of lines to format given an existing list and a new line
		/// \param prev previous lines (tail of the list)
		/// \param cur a line to add (head of the list)
		template <class Prev, class T1, class T2>
			Seq<Seq<Prev,T1>, T2>  operator , (Seq<Prev,T1> const & prev, T2 const & cur)
		{
			return Seq<Seq<Prev,T1>, T2>(prev, cur);
		}

		/// represents a block of lines of increased indentation level which may be surrounded by user defined parentheses
		template <class T>
			struct Block
		{
			/// constructs a block
			/// \param t lines in the block
			/// \param pars 2 symbol array with opening and closing symbols for the block
			explicit Block(T const& t, char const * pars = 0) : t_(t), pars_(pars) {}

			/// outputs the block into a formatter
			template <class Stream>
				friend Stream & operator << (Stream & out, Block const & x)
			{
				// putting opening symbol if exists
				if (x.pars_)	out.putch(x.pars_[0]);
				// increasing indent
				out.incindent();
				// putting block contents
				out << x.t_;
				// decreasing indent
				out.decindent();
				// putting closing symbol if exists
				if (x.pars_)	out.putch(x.pars_[1]);

				return out;
			}

		private:
			T const & t_;		//!< block contents
			char const * pars_;	//!< 2 symbols array with opening and closing symbols
		};

		/// creates a block without opening and closing symbols
		/// \param x block content
		template <class T> Block<T> block(T const & x) { return Block<T>(x); }

		/// creates a block with given opening and closing symbols
		/// \param pars 2 symbol array with opening and closing symbols for the block e.g. "{}" or "<>"
		/// \param x block content
		template <class T> Block<T> block(char const pars[2], T const & x) { return Block<T>(x, pars); }

		/// equivalent to block(x)
		template <class T1, class T2> Block<Seq<T1,T2> > operator + (Seq<T1,T2> const & x)
		{
			return Block<Seq<T1,T2> >(x);
		}

		/// equivalent to block(x)
		template <class T> Block<Block<T> > operator + (Block<T> const & x)
		{
			return Block<Block<T> >(x);
		}

		/// represents a line generator
		/// \param Range type fulfilling boost Range concept
		/// \param Func a functional object of signature void operator () (Stream &, Range::const_reference)
		template <class Range, class Func>
			struct ForEach
		{
			/// creates a generator
			/// \param range a range of source values
			/// \param func function generating lines from range's values 
			ForEach(Range const & range, Func func) : range_(range), func_(func) {}

			/// triggers line generation to the formatter
			template <class Stream>
			friend Stream & operator << (Stream & out, ForEach const & x) 
			{
				/// NB!!! We assume that each iteration have the same variables set
				out.push_scope();
                                
				std::for_each(x.range_.begin(), x.range_.end(), boost::bind(x.func_, boost::ref(out), _1));
				out.pop_scope();
				return out;
			}

		private:
			Range const & range_;	//!< a range of source values
			Func          func_;	//!< function generating lines from range's values 
		};

		/// creates a ForEach generator
		template <class Range, class Func> ForEach<Range,Func> foreach_x(Range const & rng, Func f)
		{
			return ForEach<Range,Func>(rng, f);
		}

		/// equivalent to block(x)
		template <class Range, class Func> Block<ForEach<Range,Func> > operator + (ForEach<Range,Func> const & x)
		{
			return Block<ForEach<Range,Func> >(x);
		}

		/// \brief represents a call to a line generator
		/// \param Func callable type Stream& -> void
		template <class Func> struct Call
		{
			/// creates a caller
			explicit Call(Func f) : func_(f) {}

			/// calls the line generator 
			template <class Stream>
				friend Stream& operator << (Stream & out, Call const & x)
				{
					out.push_scope();
					x.func_(out);
					out.pop_scope();
					return out;
				}

		private:
			Func    func_;	//!< callable object that generates lines to formatter
		};

		/// creates a line generator
		template <class Func> Call<Func> call(Func f) { return Call<Func>(f); }

		/// equivalent to block(x)
		template <class Func> Block<Call<Func> > operator + (Call<Func> const & x)
		{
			return Block<Call<Func> >(x);
		}
	}

	using formatter_dsl::block;
	using formatter_dsl::call;
	using formatter_dsl::foreach_x;
	using formatter_dsl::seq;
}}
