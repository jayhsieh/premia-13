#ifndef BOOST_RANGE_ADAPTOR_ZIPPED_PREPROCESSOR_DETAIL_HPP
#define BOOST_RANGE_ADAPTOR_ZIPPED_PREPROCESSOR_DETAIL_HPP

#ifndef BOOST_RANGE_MAX_ZIP_ARGUMENTS
#define BOOST_RANGE_MAX_ZIP_ARGUMENTS 5
#endif

#define BOOST_RANGE_MIN_ZIP_ARGUMENTS 2

#include <boost/config.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/preprocessor/arithmetic/dec.hpp>
#include <boost/preprocessor/arithmetic/div.hpp>
#include <boost/preprocessor/arithmetic/mul.hpp>
#include <boost/preprocessor/control.hpp>
#include <boost/preprocessor/control/while.hpp>
#include <boost/preprocessor/facilities/empty.hpp>
#include <boost/preprocessor/facilities/identity.hpp>
#include <boost/preprocessor/iteration/local.hpp>
#include <boost/preprocessor/punctuation/comma.hpp>
#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/range/adaptor/argument_fwd.hpp>
#include <boost/range/iterator_range.hpp>
#include <boost/type_traits/remove_reference.hpp>
#include <boost/utility/result_of.hpp>

namespace boost { namespace range_detail
{
	template< class ZI >
	struct zipped_range : public boost::iterator_range<boost::zip_iterator<ZI> >
	{
	private:
		typedef boost::iterator_range<boost::zip_iterator<ZI> > base;

	public:
		zipped_range( const ZI& beginIt, const ZI& endIt)
			: base( make_zip_iterator( beginIt ), make_zip_iterator( endIt ) ) {}
	};

	template <typename F, typename T, int SIZE>
	struct DetermineResultImpl;

	template <typename F, typename T>
	struct DetermineResult
		: DetermineResultImpl<F, T, tuples::length<T>::value> {};

	/////////////////////////////////////////////////////////////////////
	// DetermineResultImpl class specializations.
	/////////////////////////////////////////////////////////////////////
#define BOOST_RANGE_ZIPPED_ELEMENT_PRINT(z, n, data)                      \
	BOOST_DEDUCED_TYPENAME tuples::element<n, T>::type

#define BOOST_RANGE_ZIPPED_DETERMINE_RESULT(z, n, data)                   \
	template <typename F, typename T> struct DetermineResultImpl <F,T,n>    \
	: boost::result_of<F(                                                 \
	BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_ELEMENT_PRINT, ~))> {};

#define BOOST_PP_LOCAL_MACRO(n) BOOST_RANGE_ZIPPED_DETERMINE_RESULT(~,n,~)

#define BOOST_PP_LOCAL_LIMITS (BOOST_RANGE_MIN_ZIP_ARGUMENTS, \
	BOOST_RANGE_MAX_ZIP_ARGUMENTS)
#include BOOST_PP_LOCAL_ITERATE()

	/////////////////////////////////////////////////////////////////////
	// unpack_ function overloads.
	/////////////////////////////////////////////////////////////////////
#define BOOST_RANGE_ZIPPED_GET_PRINT(z, n, data)                          \
	get<n>(tuple)

#define BOOST_RANGE_ZIPPED_UNPACK_PRINT(z, n, data)                       \
	template <typename F, typename T> inline                                \
	BOOST_DEDUCED_TYPENAME DetermineResult<F,T>::type                       \
	unpack_(mpl::int_<n>, F f, const T& tuple)                            \
	{ return f(BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_GET_PRINT, ~)); }

#define BOOST_PP_LOCAL_MACRO(n) BOOST_RANGE_ZIPPED_UNPACK_PRINT(~,n,~)
#define BOOST_PP_LOCAL_LIMITS (BOOST_RANGE_MIN_ZIP_ARGUMENTS, \
	BOOST_RANGE_MAX_ZIP_ARGUMENTS)
#include BOOST_PP_LOCAL_ITERATE()

}//end range_detail

namespace adaptors
{
	/////////////////////////////////////////////////////////////////////
	// zip function overloads.
	/////////////////////////////////////////////////////////////////////

#define BOOST_RANGE_ZIPPED_SEQ(z, n, data)                                \
	boost::data(BOOST_PP_CAT(r,n))

#ifdef BOOST_NO_RVALUE_REFERENCES
	//We need overloads of zip for each combination of const and reference
	// arguments -- a huge explosion of code as n increases.

#define BOOST_RANGE_ZIPPED_EXP_PRED(d, data) BOOST_PP_TUPLE_ELEM(3, 0, data)

#define BOOST_RANGE_ZIPPED_EXP_OP(d, data) \
	( \
	BOOST_PP_DEC( \
	BOOST_PP_TUPLE_ELEM(3, 0, data) \
	), \
	BOOST_PP_TUPLE_ELEM(3, 1, data), \
	BOOST_PP_MUL_D( \
	d, \
	BOOST_PP_TUPLE_ELEM(3, 2, data), \
	BOOST_PP_TUPLE_ELEM(3, 1, data) \
	) \
	) \
	/**/

	// raise 'x' to the 'n'-th power -- example from pp documentation.
#define BOOST_RANGE_ZIPPED_EXP(x, n)                                      \
	BOOST_PP_TUPLE_ELEM(3, 2,                                               \
	BOOST_PP_WHILE(BOOST_RANGE_ZIPPED_EXP_PRED,                           \
	BOOST_RANGE_ZIPPED_EXP_OP, (n, x, 1)))

#define BOOST_RANGE_ZIPPED_BITSET_PRED(n, state)                          \
	BOOST_PP_TUPLE_ELEM(2,1,state)

#define BOOST_RANGE_ZIPPED_BITSET_OP(d, state)                            \
	(BOOST_PP_DIV_D(d, BOOST_PP_TUPLE_ELEM(2,0,state), 2),                  \
	BOOST_PP_DEC(BOOST_PP_TUPLE_ELEM(2,1,state)))

#define BOOST_RANGE_ZIPPED_BITSET(i, n)                                   \
	BOOST_PP_MOD(BOOST_PP_TUPLE_ELEM(2, 0,                                  \
	BOOST_PP_WHILE(BOOST_RANGE_ZIPPED_BITSET_PRED,                    \
	BOOST_RANGE_ZIPPED_BITSET_OP, (i,n))), 2)

#define BOOST_RANGE_ZIPPED_RANGE_ITERATOR(z, n, i)                        \
	BOOST_DEDUCED_TYPENAME range_iterator<BOOST_PP_CAT(R,n)                 \
	BOOST_PP_IF(BOOST_RANGE_ZIPPED_BITSET(i,n), BOOST_PP_IDENTITY(const), \
	BOOST_PP_EMPTY)()>::type

#define BOOST_RANGE_ZIPPED_ARGUMENTS(z, n, i)                             \
	BOOST_PP_CAT(R, n)                                                      \
	BOOST_PP_IF(BOOST_RANGE_ZIPPED_BITSET(i,n), const&, &)  \
	BOOST_PP_CAT(r, n)

#define BOOST_RANGE_ZIPPED_ZIP_PRINT_IMPL(z, i, n)                        \
	template <BOOST_PP_ENUM_PARAMS(n, typename R)>                          \
	inline range_detail::zipped_range<boost::tuple<                         \
	BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_RANGE_ITERATOR, i)> >             \
	zip(BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_ARGUMENTS, i))                  \
	{                                                                       \
	typedef boost::tuple<                                                 \
	BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_RANGE_ITERATOR, i)> RngTuple;   \
	return range_detail::zipped_range<RngTuple>(                          \
	RngTuple(BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_SEQ, begin)),          \
	RngTuple(BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_SEQ, end)));           \
}


#define BOOST_RANGE_ZIPPED_ZIP_PRINT(z, n, data)                          \
	BOOST_PP_REPEAT(BOOST_RANGE_ZIPPED_EXP(2,n),                            \
	BOOST_RANGE_ZIPPED_ZIP_PRINT_IMPL, n)
#else //we have rvalue references and don't need the 2^n function overloads.

	//TODO: We should also have a variadic template version -- one declaration.

#define BOOST_RANGE_ZIPPED_ARGUMENTS(z, n, i)                             \
	BOOST_PP_CAT(R, n)&& BOOST_PP_CAT(r, n)

#define BOOST_RANGE_ZIPPED_RANGE_ITERATOR(z, n, i)                        \
	BOOST_DEDUCED_TYPENAME range_iterator<                                  \
	BOOST_DEDUCED_TYPENAME boost::remove_reference<BOOST_PP_CAT(R,n)>::type \
	>::type 


#define BOOST_RANGE_ZIPPED_ZIP_PRINT(z, n, data)                          \
	template <BOOST_PP_ENUM_PARAMS(n, typename R)>                          \
	inline range_detail::zipped_range<boost::tuple<                         \
	BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_RANGE_ITERATOR, ~)> >             \
	zip(BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_ARGUMENTS, ~))                  \
	{                                                                       \
	typedef boost::tuple<                                                 \
	BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_RANGE_ITERATOR, ~)> RngTuple;      \
	return range_detail::zipped_range<RngTuple>(                          \
	RngTuple(BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_SEQ, begin)),          \
	RngTuple(BOOST_PP_ENUM(n, BOOST_RANGE_ZIPPED_SEQ, end)));           \
}

#endif

#define BOOST_PP_LOCAL_MACRO(n) BOOST_RANGE_ZIPPED_ZIP_PRINT(~,n,~)
#define BOOST_PP_LOCAL_LIMITS (BOOST_RANGE_MIN_ZIP_ARGUMENTS, \
	BOOST_RANGE_MAX_ZIP_ARGUMENTS)
#include BOOST_PP_LOCAL_ITERATE()

}}//end boost::adaptors

#endif

