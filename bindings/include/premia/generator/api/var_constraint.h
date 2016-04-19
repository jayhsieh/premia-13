#pragma once

#include <numeric>
#include <premia/import.h>

#undef max
#undef min

namespace premia {
namespace pygen  {
namespace api	 {

	/// \brief represents a range constraint on a VAR of scalar type
	/// \param Scalar type of the scalar
	template <class Scalar>
		struct Range 
		{
			/// constructs the range
			/// \param low lower bound of the range. can be std::numeric_limits<Scalar>::min() to indicate that there is no lower limit
			/// \param low_inclusive true iff the range includes its lower bound
			/// \param hi upper bound of the range. can be std::numeric_limits<Scalar>::max() to indicate that there is no upper limit
			/// \param hi_inclusive true iff the range includes its upper bound
			Range(Scalar low, bool low_inclusive, Scalar hi, bool hi_inclusive)
				:	low(low)
				,	low_inclusive(low_inclusive)
				,	hi(hi)
				,	hi_inclusive(hi_inclusive)
			{}

			/// checks whether the range has lower bound
			bool has_low() const { return low != std::numeric_limits<Scalar>::min(); }
			/// checks whether the range has upper bound
			bool has_hi() const { return hi != std::numeric_limits<Scalar>::max(); }
			
			Scalar	low;			//!< lower bound
			bool	low_inclusive;	//!< true iff the range includes its lower bound
			Scalar	hi;				//!< upper bound
			bool	hi_inclusive;	//!< true iff the range includes its upper bound
		};


	/// \brief creates a wrapper for the VAR constraint if any
	/// \param Vtype VAR's type 
	/// \param IntType type of the value held by the VAR
	/// \return pointer to created constraint wrapper, 0 if there is no constraint
	template <class IntType>
		Range<IntType> const *	getRangeConstraint(int Vtype)
		{
			static IntType Max = std::numeric_limits<IntType>::max();

			static Range<IntType>	pint(1, true, Max, false);
			static Range<IntType>   int2(2, true, Max, false);
			static Range<IntType>	rgint130(1, true, 30, true);
			static Range<IntType>	rgint13	(1, true, 3, true);
			static Range<IntType>	rgint12	(1, true, 2, true);

			switch (Vtype)
			{
			case PINT:			return &pint;
			case INT2:			return &int2;
			case RGINT130:		return &rgint130;
			case RGINT13:		return &rgint13;
			case RGINT12:		return &rgint12;
			}

			return 0;
		}

	template <>
		Range<double> const * getRangeConstraint<double>(int Vtype)
	{
		static double Max = std::numeric_limits<double>::max();
		static double Min = std::numeric_limits<double>::min();

		static Range<double> date		(0., true, Max, false);
		static Range<double> pdouble	(0., true, Max, false);
		static Range<double> sndouble	(Min, false, 0., false);
		static Range<double> rgdouble	(0., true, 1., true);
		static Range<double> rgdouble1	(1., false, Max, false);
		static Range<double> rgdoublem11(-1., true, 1., true);
		static Range<double> rgdouble12	(1., true, 2., true);
		static Range<double> rgdouble02	(0., true, 2., true);
		static Range<double> sdouble2	(2., false, Max, false);
		static Range<double> spdouble	(0., true, Max, false);
		static Range<double> rgdouble051(0.5, true, 1., true);
		static Range<double> rgdouble14 (1., true, 4., true);

		switch (Vtype)
		{
		case DATE:			return &date;
		case PDOUBLE:		return &pdouble;
		case SNDOUBLE:		return &sndouble;
		case RGDOUBLE:		return &rgdouble;
		case RGDOUBLE1:		return &rgdouble1;
		case RGDOUBLEM11:	return &rgdoublem11;
		case RGDOUBLE12:	return &rgdouble12;
		case RGDOUBLE02:	return &rgdouble02;
		case SDOUBLE2:		return &sdouble2;
		case SPDOUBLE:		return &spdouble;
		case RGDOUBLE051:	return &rgdouble051;
		case RGDOUBLE14: 	return &rgdouble14;
		}

		return 0;
	}
}

using api::Range;

}}