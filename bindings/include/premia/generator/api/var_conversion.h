#ifndef _premia_generator_api_var_conversion_h_included_
#define _premia_generator_api_var_conversion_h_included_

#include <iostream>
#include <premia/generator/api/var_wrapper.h>
#include <premia/generator/api/enum.h>

namespace premia {
namespace pygen {
namespace api   {

	/// converts a PnlVect into std::vector<double>
	std::vector<double>	 asVector(PnlVect const * v)
	{
		return std::vector<double>(v->array, v->array + v->size);
	}

	/// converts a PnlVectCompact into std::vector<double>
	std::vector<double>	 asVector(PnlVectCompact const * vc)
	{
		std::vector<double> res(vc->size);

		if (vc->convert == 'd') 
			std::fill(res.begin(), res.end(), vc->val);
		else 
			std::copy(vc->array, vc->array + vc->size, &res[0]);
		return res;
	}

	/// creates a wrapper for a VAR
	/// \param ctx a reference to the transformation context
	/// \param v source VAR
	template <class Ctx>
	Var convertToVar(Ctx & ctx, VAR const & v)
	{
		switch (true_typeV[v.Vtype])	
		{
		case ENUM:	
		    {
		        Enum const *e = &ctx.enums.insert(ctx, v.Val.V_ENUM.members);
		        return EnumValue(v.Val.V_ENUM.value, e);
		    }
		case INT:	return numeric(v.Val.V_INT, v.Vtype);
		case 3:	    return numeric(v.Val.V_LONG, v.Vtype);
		case DOUBLE:return numeric(v.Val.V_DOUBLE, v.Vtype);
		case BOOL:	throw std::logic_error("BOOL type is not supported as input parameter type");
		case FILENAME:	return  std::string(v.Val.V_FILENAME);
		case PNLVECT:	return	asVector(v.Val.V_PNLVECT);
		case PNLVECTCOMPACT:	return asVector(v.Val.V_PNLVECTCOMPACT);
		}
		throw std::logic_error("unexpected VAR type for " + std::string(v.Vname));
	}

	/// collects an entity members as a list of NamedVars
	/// \param ctx a reference to the transformation context
	/// \param vars a pointer to the array of the entity members
	/// \param stopper maximal number of members to consider
	/// \param used_pars a dictionary of used identifiers
	/// \param var_list output list of NamedVars
	template <class Ctx>
	int getVars(Ctx & ctx, VAR const* vars, int stopper, UsedIds & used_pars, VarList & var_list)
	{
		for (; stopper && vars->Vtype != PREMIA_NULLTYPE; ++vars)
		{
			switch (vars->Vtype)
			{
			case NUMFUNC_1:
				getVars(ctx, vars->Val.V_NUMFUNC_1->Par, MAX_PAR, used_pars, var_list);
				--stopper;
				break;

			case NUMFUNC_2:
				getVars(ctx, vars->Val.V_NUMFUNC_2->Par, MAX_PAR, used_pars, var_list);
				--stopper;
				break;

			case NUMFUNC_ND:
				getVars(ctx, vars->Val.V_NUMFUNC_ND->Par, MAX_PAR, used_pars, var_list);
				--stopper;
				break;

			default:
				if (vars->Vsetable == SETABLE)
				{
				    bool iterable = 
				        vars->Viter == ALLOW &&
				        vars->Vtype != ENUM &&
				        vars->Vtype != FILENAME;
				        
					var_list.push_back(
						NamedVar(vars, 
								correctIdentifierName(vars->Vname, used_pars), 
								convertToVar(ctx, *vars), 
								vars->setter != 0,
								iterable));
				}
				--stopper;
			}
		}

		return stopper;
	}

	template <class Ctx>
	int getVars(Ctx & ctx, VAR * vars, int stopper, VarList & var_list)
	{
		UsedIds used_ids;
		return getVars(ctx, vars, stopper, used_ids, var_list);
	}

}}}

#endif
