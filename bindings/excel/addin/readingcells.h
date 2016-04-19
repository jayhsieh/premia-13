#pragma once

#include "InputBuffer.h"
#include "formatting.h"
#include "db.h"
#include "execute_method.h"
#include "zcb.h"

namespace prxl {


    namespace fmt {

    /// Put a parameter to an InoutBuffer.
    /// If a parameter is a RandomGenerator, converts it.
    /// Checks that the parameter value is a number.
    /// \param key - a parameter key
    /// \param value - a parameter value
    /// \param buffer - an input buffer for storing the parameter
    template <class T> inline bool readParametersValue(_bstr_t key, _bstr_t value, T * target)
    {
        if (!lib::CheckNumericParameter(key, value))
            return false;

        return lib::putParameter(target, key, value);
    }

    inline int lookupKeyInEnum(VAR * V, _bstr_t const & value)
    {
        PremiaEnum * e = V->Val.V_ENUM.members;
        
        for (PremiaEnumMember * em = e->members; em->label; em = reinterpret_cast<PremiaEnumMember*>(reinterpret_cast<char*>(em) + e->size))
        {
            if (value == _bstr_t(em->label))
                return em->key;
        }

        return -1;
    }

    /// Loads data from a row into a double array.
    /// \param C - the first cell of the row. An empty cell is a mark for the end of the row.
    /// \param arr - an array to be loaded
    inline void loadPnlVect(Excel::RangePtr C, PnlVect * arr)
    {
        int cnt = 0;

        for (Excel::RangePtr I = C; I->Value2 != variant_t(); I = Right(I))
        {
            ++cnt;                        
        }

        pnl_vect_resize_from_double(arr, cnt, 0);

        int idx = 0;

        for (Excel::RangePtr I = C; I->Value2 != variant_t(); I = Right(I))
        {
            arr->array[idx++] = I->Value2;
        }
    }

    /// Represents a parameter to be iterated over.
    /// Stores a reference to the parameter's key and holds a functor updating an underlying premia structure.
    struct IterableParam
    {
        /// Constructor.
        /// \param T - underlying premia structure type
        /// \param reference_to_key - a reference to the parameter's key
        /// \param target - an underlying premia structure (Option, Model or PricingMethod)
        template <class T>
        IterableParam(Excel::RangePtr refernce_to_key, T * target)
            :   reference_to_key_(refernce_to_key)
            ,   action_(boost::bind(&lib::putParameter<T>, target, (_bstr_t)refernce_to_key->Value2, _1))
        {}

        /// \return reference to the parameter's key
        Excel::RangePtr getReferenceToKey() const { return reference_to_key_; }

        /// Writes new value of the parameter to the premia structure.
        /// \param v - a new parameter value
        void update(_bstr_t v)
        {
            action_(v);
        }

    private:
        /// reference to parameter's key
        Excel::RangePtr     reference_to_key_;
        /// a functor updating an underlying premia structure 
        boost::function<void (_bstr_t)> action_;
    };

    /// stores references to parameter key cell which values should be iterated over
    typedef std::list<IterableParam>  iterators_t;


/// reads an entity from its region.
/// \param T -- entity type
/// \param topleft - the top-left corner of an entity region
/// \param height - number of rows in the entity region
/// \param target - an entity to be filled with data from the region
/// \param iter - if there are iterable parameters they will be put into this variable
/// \return true if there are no errors in the parameters
template <class T>
bool readRegion(Excel::RangePtr topleft, int height, T * target, Excel::RangePtr enumRegion, iterators_t & iters)
{
    // iterate all region's row with parameters
    for (Excel::RangePtr C = Down(topleft); --height; C = Down(C))
    {
        Excel::RangePtr R = Right(C);

        // current parameter key
        _bstr_t key   = C->Value2;
        // current parameter value
        _bstr_t value = R->Value2;

        VAR * v = lib::findParameter(target, key);

        // check whether the parameter is ZCB Prices
        if (zcb::isZcbParameterInCell(key))
        {
            // if so, read the ZCB prices array
            zcb::readZcbRegion(derefCell(R));
        }
        else
        {
            // if the cell contains a reference, it means that the parameter is to be iterated over or it is a double array
            if (isRef(R))
            {
                if (v && v->Vtype == PNLVECT) 
                {
                    loadPnlVect(fmt::derefCell(R), v->Val.V_PNLVECT);
                }
				else if (v && v->Vtype == PNLVECTCOMPACT)
				{
					PnlVect *temp = pnl_vect_create(1);
					loadPnlVect(fmt::derefCell(R), temp);

					v->Val.V_PNLVECTCOMPACT->convert = 'a';
					v->Val.V_PNLVECTCOMPACT->size = temp->size;
					v->Val.V_PNLVECTCOMPACT->array = (double*)malloc(sizeof (double) * temp->size); 

					for (int i = 0; i != temp->size; ++i)
					{
						v->Val.V_PNLVECTCOMPACT->array[i] = temp->array[i];
					}

					pnl_vect_free(&temp);
				}
                else
                {
                    // NB! In case of iteration we do not fill target with all parameters
                    // leaving fields to be iterated uninitialized
                    iters.push_back(IterableParam(C, target));
                }
            }
            else
            {
                if (v && v->Vtype == ENUM)
                {
                    int e_idx = lookupKeyInEnum(v, value);

                    lib::putParameter(target, key, toString(e_idx));

                    if (lib::isEnumWithParameters(*v))
                    {
                        PremiaEnum * e = v->Val.V_ENUM.members;

                        for (PremiaEnumMember * em = e->members; em->label; 
                            em = reinterpret_cast<PremiaEnumMember*>(reinterpret_cast<char*>(em) + e->size))
                        {
                            enumRegion = Down(enumRegion);

                            for (int i = 0; i != em->nvar; ++i)
                            {
                                VAR * v = &em->Par[i];

                                // current parameter key
                                _bstr_t key   = enumRegion->Value2;

                                Excel::RangePtr R = Right(enumRegion);

                                // current parameter value
                                _bstr_t value = R->Value2;

                                if (em->key == e_idx)
                                {
                                    if (isRef(R))
                                    {
                                        if (v && v->Vtype == PNLVECT) 
                                        {
                                            loadPnlVect(fmt::derefCell(R), v->Val.V_PNLVECT);
                                        }
                                        else if (v && v->Vtype == PNLVECTCOMPACT)
                                        {
                                            PnlVect *temp = pnl_vect_create(1);
                                            loadPnlVect(fmt::derefCell(R), temp);

                                            v->Val.V_PNLVECTCOMPACT->convert = 'a';
                                            v->Val.V_PNLVECTCOMPACT->size = temp->size;
                                            v->Val.V_PNLVECTCOMPACT->array = (double*)malloc(sizeof (double) * temp->size); 

                                            for (int i = 0; i != temp->size; ++i)
                                            {
                                                v->Val.V_PNLVECTCOMPACT->array[i] = temp->array[i];
                                            }

                                            pnl_vect_free(&temp);
                                        }
                                        else
                                        {
                                            // We don't support iterable parameters in enums for the moment

                                            // NB! In case of iteration we do not fill target with all parameters
                                            // leaving fields to be iterated uninitialized
                                            //iters.push_back(IterableParam(C, target));
                                        }
                                    }
                                    else
                                    {
                                        if (v && v->Vtype == FILENAME)
                                            ;
                                        else
                                            if (lib::loadParameter(*v, value)) 
                                                return true;
                                            else
                                                ::MessageBox(0, L"Unable to load value = {" + value + L"} to parameter '" + key + L"'", L"Premia 10", MB_OK);
                                    }
                                }

                                enumRegion = Down(enumRegion);
                            }
                        }
                        enumRegion = Down(enumRegion);
                    }
                }
                else if (v && v->Vtype == FILENAME)
                    ;
                else
                    readParametersValue(key, value, target);
            }
        }
    }

    return true;
}

/// Returns the first cell of a region where the result of iteration will be put
/// \param result_base - the topleft corner of a result region
/// \return the first cell of a region where the result of iteration will be put
inline Excel::RangePtr getOutputBase(Excel::RangePtr result_base)
{
    // a reference to output region should be put in a top-right cell of a result region
    Excel::RangePtr R = Right(result_base);

    // if the cell doesn't contain a reference
    if (!isRef(R))
    {
        // issue an error
        ::MessageBox(0, address_of(R) + 
            " should contain reference to an output region for the result", L"Premia", MB_OK);

        return 0;
    }

    // otherwise return cell being referenced
    return derefCell(R);
}

/// Prints iteration labels (values of input arrays) to an output array.
/// \param inp_iter - a reference to the first cell of values array for an iterable parameter.
/// \param out_iter - a reference to the first cell where the output will be put
/// \param inc_x - column increment for two consecutive output cells 
/// \param inc_y - row increment for two consecutive output cells 
inline void printIterationLabels(Excel::RangePtr inp_iter, Excel::RangePtr out_iter, int inc_x, int inc_y)
{
    for (;inp_iter->Value2 != _variant_t(); inp_iter = Down(inp_iter), out_iter = out_iter->GetOffset(inc_y,inc_x))
    {
        Cell(out_iter).bold().bk_color(env::getCellBkColor(stResult)) = inp_iter->Value2;
    }
}

/// Calculates number of elements in an iteration
/// \param C - the first cell of input values array of a iterable parameter
/// \return number of element in the array
inline int getIterationCount(Excel::RangePtr C)
{
    int cnt = 0;

    for (; C->Value2 != _variant_t(); C = Down(C))
        ++cnt;

    return cnt;
}

/// Checks that the top-left cell of an entity region contains a string with an entity role label.
/// If not, issues an error.
/// \param T - entity type
/// \param e - an entity
/// \return true iff there are no errors
template <class T> 
    bool CheckEntityHeader(db::Entry<T> const & e)
{
    if ((_bstr_t)e.getTopLeft()->Value2 != (_bstr_t)(env::getLabel(role_of<T>())).c_str())
    {
        ::MessageBox(0, "Cell " + address_of(e.getTopLeft()) + 
            " should be equal to " + (env::getLabel(role_of<T>())).c_str(), L"Premia", MB_OK);
        return false;
    }

    return true;
}

/// Loads an entity from its region.
/// \param T - entity type
/// \param e - an entity
/// \param iters - output argument for storing information about iterable parameters 
/// \return true iff there are no errors
template <class T>
    bool load(db::Entry<T> const & e, iterators_t & iters)
{
    return
        CheckEntityHeader(e) &&
        readRegion(e.getTopLeft(), e.getHeight(), e.getData(), e.getEnumRegion(), iters) != 0;
}



}}