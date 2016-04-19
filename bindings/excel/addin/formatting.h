#pragma once

#include "result_entry.h"
#include "helpfiles.h"
#include "diff_handler.h"

namespace prxl {
    inline Role role_of(db::Result *) { return stResult; }

    /// Contains utilities for formatting Excel cells 
    namespace fmt {
        namespace zcb {
            inline bool isZcbParameter(std::string *pKey);
        }

        /// Puts a reference into a cell.
        /// \param where - a cell where the reference should be put
        /// \param what - a cell a reference to which should be put to 'where'
        inline void putRef(Excel::RangePtr where, Excel::RangePtr what)
        {
            where->FormulaLocal = _bstr_t("=REF(") + address_of(what) + ")";
        }

        /// Checks whether the cell contains a reference
        /// \param C - a cell to be tested
        /// \return true iff the cell contains a reference
        inline bool isRef(Excel::RangePtr C)
        {
            _bstr_t v = C->FormulaLocal;
            return 0 == _wcsnicmp(v.GetBSTR(), L"=REF(", 5);
        }

        /// Provided the cell contains a reference this function dereferences it
        /// \param C - a cell with reference
        /// \return a cell the reference points to
        inline Excel::RangePtr derefCell(Excel::RangePtr C)
        {
            std::string s((_bstr_t)C->FormulaLocal);

            s = s.substr(5, s.size() - 6);

            return C->Worksheet->GetRange(s.c_str());
        }

        inline void addValidation(Excel::RangePtr C, _bstr_t choices, _bstr_t initial_value, _bstr_t error_msg)
        {
            // adding list constraint to the cell
            Excel::ValidationPtr v = C->GetValidation();

            v->Delete();
            v->Add(Excel::xlValidateList, Excel::xlValidAlertStop, Excel::xlBetween, choices);
            v->IgnoreBlank = VARIANT_TRUE;
            v->InCellDropdown = VARIANT_TRUE;
            v->ErrorTitle = "Premia";
            v->ErrorMessage = error_msg;
            v->ShowInput = VARIANT_TRUE;
            v->ShowError = VARIANT_TRUE;

            C->Value2 = initial_value;
        }

        inline void putValueIntoCell(Excel::RangePtr R, _variant_t value, DiffHandler * show_errors)
        {
            if (show_errors) {
                double a = R->Value;
                double b = value;

                double eps = 1e-7;

                if (abs(a - b) > eps)
                    show_errors->showDifference(R, value);
            }
            else
                R->Value2 = value;
        }

        inline _bstr_t getChoices(VAR const & V)
        {
            PremiaEnum * e = V.Val.V_ENUM.members;
            
            _bstr_t acc = "";

            for (PremiaEnumMember * em = e->members; em->label; 
                em = reinterpret_cast<PremiaEnumMember*>(reinterpret_cast<char*>(em) + e->size))
            {
                acc += em->label;
                acc += ";";
            }

            return acc;
        }

        inline const char * getCurrentChoice(VAR const & V)
        {
            return lookup_premia_enum(const_cast<VAR*>(&V), V.Val.V_ENUM.value)->label;
        }

        struct RegionAcc
        {
            RegionAcc() : low_col_(0xffff), low_row_(0xffff), hi_col_(0), hi_row_(0) {}
    
            void add(Excel::Range * pRange, int ncols = 1, int nrows = 1)
            {
                low_row_ = low_row_ < pRange->Row                ? low_row_ : pRange->Row;
                low_col_ = low_col_ < pRange->Column             ? low_col_ : pRange->Column;
                hi_row_  = hi_row_  > pRange->Row + nrows - 1    ? hi_row_  : pRange->Row + nrows - 1;
                hi_col_  = hi_col_  > pRange->Column + ncols - 1 ? hi_col_  : pRange->Column + ncols - 1;
            }

            Excel::RangePtr getTopLeft(Excel::_Worksheet * sheet) const
            {
                return sheet->GetCells()->Get_Default(low_row_,low_col_);
            }

            Excel::RangePtr getRightBottom(Excel::_Worksheet * sheet) const
            {
                return sheet->GetCells()->Get_Default(hi_row_, hi_col_);
            }

        private:
            long low_row_, hi_row_, low_col_, hi_col_;
        };

        inline void applyEnumStyle(Excel::RangePtr R, _bstr_t trigger)
        {
            R->Font->ColorIndex = 15;

            Excel::FormatConditionsPtr  FC = R->FormatConditions;

            FC->Add(Excel::xlExpression, Excel::xlEqual, trigger);
            FC->Item(FC->Count)->Font->ColorIndex = 1;
        }

        inline void putVarIntoCell(PnlVect const * V, Excel::RangePtr R, _bstr_t style = L"", DiffHandler * show_errors = 0)
		{
            putRef(R, Right(R));

            for (int i = 0; i != V->size; ++i)
            {
                if (show_errors) show_errors->setIdx(i);

                Excel::RangePtr C = R->GetOffset(0,i + 1);

                putValueIntoCell(C, V->array[i], show_errors);

                if (style.length() > 0)
                {
                    applyEnumStyle(C, style);
                }
            }
		}

		struct CompVect
		{
			CompVect(PnlVectCompact const * src)
				:	temp_(pnl_vect_compact_to_pnl_vect(src))
			{}

			operator PnlVect const * () const 
			{
				return temp_;
			}

			PnlVect const * operator -> () const
			{
				return temp_;
			}

			~CompVect()
			{
				pnl_vect_free(&temp_);
			}

		private:
			PnlVect * temp_;
		};


        /// Formats a VAR into an Excel cell
        /// \param V - a VAR
        /// \param R - an excel cell
        inline void putVarIntoCell(VAR const & V, Excel::RangePtr R, _bstr_t style = L"", DiffHandler * show_errors = 0)
        {
			if (V.Vtype < FIRSTLEVEL && V.Vtype != ENUM && V.Vtype != FILENAME)
			{
				switch (true_typeV[V.Vtype])
				{
				case INT:
					putValueIntoCell(R, V.Val.V_INT, show_errors); break;
				case 3 /*LONG*/:
					putValueIntoCell(R, V.Val.V_LONG, show_errors); break;
				case DOUBLE:
					putValueIntoCell(R, V.Val.V_DOUBLE, show_errors); break;
				default:
	                R->Value2 = "Unsupported type";
				}
			}
			else
			{
				switch (V.Vtype)
				{
            case FILENAME:
                putValueIntoCell(R, V.Val.V_FILENAME, show_errors); break;

            case PNLVECT:
                {
					putVarIntoCell(V.Val.V_PNLVECT, R, style, show_errors);
                }
                break;

			case PNLVECTCOMPACT:
				{
					CompVect temp(V.Val.V_PNLVECTCOMPACT);

					putVarIntoCell(temp, R, style, show_errors);
				}
				break;

            case ENUM:
                {
                    addValidation(R, getChoices(V), getCurrentChoice(V), "Only values from the drop-down box is allowed");
                    break;
                }

            default:
                R->Value2 = "Unsupported type";
				}
			}
        }

        /// A wrapper over Excel::RangePtr that provides a bit more convenient operations over a cell.
        struct Cell
        {
            /// Constructs Cell from Excel::RangePtr
            /// \param C - a RangePtr
            Cell(Excel::RangePtr C) : C(C) {}

            /// Sets cell background color
            /// \param color - Excel color index
            Cell & bk_color(int color)
            {
                C->Interior->ColorIndex = color;
                return *this;
            }

            /// Sets cell font color
            /// \param color - Excel color index
            Cell & font_color(int color)
            {
                C->Font->ColorIndex = color;
                return *this;
            }

            /// Makes the cell font bold
            Cell & bold() 
            {
                C->Font->Bold = VARIANT_TRUE;
                return *this;
            }

            /// Writes a value into the cell
            /// \param v - a value to be written
            Cell & operator = (_variant_t const & v)
            {
                C->Value2 = v;
                return *this;
            }

        private:
            /// the underlying RangePtr
            Excel::RangePtr C;
        };

        /// Colorizes a cell into header color for a given role
        /// \param r - a role 
        /// \param C - a cell to be colorized
        inline void colorizeHeaderCell(Role r, Cell C)
        {
            C.bold().bk_color(env::getHeaderBkColor(r)).font_color(env::getHeaderFontColor(r));
        }

        /// Colorizes a cell, makes it bold and put some value in it
        /// \param what - a value to put into the cell
        /// \param C - a cell
        /// \param bk_color - an excel color index for background color
        /// \param font_color - an excel color index for font color
        inline void printInBoldCell(const char * what, Cell C, int bk_color, int font_color = 0)
        {
            C.bold().bk_color(bk_color).font_color(font_color) = what;
        }

        /// Prints a parameter name into a cell
        /// \param r - role of an entity the parameter belongs to
        /// \param name - a parameter name
        /// \param C - a cell where to print the parameter name
        inline void printParameterName(Role r, const char * name, Excel::RangePtr C)
        {
            printInBoldCell(name, C, env::getCellBkColor(r), env::getCellFontColor(r));
        }

        /// Add a conditional style in a cell so the cell background gets cyan as soon as its value gets equal to "#REF"
        /// \param R - a cell
        inline void addRefConditionalStyle(Excel::RangePtr R)
        {
            Excel::FormatConditionsPtr  FC = R->FormatConditions;

            //FC->Delete();
            FC->Add(Excel::xlCellValue, Excel::xlEqual, "=\"#REF\"");
            FC->Item(FC->Count)->GetInterior()->ColorIndex = 34;
        }

        /// formats an entity header
        /// \param r - entity role
        /// \param name - entity label
        /// \param helpfilename - a path to helpfile, may be equal to 0
        /// \param subaddr - if the hyperlink points to a cell this parameter shouold contain the cell address otherwise it must be 0
        /// \param depends_on - an address of a range of entity parameters values
        /// \param C - topleft cell for the header. it is moved down at the function exit 
        inline void formatHeader( Role r, const char * name, const char * helpfilename, 
            const char * subaddr, const char * depends_on, Excel::RangePtr &C )
        {
            C->Value2 = env::getLabel(r).c_str();
            colorizeHeaderCell(r, C);

            Excel::RangePtr R = Right(C);

            if (name)
            {
                R->FormulaLocal = "=PREMIAREGIONNAME(\"" + _bstr_t(name) + "\";" + depends_on + ")";

                if (helpfilename)
                    R->Worksheet->Hyperlinks->Add(R, helpfilename, subaddr ? subaddr : vtMissing);
            }
            colorizeHeaderCell(r, R);

            C = Down(C);
        }
        /// returns address of the entity range as string
        /// \param topleft - top-left corner of an entity region
        /// \param height - the region height in rows
        /// \return address of the entity range as string
        inline _bstr_t entityRangeString(Excel::RangePtr topleft, int height)
        {
            return address_of(Down(topleft)) + (height == 1 ? "" : ":" + address_of(topleft->GetOffset(height - 1, 1)));
        }

        inline void getRegion(std::vector<VAR> const & pars, Excel::RangePtr topleft, RegionAcc * reg_acc, bool * has_zcb)
        {
            Excel::RangePtr C = topleft;

            reg_acc->add(C);

            for (std::vector<VAR>::const_iterator it = pars.begin(); it != pars.end(); ++it, C = Down(C))
            {
                if (it->Vtype == PNLVECT)
                    reg_acc->add(C, it->Val.V_PNLVECT->size + 2, 1);

                if (it->Vtype == PNLVECTCOMPACT)
				{
					CompVect temp(it->Val.V_PNLVECTCOMPACT);
                    reg_acc->add(C, temp->size + 2, 1);
				}


                std::string key = it->Vname;

                *has_zcb = *has_zcb || zcb::isZcbParameter(&key);
            }

            reg_acc->add(C, 2, 1);
        }

        inline void formatRegionEnumParameters(Role r, std::vector<VAR> const & pars, Excel::RangePtr& topleft, Excel::RangePtr original)
        {
            Excel::RangePtr &C = topleft;
            Excel::RangePtr O = Right(original);

            BOOST_FOREACH(VAR const &V, pars)
            {
                if (lib::isEnumWithParameters(V))
                {
                    PremiaEnum * e = V.Val.V_ENUM.members;

                    for (PremiaEnumMember * em = e->members; em->label; 
                        em = reinterpret_cast<PremiaEnumMember*>(reinterpret_cast<char*>(em) + e->size))
                    {
                        printParameterName(r, V.Vname, C);
                        Right(C)->Value2 = em->label;

                        _bstr_t trigger = "=("+ address_of(Right(C)) +"="+ address_of(O) +")";

                        applyEnumStyle(C, trigger);
                        applyEnumStyle(Right(C), trigger);

                        C = Down(C);

                        for (int i = 0; i != em->nvar; ++i)
                        {
                            VAR const & param = em->Par[i];

                            Excel::RangePtr R = Right(C);

                            std::string par_name = param.Vname;

                            printParameterName(r, par_name.c_str(), C);
                            putVarIntoCell(param, R, trigger);

                            applyEnumStyle(C, trigger);
                            applyEnumStyle(Right(C), trigger);

                            addRefConditionalStyle(R);
                            C = Down(C);
                        }
                    }
                    C = Down(C);
                }
                O = Down(O);
            }
        }

        /// Formats entity region.
        /// \param r - entity role
        /// \param name - entity label (can be 0)
        /// \param helpfilename - help file name (can be 0)
        /// \param pars - an array of entity parameters
        /// \param topleft - top-left corner of a region'
        /// \param m - mode which the entity belongs to
        /// \param pZcbCell - reference to cell which is initialized iff the entity has ZCB prices region
        /// \return height of the region in rows
        inline int formatRegion(Role r, const char * name, const char * helpfilename,  
            std::vector<VAR> const & pars, Excel::RangePtr topleft, lib::Mode m,
            Excel::RangePtr& pZcbCell)
        {
            Excel::RangePtr C = topleft;

            // height is parameters count + 1 row for the header
            int rows_number = (int)pars.size() + 1;

            // address of entity parameters range 
            _bstr_t field = entityRangeString(C, rows_number);

            formatHeader(r, name, helpfilename, 0, field, C);

            // for each entity parameter
            for (std::vector<VAR>::const_iterator it = pars.begin(); it != pars.end(); ++it, C = Down(C))
            {
                Excel::RangePtr R = Right(C);

                std::string par_name = it->Vname;

                bool isZcb = zcb::isZcbParameter(&par_name);

                printParameterName(r, par_name.c_str(), C);
                putVarIntoCell(*it, R);
                addRefConditionalStyle(R);

                if (isZcb) 
                {
                    C = Down(C);

                    printParameterName(r, "ZCB Prices Region", C);

                    addRefConditionalStyle(pZcbCell = Right(C));

                    ++rows_number;
                }
            }

            // put border around the region
            topleft->GetResize(rows_number, 2)->BorderAround(1, Excel::xlMedium, Excel::xlColorIndexAutomatic);

            return rows_number;
        }

        //inline void getRegion()

        /// Formats a result of a pricing method.
        /// \param x - a pricing method
        /// \param topleft - where should the first field be put
        /// \param inc_x - columns stride (support for 1d and 2d iteration)
        /// \param inc_y - row stride (support for 1d and 2d iteration)
        /// \return the last filled cell
        inline Excel::RangePtr formatResult(PricingMethod * x, Excel::RangePtr topleft, DiffHandler * show_errors, int inc_x = 0, int inc_y = 1)
        {
            int count = lib::getParArraySize(x->Res);

            for (VAR * v = x->Res; count; count--, v++, topleft = topleft->GetOffset(inc_y, inc_x))
            {
                if (show_errors) show_errors->setCurrentParameter(v->Vname);
                putVarIntoCell(*v, topleft, "", show_errors);
            }

            return topleft;
        }

        /// Formats a region for results of a problem and registers the problem in a db sheet.
        /// \param x - a pricing method
        /// \param m - a mode for the method
        /// \param topleft - top-left corner of a region where to format the result
        /// \param model_id - db sheet id of a model for the problem
        /// \param option_id - db sheet id of an option for the problem
        /// \param method_id - db sheet id of a pricing method for the problem
        /// \return db sheet id of the problem
        inline db::id<db::Result> formatResultRegion(PricingMethod * x, lib::Mode m, Excel::RangePtr &topleft,
            db::id<Model> model_id, db::id<Option> option_id, db::id<PricingMethod> method_id)
        {
            Excel::RangePtr dummy;
            lib::ResultOf<PricingMethod> r(x);
            // format the region
            int height = formatRegion(stResult, 0, 0, lib::getParameters(&r), topleft, m, dummy);
            // register the problem in a db sheet
            db::id<db::Result> id = db::theResults().push(db::Entry<db::Result>(topleft, height, m, model_id, option_id, method_id));
            topleft = topleft->GetOffset(height, 0);
            return id;
        }


        /// Formats an entity region and registers the entity in a db sheet.
        /// \param x - an entity
        /// \param m - an entity mode
        /// \param zcb_cell - if the entity contains ZCB prices region then this out-parameter will be set as reference to ZCB prices region
        /// \param topleft - a top-left corner of a region where to format the entity
        /// \return sb sheet id of the entity
        template <class T>
        db::id<T> formatRegion(T * x, lib::Mode m, Excel::RangePtr &zcb_cell, Excel::RangePtr &topleft)
        {
            // format the region
            int height = formatRegion(role_of(x), db::generateNewLabel<T>(x->Name), 
                hlp::getHelpFileName(m,x, std::cerr), lib::getParameters(x), topleft, m, zcb_cell);

            // register the entity in the db sheet
            db::id<T> id = db::theStorage<T>().push(db::Entry<T>(x->Name, topleft, height, m, x));
            topleft = topleft->GetOffset(height, 0);
            return id;
        }

        template <class T>
            void formatRegionEnumParameters(T * x, lib::Mode m, db::id<T> idx, Excel::RangePtr &topleft)
        {
            db::Entry<T> & e = db::theStorage<T>().at(idx);
            e.setEnumRegion(topleft);
            e.serialize(db::theStorage<T>().locate(idx));

            Excel::RangePtr original = e.getTopLeft();

            // format the region
            formatRegionEnumParameters(role_of(x), lib::getParameters(x), topleft, Down(original));
        }

        /// Formats header for an existing entity
        /// \param T - the entity type
        /// \param idx - the entity db sheet id
        /// \param topleft - top-left corner of a region where to format the entity
        /// \return the entity db sheet id
        template <class T>
        db::id<T> formatHeaderForExistingEntity(db::id<T> idx, Excel::RangePtr & topleft)
        {
            Role r = role_of((T*)0);

            db::Entry<T> const & e = db::theStorage<T>().at(idx);

            Excel::RangePtr pos = e.getTopLeft();
            int h = e.getHeight();

            _bstr_t field = entityRangeString(pos, h);

            formatHeader(r, _bstr_t("--> ") + (_bstr_t)Right(pos)->Value2, "", 
                _bstr_t(pos->Worksheet->Name) + "!" + address_of(pos), field, topleft);

            return idx;
        }

    }}