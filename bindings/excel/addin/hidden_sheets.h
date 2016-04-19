#pragma once

#include "settings.h"

namespace prxl {
namespace db
{
    Excel::RangePtr getBaseCell(Role r, Excel::_WorkbookPtr wkb);

    /// returns an Excel sheet corresponding to the given role.
    /// \param wkb - Excel workbook where the sheet is searched
    /// \param r - sheet role 
    /// \return a pointer to the worksheet
    inline Excel::_WorksheetPtr  getSheet(Excel::_WorkbookPtr wkb, Role r)
    {
        Excel::SheetsPtr sheets = wkb->GetWorksheets();

        // iterate through all worksheets in the workbook
        for (long idx = 1; idx <= sheets->GetCount(); ++idx)
        {
            Excel::_WorksheetPtr aSheet = sheets->GetItem(idx);

            // if the sheet has the special name return reference to it
            if (aSheet->Name == _bstr_t(env::getSheetName(r).c_str()))
                return aSheet;
        }

        // otherwise, create a new sheet

        // remember currently active sheet
        Excel::_WorksheetPtr temp = wkb->ActiveSheet;

        // add a new sheet
        Excel::_WorksheetPtr aSheet = sheets->Add();

		aSheet->PutVisible(0, Excel::xlSheetHidden);		

        // and restore the old one
        temp->Activate();

        aSheet->Name = env::getSheetName(r).c_str();

        // put to the topleft cell number of db rows in the sheet
        Excel::RangePtr(aSheet->GetCells()->Get_Default(1,1))->Value2 = 0L;

        return aSheet;
    }

    /// returns the top-left cell of db sheet for the given role
    /// \param r - sheet role 
    /// \param wkb - workbook
    /// \return the top-left cell of the db sheet
    inline Excel::RangePtr getBaseCell(Role r, Excel::_WorkbookPtr wkb)
    {
        return getSheet(wkb, r)->GetCells()->Get_Default(1,1);
    }

    /// returns the leftmost cell of 'idx'-th row in the db sheet for the role 'r'
    /// \param r - sheet role
    /// \param wkb - workbook
    /// \param idx - row index
    /// \return the leftmost cell of the row
    inline Excel::RangePtr    getRow(Role r, Excel::_WorkbookPtr wkb, int idx)
    {
        Excel::RangePtr O = getBaseCell(r, wkb);

        return O->GetOffset(idx, 0);
    }

    /// returns the leftmost cell of 'idx'-th row in the db sheet for the role 'r'
    /// \param r - sheet role
    /// \param idx - row index
    /// \return the leftmost cell of the row
    inline Excel::RangePtr getRow(Role r, int idx)
    {
        return getRow(r, excelApp()->GetActiveWorkbook(), idx);
    }

    /// returns number of rows in the db sheet corresponding to the role 'r'
    /// \param r - sheet role
    /// \param wkb - workbook
    /// \return number of rows 
    inline int getRowCount(Role r, Excel::_WorkbookPtr wkb)
    {
        return getBaseCell(r, wkb)->Value2;
    }

    /// Contains a db sheet fields enumeration 
    namespace fld
    {
        /// enumerates db sheet column identifiers
        enum Offset
        {
            /// name of the entity type, e.g. BlackScholes1dim
            EntityName = 1,     
            /// sheet name of the entity region
            SheetName,          
            /// address of the top-left cell of the entity region
            Address,            
            /// number of rows in the entity region
            Height,             
            /// mode for the entity
            Mode,               
            /// model id of the problem (only for the results db)
            ModelId,            
            /// option id of the problem (only for the results db)
            OptionId,           
            /// method id of the problem (only for the results db)
            MethodId,
            /// address of the top-left cell of the entity enum parameters region
            EnumAddress
        };
    }

    /// returns a reference to the cell given by field ID in the db row 'O'
    /// \param O - a db row
    /// \param field - a column ID
    /// \return a pointer to the cell
    inline Excel::RangePtr At(Excel::RangePtr O, fld::Offset field)
    {
        return O->GetOffset(0, (int)field);
    }

    /// return number of problems registered in db
    inline int getButtonsCount()
    {
        return getRow(stResult, 0)->Value2;
    }

}}