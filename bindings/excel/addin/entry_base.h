#pragma once

#include "premia.h"
#include "excel.h"
#include "hidden_sheets.h"

namespace prxl {

/// Namespace for work with db hidden sheets
namespace db
{
    /// A cached representation of a entity instance registered in the db sheets.
    /// Fields of this class correspond to the db columns:
    /// -# Entity name
    /// -# A pointer to the top-left corner of its region
    /// -# Number of rows in its region
    /// -# A mode that the entity belongs to
    struct EntryBase
    {
        /// A constructor for just created entities
        /// \param n - Entity name
        /// \param r - A pointer to the top-left corner of the entity region
        /// \param h - Number of rows in the entity region
        /// \param m - A mode that the entity belongs to
        EntryBase(_bstr_t const & n, Excel::RangePtr r, int h, lib::Mode m)
            :   name_(n), topleft_(r), height_(h), mode_(m)
        {}

        /// A constructor that reads its fields from db sheet row.
        /// Is used when a saved document is being opened
        /// \param O - a row to be read
        EntryBase(Excel::RangePtr O)
        {
            name_ = At(O, fld::EntityName)->Value2;

            _bstr_t sheet_name = At(O, fld::SheetName)->Value2;

            Excel::_WorksheetPtr sheet = 
                excelApp()->GetActiveWorkbook()->GetWorksheets()->GetItem(sheet_name);  

            topleft_ =  sheet->GetRange(At(O, fld::Address)->Value2);

            height_ = At(O, fld::Height)->Value2;
            mode_ = lib::Mode((long)At(O, fld::Mode)->Value2);

            _bstr_t enumRegionString = At(O, fld::EnumAddress)->Value2;
            if (enumRegionString.length() > 0)
                enumRegion_ = sheet->GetRange(enumRegionString);
        }

        /// Flushes entity fields into a db sheet row
        /// \param row - a db row where fields of the entity should be written
        void serialize(Excel::RangePtr row) const
        {
            At(row, fld::EntityName)->Value2 = name_;
            At(row, fld::SheetName)->Value2 = topleft_->Worksheet->Name;
            At(row, fld::Address)->Value2 = address_of(topleft_);
            At(row, fld::Height)->Value2 = height_;
            At(row, fld::Mode)->Value2 = mode_;
            
            if (enumRegion_)
                At(row, fld::EnumAddress)->Value2 = address_of(enumRegion_);
        }

        /// returns a string describing entity type, e.g. BlackScholes1dim
        _bstr_t         getName()   const { return name_; }
        /// returns RangePtr to the top-left cell of the entity region
        Excel::RangePtr getTopLeft()const { return topleft_; }
        /// returns number of rows in the entity region
        int             getHeight() const { return height_; }
        /// returns mode corresponding to the entity
        lib::Mode       getMode()   const { return mode_; }

        Excel::RangePtr getEnumRegion() const { return enumRegion_; }
        void setEnumRegion(Excel::RangePtr topleft) { enumRegion_ = topleft; }
        
        /// returns entity instance label, e.g. BlackScholes1dim:2
        _bstr_t         getLabel() const
        {
            return (_bstr_t)Right(getTopLeft())->Value2;
        }

    private:
        /// Entity name
        _bstr_t         name_;
        /// A pointer to the top-left corner of the entity region
        Excel::RangePtr topleft_;
        /// Number of rows in the entity region
        int             height_;
        /// A mode that the entity belongs to
        lib::Mode       mode_;
        /// a pointer to entity enum parameters region
        Excel::RangePtr enumRegion_;
    };

    /// returns an entity instance label, e.g. BlackScholes1dim:2
    /// \param e - an entity whose label is retrieved
    inline _bstr_t getLabel(EntryBase const & e)
    {
        return e.getLabel();
    }

    template <class T>
    struct id 
    {
        explicit id(int no) : no_(no) {}

        int getValue() const { return no_; }

        static id notExist() { return id(-1); }

        friend bool operator == (id<T> a, id<T> b) { return a.getValue() == b.getValue(); }

    private:
        int no_;
    };

}}