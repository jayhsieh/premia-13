#pragma once

namespace prxl {
namespace fmt {

    struct DiffHandler
    {
        void problemBegin(db::id<db::Result> idx)
        {
            id_ = idx;
            has_errors_ = false;

            idx_1_ = idx_2_ = -1;

            (base_ = row_)->Value2 = idx.getValue();

            row_->Worksheet->Hyperlinks->Add(row_, "", full_address_of(db::theResults().at(idx).getTopLeft()));

            row_ = Down(row_);
        }

        void problemEnd()
        {
            base_->Interior->ColorIndex = has_errors_ ? 3 : 6;
        }

        void setCurrentParameter(const char * parameter_name)
        {
            current_parameter_ = parameter_name;
        }

        void setIdx(int idx_1){ idx_1_ = idx_1; }
        void setIdx(int idx_1, int idx_2) { idx_2_ = idx_2; idx_1_ = idx_1; }

        void showDifference(Excel::RangePtr original, _variant_t const & newValue)
        {
            has_errors_ = true;

            _bstr_t name = current_parameter_;

            if (idx_1_ != -1)
                ((name += "[") += toString(idx_1_)) += "]";

            if (idx_2_ != -1)
                ((name += "[") += toString(idx_2_)) += "]";

            row_->GetOffset(0,1)->Value2 = name;
            row_->Worksheet->Hyperlinks->Add(row_->GetOffset(0,1), "", full_address_of(original));

            row_->GetOffset(0,2)->Value2 = original->Value2;
            row_->Worksheet->Hyperlinks->Add(row_->GetOffset(0,2), "", full_address_of(original));

            row_->GetOffset(0,3)->Value2 = newValue;

            row_ = Down(row_);
        }

        DiffHandler(Excel::RangePtr row) : row_(row), id_(db::id<db::Result>::notExist()) 
        {
            row->Value2 = "Problem #";
            row->GetOffset(0,1)->Value2 = "Parameter:";
            row->GetOffset(0,2)->Value2 = "Old value:";
            row->GetOffset(0,3)->Value2 = "New value:";

            row_ = Down(row_);
        }

        ~DiffHandler()
        {
            row_->Worksheet->Columns->AutoFit();
        }

    private:
        Excel::RangePtr     row_, base_;
        _bstr_t             current_parameter_;
        db::id<db::Result>  id_;
        bool                has_errors_;
        int                 idx_1_, idx_2_;
    };


}
}