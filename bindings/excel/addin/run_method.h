#pragma once

#include "InputBuffer.h"
#include "formatting.h"
#include "readingcells.h"

namespace prxl
{
    inline void runMethod(db::id<db::Result> idx, fmt::DiffHandler * show_errors = 0)
    {
        // retrieving the problem's description object
        db::Entry<db::Result> const & res = db::theResults().at(idx);

        if (show_errors) show_errors->problemBegin(idx);

        // then, using it, retrieve objects describing the problem's parts: a model, an option and a method
        db::Entry<Model> const & mod = db::theModels().at(res.getModelId());
        db::Entry<Option> const & opt = db::theOptions().at(res.getOptionId());
        db::Entry<PricingMethod> const & met = db::theMethods().at(res.getMethodId());

        // getting raw Premia structures
        lib::Mode  mode = res.getMode();
        Model * model = mod.getData();
        Option * option = opt.getData();
        PricingMethod * method = met.getData();
        Pricing * pricing = lib::getPricing(mode, method);

        // iterators will hold information about parameters to be iterated over
        fmt::iterators_t iters;

        // number of rows in the result's region
        int result_height = res.getHeight();

        // top-left corner of the result's region
        Excel::RangePtr result_topleft = res.getTopLeft();

        // trying to load the model parameters, the option parameters and the method parameters
        if (lib::init(model) && fmt::load(mod, iters) 
         && lib::init(option, model) && fmt::load(opt, iters) 
         && lib::init(method, model) && fmt::load(met, iters))
        {
            // how many iterable parameters do we have?
            switch (iters.size())
            {
            case 0:
                {
                    if (lib::executeMethod(model, option, method, pricing))
                        fmt::formatResult(method, Down(Right(result_topleft)), show_errors);
                }
                break;
            case 1:
                {
                    // obtain place where the result should be written
                    Excel::RangePtr output_cols = fmt::getOutputBase(result_topleft);

                    if (!output_cols)
                        return;

                    // for each result field
                    for (int res_idx = 0; res_idx != result_height - 1; ++res_idx)
                    {
                        // print its name to the corresponding column header
                        fmt::printInBoldCell(
                            (_bstr_t)result_topleft->GetOffset(res_idx + 1, 0)->Value2, 
                            output_cols->GetOffset(0, res_idx),
                            env::getCellBkColor(stResult), env::getCellFontColor(stResult));

                        // and put a reference to it into the result domain cells
                        fmt::putRef(result_topleft->GetOffset(res_idx + 1, 1), output_cols->GetOffset(1, res_idx));
                    }

                    // performing 1d iteration
                    fmt::IterableParam &iter = iters.front();

                    // the first cell of the input array
                    Excel::RangePtr input_array = fmt::derefCell(Right(iter.getReferenceToKey()));

                    // iterating over the input array until an empty cell
                    for (int iteration_no = 0; input_array->Value2 != _variant_t(); 
                        ++iteration_no, input_array = Down(input_array))
                    {
                        // read the current cell and update the underlying Premia structure
                        iter.update((_bstr_t)input_array->Value2);

                        if (show_errors) show_errors->setIdx(iteration_no);

                        if (lib::executeMethod(model, option, method, pricing))
                            fmt::formatResult(method, output_cols->GetOffset(1 + iteration_no, 0), show_errors, 1, 0);
                    }
                }
                break;
            case 2:
                {
                    // the top-left cell of a region for outputting the results
                    Excel::RangePtr output_cols = fmt::getOutputBase(result_topleft);

                    if (!output_cols)
                        return;

                    // the first iterable parameter
                    fmt::IterableParam  &iter_1 = iters.front();
                    // the first cell of its input range
                    Excel::RangePtr input_1 = fmt::derefCell(Right(iter_1.getReferenceToKey()));
                    // the number of values in the iteration
                    int input_1_size = fmt::getIterationCount(input_1);

                    // the second iterable parameter
                    fmt::IterableParam &iter_2 = *++iters.begin();
                    // the first cell of its input range
                    Excel::RangePtr input_2 = fmt::derefCell(Right(iter_2.getReferenceToKey()));

                    // row difference between two consecutive output matrices  
                    int stride = input_1_size + 3;

                    // for each result field
                    for (int res_idx = 0; res_idx != result_height - 1; ++res_idx)
                    {
                        // calculate its base row of the output matrix
                        int outp_row = res_idx * stride;

                        // print the name of field into the top-left cell
                        fmt::printInBoldCell(
                            (_bstr_t)result_topleft->GetOffset(res_idx + 1, 0)->Value2, 
                            output_cols->GetOffset(outp_row, 0),
                            env::getCellBkColor(stResult), env::getCellFontColor(stResult));

                        // put reference to the output matrix in the result domain cell
                        fmt::putRef(result_topleft->GetOffset(res_idx + 1, 1), output_cols->GetOffset(1 + outp_row, 1));

                        // print labels with input values being iterated to the sides of the output array
                        fmt::printIterationLabels(input_1->GetOffset(0,0), output_cols->GetOffset(1 + outp_row, 0), 0, 1);
                        fmt::printIterationLabels(input_2->GetOffset(0,0), output_cols->GetOffset(0 + outp_row, 1), 1, 0);
                    }

                    // iterating over the first parameter
                    for (int iteration_no_1 = 0; input_1->Value2 != _variant_t(); 
                        ++iteration_no_1, input_1 = Down(input_1))
                    {
                        // updating the first parameter underlying structure
                        iter_1.update((_bstr_t)input_1->Value2);

                        // the first cell of the second parameter input array
                        Excel::RangePtr input = input_2->GetOffset(0,0);

                        // iterating over the second parameter
                        for (int iteration_no_2 = 0; input->Value2 != _variant_t(); 
                            ++iteration_no_2, input = Down(input))
                        {
                            // updating the second parameter underlying structure
                            iter_2.update((_bstr_t)input->Value2);

                            if (show_errors)
                                show_errors->setIdx(iteration_no_1, iteration_no_2);

                            if (lib::executeMethod(model, option, method, pricing))
                                fmt::formatResult(method, 
                                output_cols->GetOffset(1 + iteration_no_1, 1 + iteration_no_2), show_errors, 0, stride);
                        }
                    }
                }
                break;
            default:
                ::MessageBox(0, L"Iteration over than 2 parameters is not supported", L"Premia", MB_OK);
            }
        }

		result_topleft->Worksheet->Columns->AutoFit();

        //topleft->Worksheet->Columns->AutoFit();

        if (show_errors)
            show_errors->problemEnd();
    }
}