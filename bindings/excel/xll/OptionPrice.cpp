#include "stdafx.h"

#include "OptionPrice.h"

std::string Ref(const CellMatrix & referencee)
{
    return "#REF";
}

CellMatrix // finds an option price
OptionPrice(const std::string & requested_field,  // requested field
            const std::string & model_id, // model instance
            const std::string & option_id, // option instance
            const std::string & method_id  // pricing method instance
            )
{
    Excel_2::IOptionPricerPtr  pricer(__uuidof(Excel_2::COptionPricer));

    pricer->Init(requested_field.c_str(), model_id.c_str(), option_id.c_str(), method_id.c_str());

    long n_1, n_2;

    pricer->GetExtents(&n_1, &n_2);

    CellMatrix  res(n_1, n_2);

    for (long i_1 = 0; i_1 != n_1; ++i_1)
        for (long i_2 = 0; i_2 != n_2; ++i_2)
            res(i_1,i_2) = pricer->Run(i_1,i_2);

    return res;
}

CellMatrix // finds an option price in 1d iteration
OptionPrice_1(const std::string & requested_field,  // requested field
              const std::string & model_id, // model instance
              const std::string & option_id, // option instance
              const std::string & method_id,  // pricing method instance
              int it_1    // iteration index
              )
{
    Excel_2::IOptionPricerPtr  pricer(__uuidof(Excel_2::COptionPricer));

    pricer->Init(requested_field.c_str(), model_id.c_str(), option_id.c_str(), method_id.c_str());

    return CellMatrix(pricer->Run(it_1, 0));
}

CellMatrix // finds an option price in 2d iteration
OptionPrice_2(const std::string & requested_field,  // requested field
              const std::string & model_id, // model instance
              const std::string & option_id, // option instance
              const std::string & method_id,  // pricing method instance
              int it_1,    // iteration_1 index
              int it_2     // iteration_2 index
              )
{
    Excel_2::IOptionPricerPtr  pricer(__uuidof(Excel_2::COptionPricer));

    pricer->Init(requested_field.c_str(), model_id.c_str(), option_id.c_str(), method_id.c_str());

    return CellMatrix(pricer->Run(it_1, it_2));
}


std::string // returns ID corresponding to the region
PremiaRegionName(const std::string & name, // name to be used
                 const CellMatrix & region // region to be dependent upon
                 )
{
    return name;
}
