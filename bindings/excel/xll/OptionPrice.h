#pragma once

#include <xlw/CellMatrix.h>
 
//<xlw:libraryname=Premia
  
std::string  // tells to Premia that this parameter should be iterated over
Ref(CellMatrix const & referencee   // reference to the first cell of the iteration range
    );

CellMatrix // finds an option price
OptionPrice(const std::string & requested_field,  // requested field
            const std::string & model_id, // model instance
            const std::string & option_id, // option instance
            const std::string & method_id  // pricing method instance
            );

CellMatrix // finds an option price in 1d iteration
OptionPrice_1(const std::string & requested_field,  // requested field
              const std::string & model_id, // model instance
              const std::string & option_id, // option instance
              const std::string & method_id,  // pricing method instance
              int it_1    // iteration index
              );

CellMatrix // finds an option price in 2d iteration
OptionPrice_2(const std::string & requested_field,  // requested field
              const std::string & model_id, // model instance
              const std::string & option_id, // option instance
              const std::string & method_id,  // pricing method instance
              int it_1,    // iteration_1 index
              int it_2     // iteration_2 index
              );

std::string // returns ID corresponding to the region
PremiaRegionName(const std::string & name, // name to be used
                 const CellMatrix & region // region to be dependent upon
                 );
