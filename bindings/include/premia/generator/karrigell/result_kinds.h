#ifndef _premia_generator_karrigell_result_kinds_h_included
#define _premia_generator_karrigell_result_kinds_h_included

#include <boost/assign/list_of.hpp>

namespace premia {
namespace pygen {
namespace karrigell {

   struct ResultKinds
   {
      ResultKinds() : next_id_(0)
      {}
      
      ResultKinds& operator() (bool show_by_default, std::list<const char*> const &elems)
      {
         BOOST_FOREACH(const char* e, elems)
         {
            result_to_kind_[e] = Kind(show_by_default, next_id_);
         }
         ++next_id_;
         return *this;
      }
      
      typedef unsigned                    KindId;
      typedef std::pair<bool, KindId>     Kind;
      typedef std::map<std::string, Kind> Result2Kind;

      Kind lookup(std::string const &res_name) const
      {
         Result2Kind::const_iterator it = result_to_kind_.find(res_name);
         if (it == result_to_kind_.end())
            throw std::logic_error(("cannot find kind for result " + res_name).c_str());
         return it->second;   
      }      

      static void printKind(Formatter &out, std::pair<std::string, Kind> const &k)
      {
         out("LABEL", k.first)("ID", k.second.second)("SHOW", k.second.first ? "True" : "False")
            << "'%LABEL%' : (%SHOW%, %ID%),";
      }
     
     template <class Ctx> 
	     friend Ctx& operator << (Ctx & ctx, ResultKinds const & kinds)
	   {
	      Formatter out(ctx.packageDir() / "result_kinds.py");   
	      
	      out << (seq, "result_kinds = {", +foreach_x(kinds.result_to_kind_, printKind), "}");
	      
	      return ctx;
	   }
      
   private:      
      Result2Kind    result_to_kind_;
      KindId         next_id_;
   };
   
   struct ResultKindsInitialized : ResultKinds
   {
      ResultKindsInitialized()
      {
         using boost::assign::list_of;
         
         (*this)
         (true, 
            list_of("Price")
                   ("OptimalPrice")
                   ("Price_Forward")
                   ("Price_Backward")
                   ("Forward_Price")
                   ("Backward_Price")
                   ("Price_by_2nd_Strategy")
                   ("Price_by_1st_Strategy")
                   ("Price_of_chosen_state")
                   ("Inf_Price")
                   ("Sup_Price")
                   ("Lower_Price")
                   ("Price_")
                   ("At_the_Money_Price")
                   ("At_the_Strike_Price")
         )
         (true,
            list_of("Inf_Delta1")
                   ("Inf_Delta2")
                   ("Delta")
                   ("Deltas")
                   ("Delta_")
                   ("Delta1")
                   ("Delta2")
                   ("Delta_of_chosen_state")
                   ("Lower_Delta")
                   ("Inf_Delta")
                   ("OptimalDelta")
                   ("Sup_Delta")
                   ("Sup_Delta2")
                   ("Sup_Delta1")
         )
         (true,
            list_of("Price_bp_")
         )
         (true, 
            list_of("Delta_bp_")
         )
         (true, 
            list_of("Asymptotics_for_Implied_Volatility_")
         )
         (true, 
            list_of("D_leg")
         )
         (true,
            list_of("P_leg")
         )
         (true,
            list_of("Price_in_annual_volatility_points")
         )
         (true,
            list_of("Price_in_10000_variance_points")
         )
         (false, 
            list_of("PriceError")
                   ("Price_Error")
                   ("MC_Error")
                   ("Error_Indicator")
                   ("ErrorPrice")
                   ("Error_Price")
                   ("Error")
         )
         (false,
            list_of("Error_Delta2")
                   ("Error_Delta1")
                   ("Error_Delta")
                   ("ErrorDelta")
                   ("DeltaError")
                   ("Delta_Error")
         )
         (true,
            list_of("Value_At_Risk")   
         )
         (true,
            list_of("Variance")
         )
         (true, 
            list_of("HedgeNow")
         )
         (true,
            list_of("Implied_Volatility_for_Small_Time")
         )
         (true,
            list_of("Fair_strike_value_in_annual_volatility_points")
                   ("Fair_strike_in_annual_volatility_points")
         )
         (true, 
            list_of("CDS_Spread")
         )
         (true,
            list_of("CDS_Spread_StdDev")
         )
         (true,
            list_of("Default_Leg")
         )
         (true,
            list_of("Premium_Leg")
         )
         (true, 
            list_of("Strikes")
                   ("Fair_strike_for_variance_swap")
         )
         (true,
            list_of("Strikes_Weights")
         )
         (true,
            list_of("Conditional_Tail_Expectation_")
         )         
         ;
      }
   };
 
}}}


#endif
