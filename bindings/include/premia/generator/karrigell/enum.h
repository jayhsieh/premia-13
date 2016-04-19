#pragma once

#include <premia/generator/utils/formatter_dsl.h>
#include <premia/generator/api/enum.h>
#include <premia/generator/karrigell/ctx.h>
#include <premia/generator/karrigell/var.h>

namespace premia {
namespace pygen {
namespace karrigell {
	
	namespace print 
	{
	}
	
		inline void printEnumChoices(Formatter &out, api::Enum::Members::const_reference p)
		{
		    out("LABEL", p.second.quoted_original_label) <<  "'%LABEL%',";    
		}
	
	/// generates a wrapper for a model and a file with enumeration of all pricings for the model
	inline void generateEnum (Ctx const & ctx, Enum const & e)
	{
	    bool has_params = false;
	    
	    BOOST_FOREACH(api::Enum::Members::const_reference p, e.members)
	    { 
	        if (!p.second.params.empty())
	            has_params = true;
	    }

		Formatter(ctx.filename(e))
		   ("CHANGE", has_params ? "True" : "False")
		   ("NAME", e.label)
		   << (seq, 
		      call(print::commonHeader),	
		     "from kspremia.field_base import *",
		     "class %NAME%(FieldBase):", +(seq, 	
		     	  "def __init__(self, propertyName, friendlyName, fullName):",
		     	  "	super(%NAME%, self).__init__(propertyName, friendlyName, fullName)",
		     	  "",   
		          "labels = [", 
		            +foreach_x(e.members, printEnumChoices), 
		          "]",
		     	  "",	
		     	  "should_be_reloaded = %CHANGE%",
		     	  "",
                  call(boost::bind(print::enumInits, _1, boost::cref(e))),
                  "",
                  call(boost::bind(print::enumChoicesEx, _1, boost::cref(e))),
                  "",
				  "def process(self, pv):", +(seq,
					"member = getattr(pv.entity, self.propertyName)",
					"return pv.enumVisitor()(pv, member, self).processX()"
		   )));
	}

	inline Ctx& operator << (Ctx & ctx, Enums const & enums)
	{
		ctx.create(ctx.enumDir());
		
		BOOST_FOREACH(Enums::EnumsToLabels::const_reference r, enums.enums())
		{
		    generateEnum(ctx, r.second);
		}

		return ctx;
	}

}}}
