//// 
//// created by xlwplus
////

#include <xlw/pragmas.h>
#include <xlw/MyContainers.h>
#include <xlw/CellMatrix.h>
#include "optionprice.h"
#include <xlw/xlw.h>
#include <xlw/XlFunctionRegistration.h>
#include <stdexcept>
#include <xlw/XlOpenClose.h>
#include <ctime>
namespace {
const char* LibraryName = "Premia";
};

// dummy function to force linkage
namespace {
void DummyFunction()
{
xlAutoOpen();
xlAutoClose();
}
}

// registrations start here


namespace
{
XLRegistration::Arg
RefArgs[]=
{
{ "referencee"," reference to the first cell of the iteration range "}
};
  XLRegistration::XLFunctionRegistrationHelper
registerRef("xlRef",
"Ref",
" tells to Premia that this parameter should be iterated over ",
LibraryName,
RefArgs,
"P"
);
}



extern "C"
{
LPXLOPER EXCEL_EXPORT
xlRef(
LPXLOPER referenceea)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper referenceeb(
	(referenceea));
CellMatrix referencee(
	referenceeb.AsCellMatrix("referencee"));

std::string result(
	Ref(
		referencee)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
OptionPriceArgs[]=
{
{ "requested_field"," requested field "},
{ "model_id"," model instance "},
{ "option_id"," option instance "},
{ "method_id"," pricing method instance "}
};
  XLRegistration::XLFunctionRegistrationHelper
registerOptionPrice("xlOptionPrice",
"OptionPrice",
" finds an option price ",
LibraryName,
OptionPriceArgs,
"RRRR"
);
}



extern "C"
{
LPXLOPER EXCEL_EXPORT
xlOptionPrice(
LPXLOPER requested_fielda,
LPXLOPER model_ida,
LPXLOPER option_ida,
LPXLOPER method_ida)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper requested_fieldb(
	(requested_fielda));
std::string requested_field(
	requested_fieldb.AsString("requested_field"));

XlfOper model_idb(
	(model_ida));
std::string model_id(
	model_idb.AsString("model_id"));

XlfOper option_idb(
	(option_ida));
std::string option_id(
	option_idb.AsString("option_id"));

XlfOper method_idb(
	(method_ida));
std::string method_id(
	method_idb.AsString("method_id"));

CellMatrix result(
	OptionPrice(
		requested_field,
		model_id,
		option_id,
		method_id)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
OptionPrice_1Args[]=
{
{ "requested_field"," requested field "},
{ "model_id"," model instance "},
{ "option_id"," option instance "},
{ "method_id"," pricing method instance "},
{ "it_1"," iteration index "}
};
  XLRegistration::XLFunctionRegistrationHelper
registerOptionPrice_1("xlOptionPrice_1",
"OptionPrice_1",
" finds an option price in 1d iteration ",
LibraryName,
OptionPrice_1Args,
"RRRRB"
);
}



extern "C"
{
LPXLOPER EXCEL_EXPORT
xlOptionPrice_1(
LPXLOPER requested_fielda,
LPXLOPER model_ida,
LPXLOPER option_ida,
LPXLOPER method_ida,
double it_1a)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper requested_fieldb(
	(requested_fielda));
std::string requested_field(
	requested_fieldb.AsString("requested_field"));

XlfOper model_idb(
	(model_ida));
std::string model_id(
	model_idb.AsString("model_id"));

XlfOper option_idb(
	(option_ida));
std::string option_id(
	option_idb.AsString("option_id"));

XlfOper method_idb(
	(method_ida));
std::string method_id(
	method_idb.AsString("method_id"));

int it_1(
	static_cast<int>(it_1a));

CellMatrix result(
	OptionPrice_1(
		requested_field,
		model_id,
		option_id,
		method_id,
		it_1)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
OptionPrice_2Args[]=
{
{ "requested_field"," requested field "},
{ "model_id"," model instance "},
{ "option_id"," option instance "},
{ "method_id"," pricing method instance "},
{ "it_1"," iteration_1 index "},
{ "it_2"," iteration_2 index "}
};
  XLRegistration::XLFunctionRegistrationHelper
registerOptionPrice_2("xlOptionPrice_2",
"OptionPrice_2",
" finds an option price in 2d iteration ",
LibraryName,
OptionPrice_2Args,
"RRRRBB"
);
}



extern "C"
{
LPXLOPER EXCEL_EXPORT
xlOptionPrice_2(
LPXLOPER requested_fielda,
LPXLOPER model_ida,
LPXLOPER option_ida,
LPXLOPER method_ida,
double it_1a,
double it_2a)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper requested_fieldb(
	(requested_fielda));
std::string requested_field(
	requested_fieldb.AsString("requested_field"));

XlfOper model_idb(
	(model_ida));
std::string model_id(
	model_idb.AsString("model_id"));

XlfOper option_idb(
	(option_ida));
std::string option_id(
	option_idb.AsString("option_id"));

XlfOper method_idb(
	(method_ida));
std::string method_id(
	method_idb.AsString("method_id"));

int it_1(
	static_cast<int>(it_1a));

int it_2(
	static_cast<int>(it_2a));

CellMatrix result(
	OptionPrice_2(
		requested_field,
		model_id,
		option_id,
		method_id,
		it_1,
		it_2)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

namespace
{
XLRegistration::Arg
PremiaRegionNameArgs[]=
{
{ "name"," name to be used "},
{ "region"," region to be dependent upon "}
};
  XLRegistration::XLFunctionRegistrationHelper
registerPremiaRegionName("xlPremiaRegionName",
"PremiaRegionName",
" returns ID corresponding to the region ",
LibraryName,
PremiaRegionNameArgs,
"RP"
);
}



extern "C"
{
LPXLOPER EXCEL_EXPORT
xlPremiaRegionName(
LPXLOPER namea,
LPXLOPER regiona)
{
EXCEL_BEGIN;

	if (XlfExcel::Instance().IsCalledByFuncWiz())
		return XlfOper(true);

XlfOper nameb(
	(namea));
std::string name(
	nameb.AsString("name"));

XlfOper regionb(
	(regiona));
CellMatrix region(
	regionb.AsCellMatrix("region"));

std::string result(
	PremiaRegionName(
		name,
		region)
	);
return XlfOper(result);
EXCEL_END
}
}



//////////////////////////

