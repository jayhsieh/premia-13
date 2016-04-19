// AppEventsHandler.h : déclaration de CAppEventsHandler

#pragma once
#include "resource.h"       // symboles principaux



#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Les objets COM monothread ne sont pas correctement pris en charge par les plates-formes Windows CE, notamment les plates-formes Windows Mobile qui ne prennent pas totalement en charge DCOM. Définissez _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA pour forcer ATL à prendre en charge la création d'objets COM monothread et permettre l'utilisation de leurs implémentations. Le modèle de thread de votre fichier rgs a été défini sur 'Libre' car il s'agit du seul modèle de thread pris en charge par les plates-formes Windows CE non-DCOM."
#endif

namespace prxl {

//DEFINE_GUID ( DIID_AppEvents, 0x00024413,0x0000,0x0000, 0xc0,0x00,0x00,0x00,0x00,0x00,0x00,0x46 );
/// DIID of Excel::IAppEvents interface
extern "C" const GUID __declspec(selectany) DIID_AppEvents =
   {0x00024413,0x0000,0x0000,{0xc0,0x00,0x00,0x00,0x00,0x00,0x00,0x46}};

/// used to connect to listening Excel dispinterface events
_ATL_FUNC_INFO afiExcelAppDocumentOpen					= { CC_STDCALL, VT_EMPTY, 1, VT_DISPATCH };
#define DISPID_WORKBOOKOPEN			0x0000061f

static _ATL_FUNC_INFO afiExcelAppDocumentBeforeClose	= { CC_STDCALL, VT_EMPTY, 2, { VT_BOOL | VT_BYREF, VT_DISPATCH} };
#define DISPID_WORKBOOKBEFORECLOSE	0x00000622
    
/// Listener to excel's events about opening saved workbooks.
/// Since saved workbooks may contain Premia data 
/// we must find all Premia's Compute buttons on the sheets
/// and connect them to CButtonHandlers
[
	coclass,
	threading(apartment),
	vi_progid("Excel_2.AppEventsHandler"),
	progid("Excel_2.AppEventsHandler.1"),
	version(1.0),
	uuid("20F7A19F-FE41-4C2F-9C55-4DD14355E652"),
	helpstring("AppEventsHandler Class")
]
class ATL_NO_VTABLE CAppEventsHandler :
    public IDispEventSimpleImpl< 1, CAppEventsHandler, &DIID_AppEvents >
{
public:

    /// starts listening to Excel workbook events
    /// \param pUnk - pointer to Excel application
    void Advise(IUnknown * pUnk)
    {
        HRESULT hr = DispEventAdvise( pUnk );
    }

    BEGIN_SINK_MAP( CAppEventsHandler )
        SINK_ENTRY_INFO( 1, DIID_AppEvents, DISPID_WORKBOOKOPEN,			OnDocumentOpen,			&afiExcelAppDocumentOpen )
        SINK_ENTRY_INFO( 1, DIID_AppEvents, DISPID_WORKBOOKBEFORECLOSE,	    OnBeforeClose,			&afiExcelAppDocumentBeforeClose )
    END_SINK_MAP( )


	DECLARE_PROTECT_FINAL_CONSTRUCT()

private:

    void __stdcall OnBeforeClose(VARIANT_BOOL *, IDispatch * pDispDoc)
    {
        db::cleanUpStorage<Model>();
        db::cleanUpStorage<Option>();
        db::cleanUpStorage<PricingMethod>();
        db::cleanUpStorage<db::Result>();
    }

    /// Document opening event handler
    /// \param pDispDoc - IDispatch of a document being opened
    void __stdcall OnDocumentOpen			( IDispatch *pDispDoc )
    {
        Excel::_WorkbookPtr wkb = pDispDoc;

        Excel::SheetsPtr sheets = wkb->GetWorksheets();

        // iterate through all worksheets
        for (long idx = 1; idx <= sheets->GetCount(); ++idx)
        {
            Excel::_WorksheetPtr aSheet = sheets->GetItem(idx);

            Excel::OLEObjectsPtr oles = aSheet->OLEObjects();

            // iterate all OLE objects contained in the sheet
            for (long ole_no = 1; ole_no <= oles->GetCount(); ++ole_no)
            {
                Excel::_OLEObjectPtr obj = oles->Item(ole_no);

                // if it is button
                if (obj->GetprogID() == _bstr_t("Forms.CommandButton.1"))
                {
                    // check if it is our button
                    // by looking up in the hidden sheet its position 
                    for (db::Storage<db::Result>::iterator button_id(db::theResults()); button_id ; ++button_id)
                    {
                        Excel::RangePtr button_pos = db::getButtonPosition(*button_id);

                        if (button_pos->Worksheet->Name == aSheet->Name)
                        {
                            if (button_pos->Top == obj->TopLeftCell->Top)
                            {
                                if (button_pos->Left == obj->TopLeftCell->Left)
                                {
                                    // if it is our button let's start listening its events
                                    MSForms::ICommandButtonPtr button = obj->Object;

                                    com_object<CButtonHandler>()->Advise(button, *button_id);
                                }
                            }
                        }
                    }
                }
            }
        }
    }

};

}