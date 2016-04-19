// RunTestsButton.h : déclaration de CRunTestsButton

#pragma once
#include "resource.h"       // symboles principaux



#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Les objets COM monothread ne sont pas correctement pris en charge par les plates-formes Windows CE, notamment les plates-formes Windows Mobile qui ne prennent pas totalement en charge DCOM. Définissez _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA pour forcer ATL à prendre en charge la création d'objets COM monothread et permettre l'utilisation de leurs implémentations. Le modèle de thread de votre fichier rgs a été défini sur 'Libre' car il s'agit du seul modèle de thread pris en charge par les plates-formes Windows CE non-DCOM."
#endif


namespace prxl {

// CRunTestsButton

[
	coclass,
	threading(apartment),
	vi_progid("Excel_2.RunTestsButton"),
	progid("Excel_2.RunTestsButton.1"),
	version(1.0),
	uuid("B71D11B6-7ECA-45C3-B505-62128BF1409B"),
	helpstring("RunTestsButton Class")
]
class ATL_NO_VTABLE CRunTestsButton :
    public IDispEventSimpleImpl<1,CRunTestsButton, &__uuidof(Office::_CommandBarButtonEvents)>
{
public:
	CRunTestsButton()
	{
	}

    void Advise(Office::_CommandBarButtonPtr btn)
    {
        DispEventAdvise((IDispatch*)(btn_ = btn),&__uuidof(Office::_CommandBarButtonEvents));
    }

    void UnAdvise()
    {
        //     DispEventUnadvise((IDispatch*)btn_);
        btn_ = 0;
    }

    /// This method is called when a user clicks the menu that show Premia wizard.
    /// It shows the wizard and after exiting it creates problem regions
    void __stdcall OnClickButton(IDispatch * /*Office::_CommandBarButton**/ Ctrl,  VARIANT_BOOL * CancelDefault)
    {
        USES_CONVERSION;
        CComQIPtr<Office::_CommandBarButton> pCommandBarButton(Ctrl);
        
        if (pCommandBarButton->Caption == _bstr_t("Run all tests"))
        {
            Excel::_WorksheetPtr sheet = g_app->ActiveWorkbook->GetWorksheets()->Add();

            Excel::RangePtr row = sheet->GetCells()->Get_Default(1,1);

            fmt::DiffHandler diff(row);

            for (db::Storage<db::Result>::iterator it(db::theResults()); it; ++it)
            {
                runMethod(*it, &diff);
            }
        }
    }

    BEGIN_SINK_MAP(CRunTestsButton)
        SINK_ENTRY_INFO(1, __uuidof(Office::_CommandBarButtonEvents),/*dispid*/ 0x01, OnClickButton, &OnClickButtonInfo)
    END_SINK_MAP()

	DECLARE_PROTECT_FINAL_CONSTRUCT()

	HRESULT FinalConstruct()
	{
		return S_OK;
	}

	void FinalRelease()
	{
	}

private:
    Office::_CommandBarButtonPtr btn_;
};

}