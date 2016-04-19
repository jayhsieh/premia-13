// PremiaAddin.h : déclaration de CPremiaAddin

#pragma once
#include "resource.h"       // symboles principaux


#include "dialog.h"
#include "db.h"
#include "formatting.h"
#include "zcb.h"
#include "ButtonHandler.h"
#include "AppEventsHandler.h"
#include "nd.h"
#include "run_method.h"
#undef DATE

#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Les objets COM monothread ne sont pas correctement pris en charge par les plates-formes Windows CE, notamment les plates-formes Windows Mobile qui ne prennent pas totalement en charge DCOM. Définissez _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA pour forcer ATL à prendre en charge la création d'objets COM monothread et permettre l'utilisation de leurs implémentations. Le modèle de thread de votre fichier rgs a été défini sur 'Libre' car il s'agit du seul modèle de thread pris en charge par les plates-formes Windows CE non-DCOM."
#endif

/// event callback metainformation
_ATL_FUNC_INFO OnClickButtonInfo = {CC_STDCALL,VT_EMPTY,2,{VT_DISPATCH,VT_BYREF | VT_BOOL}};

#include "WizardButtonHandler.h"
#include "RunTestsButton.h"

/// Namespace for Premia Excel interface names
namespace prxl {



/// The main coclass of the add-in. 
/// the instance of the class is created by Excel when it starts.
/// ProgId of the class is written down as a subbranch to HKCU\Software\Microsoft\Office\Excel\Addins registry branch
/// This class does following things:
/// 1. At the add-in loading it adds a menu to the main Excel menu and begins listening the menu's events.
/// 2. When the menu is chosen it shows the wizard dialog for selecting an option, a model and a method
/// 3. If a user selects a model, an option and a method 
///    a) it fills a range in active worksheet with data describing the selected model, option and method
///    b) and adds 'Compute' button beneath the region
/// 4. At the add-in loading it connects to listening events about opening existing documents
[
    coclass,
    threading(apartment),
    vi_progid("Excel_2.PremiaAddin"),
    progid("Excel_2.PremiaAddin.1"),
    version(1.0),
    uuid("3FB75F7F-7EBE-4916-ADD0-9C9D459C2367"),
    helpstring("PremiaAddin Class")
]
class ATL_NO_VTABLE CPremiaAddin :
    public IDispatchImpl<_IDTExtensibility2, &__uuidof(_IDTExtensibility2), &LIBID_AddInDesignerObjects, /* wMajor = */ 1>
{
    // _IDTExtensibility2 Methods
public:

    static Excel::AddInPtr findAddin(Excel::_Application * spApp, _bstr_t const & addin_path)
    {
        Excel::AddInsPtr addins = spApp->AddIns;

        for (long i = 0; i != addins->GetCount(); ++i)
        {
            _bstr_t fullname = addins->GetItem(i + 1)->GetFullName();

            if (_wcsicmp(fullname, addin_path) == 0)
                return addins->GetItem(i + 1);
        }

        return addins->Add(addin_path);
    }

    /// the method is called by Excel when the add-in is being loaded.
    /// -# Registers XLL with Premia functions
    /// -# Adds a menu for starting Premia wizard
    STDMETHOD(OnConnection)(LPDISPATCH Application, ext_ConnectMode ConnectMode, LPDISPATCH AddInInst, SAFEARRAY * * custom)
    {
        InitVar();
        InitErrorMsg();

        // obtaining _Application interface from IDispatch
        CComQIPtr<Excel::_Application>  spApp(Application);

        // registering XLL with Premia functions
        // the XLL must be located in directory given by XLL registry entry
        spApp->RegisterXLL(env::getXLLPath() + "\\premia.xll", 0);

        // storing a pointer to Excel instance in the global variable
        g_app = spApp;

        addPremiaMenu(spApp);

        // registering XLA with old Premia8 interface
        findAddin(spApp, env::getXLAPath() + "\\premia.xla")->PutInstalled(VARIANT_TRUE);

        // switching the current directory to the Unix directory of Premia installation
        // TODO: Do this only when initialyield.dat file is touched
        // TODO: Switch back the current directory to the initial one
        SetCurrentDirectory(env::getDataPath());

        strcpy(premia_data_dir, env::getDataPath());
        path_sep = "\\";

        /// start listening events about a document opening 
        /// in order to connect to existing there Compute buttons
        com_object<CAppEventsHandler>()->Advise(spApp);

        return S_OK;
    }




    /// This method is called when the add-in is being disconnected.
    STDMETHOD(OnDisconnection)(ext_DisconnectMode RemoveMode, SAFEARRAY * * custom)
    {
        wizardMenu_->UnAdvise();
        runTests_->UnAdvise();
        return S_OK;
    }

    /// Not implemented
    STDMETHOD(OnAddInsUpdate)(SAFEARRAY * * custom)
    {
        return E_NOTIMPL;
    }

    /// Not implemented
    STDMETHOD(OnStartupComplete)(SAFEARRAY * * custom)
    {
        return E_NOTIMPL;
    }

    /// This method is called when the add-in is being disconnected.
    STDMETHOD(OnBeginShutdown)(SAFEARRAY * * custom)
    {
        // since storages contain references to cells we must free them
        db::storageMap<Model>().clear();
        db::storageMap<Option>().clear();
        db::storageMap<PricingMethod>().clear();
        db::storageMap<db::Result>().clear();
        return E_NOTIMPL;
    }

    /// A constructor.
    CPremiaAddin()
    {
        ::InitCommonControls();
    }



    DECLARE_PROTECT_FINAL_CONSTRUCT()

private:

    com_object<CWizardButtonHandler>    wizardMenu_;
    com_object<CRunTestsButton>         runTests_;


    /// Adds Premia wizard menu to Excel main menu bar
    void addPremiaMenu(Excel::_Application * spApp)
    {
        CComPtr<Office::_CommandBars>   spBars = spApp->CommandBars;
        CComPtr<Office::CommandBar>     spActiveBar = spBars->ActiveMenuBar;

        // position it below all toolbands
        //MsoBarPosition::msoBarTop = 1
        CComVariant vPos(1); 

        CComVariant vTemp(VARIANT_TRUE); // menu is temporary        
        CComVariant vEmpty(DISP_E_PARAMNOTFOUND, VT_ERROR);  

        CComQIPtr<Office::CommandBarControl> premiaMenu = 
            spActiveBar->GetControls()->Add(_variant_t(Office::msoControlPopup), vEmpty, vEmpty, vEmpty, vTemp);

        premiaMenu->Caption = _bstr_t("Premia");

        CComQIPtr<Office::CommandBarPopup> ppPremiaMenu = premiaMenu->GetControl();

        CComPtr < Office::CommandBarControls> spCmdBarCtrls = ppPremiaMenu->GetControls();

        CComVariant vMenuType(1); // type of control - menu
        CComVariant vMenuEmpty(DISP_E_PARAMNOTFOUND, VT_ERROR);
        CComVariant vMenuShow(VARIANT_TRUE); // menu should be visible
        CComVariant vMenuTemp(VARIANT_TRUE); // menu is temporary        


        Office::_CommandBarButtonPtr wizardMenu;

        // now create the actual menu item and add it
        wizardMenu = spCmdBarCtrls->Add(vMenuType, vMenuEmpty, vMenuEmpty, vMenuEmpty, vMenuTemp); 
        ATLASSERT(wizardMenu);

        wizardMenu->PutCaption(L"Premia Wizard...");
        wizardMenu->PutEnabled(VARIANT_TRUE);
        wizardMenu->PutVisible(VARIANT_TRUE); 

        wizardMenu_->Advise(wizardMenu);
    }
};

}