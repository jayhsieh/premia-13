// ButtonHandler.h : déclaration de CButtonHandler

#pragma once
#include "resource.h"       // symboles principaux

#include <direct.h>
#include "run_method.h"

#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Les objets COM monothread ne sont pas correctement pris en charge par les plates-formes Windows CE, notamment les plates-formes Windows Mobile qui ne prennent pas totalement en charge DCOM. Définissez _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA pour forcer ATL à prendre en charge la création d'objets COM monothread et permettre l'utilisation de leurs implémentations. Le modèle de thread de votre fichier rgs a été défini sur 'Libre' car il s'agit du seul modèle de thread pris en charge par les plates-formes Windows CE non-DCOM."
#endif

namespace prxl {
extern _ATL_FUNC_INFO OnClickButtonInfo_2;


/// This class handles clicks on Compute button.
/// It holds ID of a button it is connected to
[
	coclass,
	threading(apartment),
	vi_progid("Excel_2.ButtonHandler"),
	progid("Excel_2.ButtonHandler.1"),
	version(1.0),
	uuid("A37FEBB2-19EC-43D6-B828-38555D67BF08"),
	helpstring("ButtonHandler Class")
]
class ATL_NO_VTABLE CButtonHandler :
    public IDispEventSimpleImpl<1,CButtonHandler, &__uuidof(MSForms::CommandButtonEvents)>
{
    /// id of a button the handler is connected to
    db::id<db::Result> id_;
public:

    CButtonHandler() : id_(db::id<db::Result>::notExist()) {}

    /// Handles events from the button
    STDMETHOD(Invoke)(DISPID dispidMember, REFIID /*riid*/,
        LCID lcid, WORD /*wFlags*/, DISPPARAMS* pdispparams, VARIANT* pvarResult,
        EXCEPINFO* /*pexcepinfo*/, UINT* /*puArgErr*/)
    {
        if (dispidMember == -600)
        {
            // switching the current directory to the Unix directory of Premia installation
            // TODO: Do this only when initialyield.dat file is touched
            // TODO: Switch back the current directory to the initial one
            SetCurrentDirectory(env::getDataPath());

            runMethod(id_);

        }
        return S_OK;
    }

    BEGIN_SINK_MAP(CButtonHandler)
//        SINK_ENTRY_INFO(1, __uuidof(MSForms::CommandButtonEvents),/*dispid*/ -600, OnClickButton, &OnClickButtonInfo_2)
    END_SINK_MAP()

#undef LONG

    /// Starts listening to the button events
    /// \param button - a button which events are to be listened by the handler
    /// \param id - id of the button
    void Advise(MSForms::ICommandButtonPtr button, db::id<db::Result> id)
    {
        id_ = id;
        DispEventAdvise(button,&__uuidof(MSForms::CommandButtonEvents));
    }


	DECLARE_PROTECT_FINAL_CONSTRUCT()

};

}