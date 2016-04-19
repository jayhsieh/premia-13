// WizardButtonHandler.h : déclaration de CWizardButtonHandler

#pragma once
#include "resource.h"       // symboles principaux



#if defined(_WIN32_WCE) && !defined(_CE_DCOM) && !defined(_CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA)
#error "Les objets COM monothread ne sont pas correctement pris en charge par les plates-formes Windows CE, notamment les plates-formes Windows Mobile qui ne prennent pas totalement en charge DCOM. Définissez _CE_ALLOW_SINGLE_THREADED_OBJECTS_IN_MTA pour forcer ATL à prendre en charge la création d'objets COM monothread et permettre l'utilisation de leurs implémentations. Le modèle de thread de votre fichier rgs a été défini sur 'Libre' car il s'agit du seul modèle de thread pris en charge par les plates-formes Windows CE non-DCOM."
#endif

namespace prxl {
    // CWizardButtonHandler

[
	coclass,
	threading(apartment),
	vi_progid("Excel_2.WizardButtonHandler"),
	progid("Excel_2.WizardButtonHandler.1"),
	version(1.0),
	uuid("F596A33E-C1D9-4BE0-9D88-C38A2BDA78E4"),
	helpstring("WizardButtonHandler Class")
]
class ATL_NO_VTABLE CWizardButtonHandler 
    : public IDispEventSimpleImpl<1,CWizardButtonHandler, &__uuidof(Office::_CommandBarButtonEvents)>
{
public:
	CWizardButtonHandler()
	{
	}

    template <class T>
    db::id<T> formatRegion(lib::Mode m, Excel::RangePtr &cell, Excel::RangePtr &zcb_cell)
    {
        return 
            wzd_.getPanel<T>().getCreationMode() == wzd::UseExisting
            ?   fmt::formatHeaderForExistingEntity<T>(wzd_.getPanel<T>().existing().getId(), cell)
            :   fmt::formatRegion(wzd_.get<T>(),m, zcb_cell,cell);
    }

    template <class T>
        void formatRegionEnumParameters(lib::Mode m, db::id<T> idx, Excel::RangePtr &cell)
    {
        if (wzd_.getPanel<T>().getCreationMode() == wzd::CreateNew)
            fmt::formatRegionEnumParameters(wzd_.get<T>(), m, idx, cell);
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

        if (pCommandBarButton->Caption == _bstr_t(L"Premia Wizard...") && g_app->ActiveWorkbook != 0)
        {
            // make the wizard modal and 
            // if a user presses Ok button then...
            if (wzd_.DoModal() == IDOK)
            {
                setUpDimensionsIfAny(wzd_);

                // switching the current directory to the Unix directory of Premia installation
                // TODO: Do this only when initialyield.dat file is touched
                // TODO: Switch back the current directory to the initial one
                SetCurrentDirectory(env::getDataPath());

                // start filling a range starting from the active cell in Excel
                // this variable points to the current cell for formatting
                // so after calling formatRegion it points to the next cell below a just created region
                Excel::RangePtr  cell = excelApp()->GetActiveCell();

                // retrieve Premia mode 
                lib::Mode m = wzd_.getMode();

                // if an option selected has an array parameter describing zcb interest rate structure
                // then a special region for that array is created and reference to it is written 
                // to this variable
                Excel::RangePtr zcb_cell;

                Model * mod = wzd_.get<Model>();
                Option * opt = wzd_.get<Option>();
                PricingMethod * met = wzd_.get<PricingMethod>();

                mod->Init(mod);
                opt->Init(opt,mod);
                met->Init(met, reinterpret_cast<Option*>(mod));

                // if a user selected to use existing model instance then just format the header for it
                // otherwise put also its parameters into the cell region
                db::id<Model> model_id = formatRegion<Model>(m, cell, zcb_cell);
                db::id<Option> option_id = formatRegion<Option>(m, cell, zcb_cell);
                db::id<PricingMethod> method_id = formatRegion<PricingMethod>(m, cell, zcb_cell);

                db::id<db::Result> result_id = fmt::formatResultRegion(wzd_.get<PricingMethod>(),m, cell, model_id, option_id, method_id);

                // create Compute button 
                // then create CButtonHandler and subscribe it to listening the button's events
                com_object<CButtonHandler>()->Advise(createComputeButton(Right(cell)), result_id);

                cell = cell->GetOffset(3, 0);

                formatRegionEnumParameters<Model>(m, model_id, cell);
                formatRegionEnumParameters<Option>(m, option_id, cell);
                formatRegionEnumParameters<PricingMethod>(m, method_id, cell);

                // resize columns to fit the content
                cell->Worksheet->Columns->AutoFit();

                // if there is an array parameter zcb create a zcb region
                // and put a reference to it to the value cell of the parameter
                //if (zcb_cell)
                //{
                //    fmt::zcb::createZcbRegion(zcb_cell, cell->GetOffset(4,0));
                //}

            }
        }
    }

    BEGIN_SINK_MAP(CWizardButtonHandler)
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

    /// Creates a Compute button in a given cell.
    /// \param cell - a place where to create the button
    /// \return a pointer to a created button
    MSForms::ICommandButtonPtr  createComputeButton(Excel::Range * cell)
    {
        Excel::OLEObjectsPtr oles =  cell->Worksheet->OLEObjects();

        Excel::_OLEObjectPtr obj = oles->Add("Forms.CommandButton.1", 
            vtMissing, vtMissing, vtMissing, vtMissing, vtMissing, vtMissing, 
            cell->Left, cell->Top, 64, 20);

        MSForms::ICommandButtonPtr button = obj->Object;

        button->Caption = "Compute";

        return button;
    }


private:
    /// Premia wizard dialog
    wzd::CWizard  wzd_;
    Office::_CommandBarButtonPtr btn_;
};

}