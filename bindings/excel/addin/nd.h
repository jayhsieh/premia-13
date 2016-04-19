#pragma once

namespace prxl
{
    /// Dialog box allowing to enter number of dimensions of a model
    class NumberOfDimensionsDialog : public CSimpleDialog<IDD_DIALOG2>
    {
        /// the edit-box
        WTL::CEdit  m_ndim;
        /// number of dimensions
        int value_;

    public:

        /// A constructor.
        /// \param initial_value - initial value for the number of dimensions
        NumberOfDimensionsDialog(int initial_value) : value_(initial_value) {}

        enum { IDD = IDD_DIALOG2 };

        BEGIN_MSG_MAP(NumberOfDimensionsDialog)
            MESSAGE_HANDLER(WM_INITDIALOG, OnInitDialog)
            COMMAND_HANDLER(IDC_NDIM, EN_CHANGE, OnChange)
            CHAIN_MSG_MAP(CSimpleDialog<IDD_DIALOG2>)
        END_MSG_MAP()

        /// Initializes dialog controls.
        /// Called on dialog creation.
        LRESULT OnInitDialog(UINT, WPARAM, LPARAM, BOOL&)
        {
            m_ndim = GetDlgItem(IDC_NDIM);

            m_ndim.SetFocus();

            m_ndim.AppendText(toString(value_));

            return 0;
        }

        /// Updates value_. Called whenever the edit-box string is changed.
        LRESULT OnChange(WORD d1, WORD d2, HWND d3, BOOL& d4)
        {
            TCHAR buffer[10]; 
            buffer[0] = 9;
            m_ndim.GetLine(0, buffer);
            value_ = _wtoi(buffer);
            return 0;
        }

        /// \return number of dimensions
        int get() 
        {
            return value_;
        }
    };

    /// Detects whether the model is n-dimensional
    /// \param model - a pointer to the model
    /// \return number of dimensions if the model is n-dimensional and -1 otherwise
    inline int getDimensions(Model * model)
    {
        VAR * v = lib::findParameter(model, "Model Size");

        return v ? v->Val.V_INT : -1;
    }

    /// Resizes all double arrays in 'vars' to have the same length
    /// \param vars - array of VARs to be inspected
    /// \param stopper - maximum number of VARs to be processed
    /// \param nDim - the length that all double arrays will have
    inline void setUpDimensions(VAR * vars, int stopper, int nDim)
    {
        for (; stopper && vars->Vtype != PREMIA_NULLTYPE; ++vars, --stopper)
        {
            switch (vars->Vtype)
            {
            case PNLVECT:
                {
                    pnl_vect_resize_from_double(vars->Val.V_PNLVECT, nDim, 0);
                    break;
                }
            case PNLVECTCOMPACT:
                {
					double e = pnl_vect_compact_get(vars->Val.V_PNLVECTCOMPACT, 0);
                    pnl_vect_compact_resize(vars->Val.V_PNLVECTCOMPACT, nDim, e);
                    break;
                }
            case NUMFUNC_1:
                {                
                    setUpDimensions(vars->Val.V_NUMFUNC_1->Par, MAX_PAR, nDim);
                    break;
                }
            case NUMFUNC_2:
                {
                    setUpDimensions(vars->Val.V_NUMFUNC_2->Par, MAX_PAR, nDim);
                    break;
                }
            }
        }
    }

    /// Changes model dimension
    /// \param model - model dimension of which is to be changed
    /// \param nDim - number of dimensions
    inline void setUpDimensions(Model * model, int nDim)
    {
        if (VAR * v = lib::findParameter(model, "Model Size"))
        {
            v->Val.V_INT = nDim;
            setUpDimensions(reinterpret_cast<VAR*>(model->TypeModel), model->nvar, nDim);
        }
    }

    /// Changes option dimension
    /// \param option - option dimension of which is to be changed
    /// \param nDim - number of dimensions
    inline void setUpDimensions(Option * option, int nDim)
    {
        setUpDimensions(reinterpret_cast<VAR*>(option->TypeOpt), option->nvar_setable, nDim);
    }

    /// Changes pricing method dimension
    /// \param method - method dimension of which is to be changed
    /// \param nDim - number of dimensions
    inline void setUpDimensions(PricingMethod * method, int nDim)
    {
        setUpDimensions(method->Par, MAX_PAR, nDim);
    }

    /// If a model in the wizard is n-dimensional then it set appropriate number of dimensions for 
    /// the model, option and pricing method at hand.
    /// \param w - wizard dialog
    inline void setUpDimensionsIfAny(wzd::CWizard &w)
    {
        int ndim; 

        if (w.getPanel<Model>().getCreationMode() == wzd::CreateNew && (ndim = getDimensions(w.get<Model>()) != -1))
        {
            NumberOfDimensionsDialog  ndim_dialog(ndim);

            ndim_dialog.DoModal();

            ndim = ndim_dialog.get();

            setUpDimensions(w.get<Model>(), ndim);

            if (w.getPanel<Option>().getCreationMode() == wzd::CreateNew)
                setUpDimensions(w.get<Option>(), ndim);

            // TODO: add handling for incompatibilities between existing model-option-method dimensions

            if (w.getPanel<PricingMethod>().getCreationMode() == wzd::CreateNew)
                setUpDimensions(w.get<PricingMethod>(), ndim);
        }
    }
}