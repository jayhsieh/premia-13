add_library (libsg1d OBJECT
	sg1d.c
	sg1d_stdi/sg1d_stdi.c
	sg1d_stdi/QuadraticModel.c
	sg1d_stdi/cf_quadratic1d_zcbond.c
	sg1d_stdi/cf_quadratic1d_zbputeuro.c
	sg1d_stdi/cf_quadratic1d_zbcalleuro.c
	sg1d_stdi/cf_quadratic1d_cap.c
	sg1d_stdi/cf_quadratic1d_floor.c
	sg1d_stdi/cf_quadratic1d_payerswaption.c
	sg1d_stdi/cf_quadratic1d_receiverswaption.c
	sg1d_stdi/tr_quadratic1d_zcbond.c
	sg1d_stdi/tr_quadratic1d_zbo.c
	sg1d_stdi/tr_quadratic1d_capfloor.c
	sg1d_stdi/tr_quadratic1d_swaption.c
	sg1d_stdi/tr_quadratic1d_bermudianswaption.c
	)
