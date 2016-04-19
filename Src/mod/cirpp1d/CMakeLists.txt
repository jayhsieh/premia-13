add_library (libcirpp1d OBJECT
	cirpp1d_stdi/cf_cirpp1d_cap.c
	cirpp1d_stdi/cf_cirpp1d_floor.c
	cirpp1d_stdi/cf_cirpp1d_payerswaption.c
	cirpp1d_stdi/cf_cirpp1d_receiverswaption.c
	cirpp1d_stdi/cf_cirpp1d_zbcalleuro.c
	cirpp1d_stdi/cf_cirpp1d_zbputeuro.c
	cirpp1d_stdi/cf_cirpp1d_zcbond.c
	cirpp1d_stdi/cirpp1d_stdi.c
	cirpp1d_stdi/fd_cirpp1d_zbo.c
	cirpp1d_stdi/fd_cirpp1d_zcbond.c
	cirpp1d_stdi/tr_cirpp1d_zcbond.c
	cirpp1d_stdi/tr_cirpp1d_zbo.c
	cirpp1d_stdi/tr_cirpp1d_capfloor.c
	cirpp1d_stdi/tr_cirpp1d_swaption.c
	cirpp1d.c
	)
