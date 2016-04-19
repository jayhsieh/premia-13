add_library (libcir1d OBJECT
	cir1d_stdi/cf_cir1d_zbcalleuro.c
	cir1d_stdi/cf_cir1d_zbputeuro.c
	cir1d_stdi/cf_cir1d_zcbond.c
	cir1d_stdi/cir1d_stdi.c
	cir1d_stdi/fd_explicit_cir1d_capfloor.c
	cir1d_stdi/fd_explicit_cir1d_swaption.c
	cir1d_stdi/fd_explicit_cir1d_zbo.c
	cir1d_stdi/fd_explicit_cir1d_zcbond.c
	cir1d_stdi/fd_gauss_cir1d_capfloor.c
	cir1d_stdi/fd_gauss_cir1d_swaption.c
	cir1d_stdi/fd_gauss_cir1d_zbo.c
	cir1d_stdi/fd_gauss_cir1d_zcbond.c
	cir1d_stdi/mc_teichmannbayer.c
	cir1d.c
	)
