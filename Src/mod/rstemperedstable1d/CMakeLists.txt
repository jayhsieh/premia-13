add_library (librstemperedstable1d OBJECT
	rstemperedstable1d.c
	rstemperedstable1d_lim/rstemperedstable1d_lim.c
	rstemperedstable1d_lim/ap_fastwhbar_rsts.c
	rstemperedstable1d_std/rstemperedstable1d_std.c
	rstemperedstable1d_std/ap_fastwhdig_rsts.c
	rstemperedstable1d_std/ap_fastwh_rsts.c
	)
