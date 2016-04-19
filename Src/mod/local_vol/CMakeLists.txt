add_library (liblocal_vol OBJECT
	local_vol.c
	local_vol_callable/highly_path_dependent_protection.c
	local_vol_callable/no_call_protection.c
	local_vol_callable/intermittent_protection.c
	local_vol_callable/path_dependent_protection.c
	local_vol_callable/local_vol_callable.c
	local_vol_callable/mc_ribsde.c
	local_vol_callable/standard_protection.c
	)
