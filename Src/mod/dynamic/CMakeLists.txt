add_library (libdynamic OBJECT
	dynamic.c
	dynamic_stdc/dynamic_stdc.c
	dynamic_stdc/herbertsson_cds.c
	dynamic_stdndc/dynamic_stdndc.c
	dynamic_stdndc/RogersDiGraziano.c
	dynamic_stdndc/LCF_cdo_hedging.c
	dynamic_stdndc/FB_cdo_hedging.c
	dynamic_stdndc/EFVH_pricing_cdo.c
	dynamic_stdndc/Herbertsson_cdo.c
	)
