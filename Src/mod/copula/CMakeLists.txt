add_library (libcopula OBJECT
	copula.c
	copula_stdndc/copula_stdndc.c
	copula_stdndc/stein.c
	copula_stdndc/hull_white.c
	copula_stdndc/laurent_gregory.c
	copula_stdndc/price_cdo.c
	copula_stdndc/monte_carlo.c
	copula_stdndc/saddlepoint.c
	)
