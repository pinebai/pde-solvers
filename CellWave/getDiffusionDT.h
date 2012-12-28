real
getDiffusionDT(const real & cfl, 
	       int numberOfComponents,
	       realArray &nuArray,
	       CompositeGrid & cg,
	       const real alpha0 = -2.,
	       const real beta0  = 1. );

void
getDiffusionDT(const real & cfl, 
	       int numberOfComponents,
	        realArray &nuArray,
	       realArray &dtArray,
	       CompositeGrid & cg,
	       const real alpha0 = -2.,
	       const real beta0  = 1. );

void
getDiffusionDT(const real & cfl, 
	       int numberOfComponents,
	         realArray &nuArray,
	       realArray &dtArray,
	       realArray &dtGridArray,
	       CompositeGrid & cg,
	       const real alpha0 = -2.,
	       const real beta0  = 1. );

