#ifndef INTERPOLATE_FLUX_BOUNDARY_H 
#define INTERPOLATE_FLUX_BOUNDARY_H

// define some interpolation functions

class OGFunction;
// extern intArray Overture::nullIntegerDistributedArray();
  
int
interpolatePointsPF(int iThisGrid,
		    const realArray & positionToInterpolate,
		  const realCompositeGridFunction & u,
		  realArray & uInterpolated, 
		  intArray & indexGuess=Overture::nullIntegerDistributedArray(),
                  intArray & interpoleeGrid=Overture::nullIntegerDistributedArray(),
                  intArray & wasInterpolated=Overture::nullIntegerDistributedArray());



#endif
