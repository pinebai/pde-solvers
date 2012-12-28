#ifndef INTERPOLATE_POINTS_CG_H
#define INTERPOLATE_POINTS_CG_H


// class that can be used to interpolate arbitrary points on a Composite grid.


class InterpolatePoints
{
 public:

  enum InterpolationStatusEnum
  {
    notInterpolated=0,
    interpolated,
    extrapolated
  };

  InterpolatePoints();
  ~InterpolatePoints();
  
  int interpolatePoints(const realArray & positionToInterpolate,
		  const realCompositeGridFunction & u,
		  realArray & uInterpolated, 
		  const Range & R0=nullRange,           
		  const Range & R1=nullRange,
		  const Range & R2=nullRange,
		  const Range & R3=nullRange,
		  const Range & R4=nullRange );

  int interpolationCoefficients(const realArray & positionToInterpolate,
		  const CompositeGrid &cg,
		  realArray & uInterpolationCoeff);

  int interpolateAllPoints(const realCompositeGridFunction & uFrom,
			   realCompositeGridFunction & uTo );
  
  int interpolateAllPoints(const realCompositeGridFunction & uFrom,
			   realMappedGridFunction & uTo );
  

  // Here is a way to interpolate in two steps, this will save the interpolation info
  int buildInterpolationInfo(const realArray & positionToInterpolate,CompositeGrid & cg, 
			     int ignoreGrid);
  int buildInterpolationInfo(const realArray & positionToInterpolate,CompositeGrid & cg, 
			     IntegerArray &ignoreGrid = Overture::nullIntArray() );
  int interpolatePoints(const realCompositeGridFunction & u,
			realArray & uInterpolated, 
			const Range & R0=nullRange,           
			const Range & R1=nullRange,
			const Range & R2=nullRange,
			const Range & R3=nullRange,
			const Range & R4=nullRange );

  // return the status array for the last interpolation.
  const IntegerArray & getStatus() const;

  // return the index values and interpoleeGrid for the last interpolation.
  int getInterpolationInfo(CompositeGrid & cg, IntegerArray & indexValues, IntegerArray & interpoleeGrid) const;

 protected:

  IntegerArray *indirection;
  IntegerArray *interpolationLocation;
  IntegerArray *interpolationLocationPlus;
  RealArray *interpolationCoordinates;

  IntegerArray numberOfInterpolationPoints;
  IntegerArray status;

};

#endif
