#include "Overture.h"
#include <stdio.h>
#include <string>
#include "InterpolatePoints.h"

int main(int argc, char **argv)
{
  Overture::start( argc, argv );

  CompositeGrid cg;
  std::string gridName="Grids/cell20um.hdf";
  int npoints=4;
  int ignoreGrid=-1;

  if( argc>1 ) gridName = argv[1];
  else { 
    printf("usage: %s <grid name.hdf> <#points>\nUsing default grid %s\n",
	   argv[0], gridName.c_str());
  }
  if( argc>2 ) npoints=atoi(argv[2]);
  //if( argc>3 ) ignoreGrid=atoi(argv[3]);

  getFromADataBase( cg, gridName.c_str() );
  cg.update();
  int ndim   = cg.numberOfDimensions();

  printf("..Interpolating from %d dimensional grid '%s'\n",
	 ndim, gridName.c_str());
  printf("  at %d coordinate points. ", npoints);
  //  if( ignoreGrid>0 ) printf("Ignoring grid %d.", ignoreGrid);
  //printf("\n");
	 
  InterpolatePoints interp;
  realArray posToInterp( npoints, ndim );
  IntegerArray ignoreGrids(npoints);
  for(int i=0; i<npoints; ++i) {
    real dx=(2./npoints)*i;
    posToInterp(i, 0) = 2.+dx;
    posToInterp(i, 1) = 2.+dx;
    posToInterp(i, 2) = 0;
    ignoreGrids(i) = (i<5);
  }


  printf("..Passing in xyz array of dimensions [%d, %d]\n",
	 posToInterp.getLength(0), posToInterp.getLength(1));

  printf("..before 'buildInterpolationInfo\n");fflush(0);
  interp.buildInterpolationInfo(posToInterp, cg, ignoreGrids );
  printf("..after  'buildInterpolationInfo\n");fflush(0);

  //..output the interp info
  IntegerArray indexValues, interpoleeGrid;
  interp.getInterpolationInfo(cg, indexValues, interpoleeGrid);
  indexValues.display("IndexValues");
  interpoleeGrid.display("InterpoleeGrid");

  Overture::finish();
}
