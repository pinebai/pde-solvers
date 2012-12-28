#include "Overture.h"

#include <stdio.h>
#include <iostream>

#include "InterpGridArray.h"
  
int main(int argc, char **argv)
{
  Overture::start( argc, argv );
  const int nGrids=5;
  const int nAxis=3;
  const int nSides=2;

  printf("........testInterpGridArray -- test classes under fluxbc.......\n");

  CompositeGrid cg;
  aString gridName="Grids/cell20um.hdf";
  if (argc>1) {
    gridName = argv[1];
  } else {
    printf("usage: %s <grid.hdf>\n", argv[0]);
  }

  getFromADataBase( cg, gridName );

  GridAxisSide dimensions;
  dimensions.setGrid( nGrids );
  dimensions.setAxis( nAxis );
  dimensions.setSide( nSides );

  InterpGridArray interpGridEdges;
  interpGridEdges.resize( dimensions );

  for( int ig=0; ig< nGrids; ++ig ) {
    for( int iaxis=0; iaxis< nAxis; iaxis++ ) { 
      for( int iside=0; iside< nSides; ++iside ) {
        GridAxisSide graxsi;
	const int npoints=3;
	realArray posToInterp(npoints, cg.numberOfDimensions());
	posToInterp = 0.;
        graxsi.setGrid( ig );
        graxsi.setAxis( iaxis );
        graxsi.setSide( iside );
        //interpGridEdges.get( graxsi ).buildInterpolationInfo(posToInterp, cg );
	InterpolatePoints &interp=interpGridEdges.get( graxsi );
	interp.buildInterpolationInfo(posToInterp, cg );

        printf("  setup grid %d, axis %d, side %d\n",
               ig, iaxis, iside );
      }
    }
  } // end for

  //..now print them
  for( int ig=0; ig< nGrids; ++ig ) {
    for( int iaxis=0; iaxis< nAxis; iaxis++ ) { 
      for( int iside=0; iside< nSides; ++iside ) {
        GridAxisSide graxsi;
        graxsi.setGrid( ig );
        graxsi.setAxis( iaxis );
        graxsi.setSide( iside );
        InterpolatePoints &interp = interpGridEdges.get( graxsi );
	//.. and HERE, use the interp to set flux bcs...
      }
    }
  }
  printf("..done..\n");

  Overture::finish();
}
