#include <iostream>
#include <stdio.h>
#include "Overture.h"

#include "Display.h"

int main(int argc, char **argv)
{
  Overture::start( argc, argv);
  aString nameOfOGFile="Grids/test2d-cell4.hdf";
  if( argc >1 ) {
    nameOfOGFile = argv[1];
  }
  std::cout << argv[0] << " using grid=" 
	    << nameOfOGFile << "..." << std::endl;

  Display disp; //David's display

  CompositeGrid cg;
  getFromADataBase( cg, nameOfOGFile );
  cg.update( MappedGrid::THEvertex | MappedGrid::THEcenter 
            | MappedGrid::THEvertexBoundaryNormal);

  for( int ig=0; ig< cg.numberOfComponentGrids() ; ++ig ) {
    MappedGrid &mg = cg[ig];
    for( int side=0; side<2; ++side ) {
      for( int axis=0; axis<cg.numberOfDimensions(); ++axis) {
	realArray & normal = mg.vertexBoundaryNormal(side,axis);
	printf("--grid %d, side %d, axis %d, bc %d--\n", 
	       ig, side, axis, mg.boundaryCondition()(side,axis));
	disp.display( normal );
      }
    }
  }

  Overture::finish();
}
