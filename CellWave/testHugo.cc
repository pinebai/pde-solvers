#include <vector>
#include <stdio.h>
#include <iostream>

class GridAxisSide;
class GridAxisSide {
public:
  GridAxisSide() : igrid( -1 ), iaxis(-1), iside( -1 ) {};
  ~GridAxisSide() {};

  void setGrid( int j ) { igrid=j; };
  void setAxis( int j ) { iaxis=j; };
  void setSide( int j ) { iside=j; };

  int  getGrid() const { return igrid; };
  int  getAxis() const { return iaxis; };
  int  getSide() const { return iside; };

  int getSize()  const 
  {
    int  igrid, iside, iaxis;
    igrid = this->getGrid();
    iaxis = this->getAxis();
    iside = this->getSide();

    int size= igrid*iside*iaxis;

    return( size );
  };

  int getOffset( GridAxisSide &dimensions ) const  {
    int  igrid, iside, iaxis;
    igrid = this->getGrid();
    iaxis = this->getAxis();
    iside = this->getSide();

    int  ngrids, nsides, naxes;
    ngrids = dimensions.getGrid();
    naxes = dimensions.getAxis();
    nsides = dimensions.getSide();

    int offset= iside+(nsides*iaxis) + naxes*nsides*igrid;

    return offset;
  };

  int  igrid, iside, iaxis;

};

class Hugo {
public:
  Hugo()    : myNameIs("Hugo")
     { };
  ~Hugo() {};

  void setGridAxisSide( const GridAxisSide &graxsi )
  { gridAxisSide = graxsi; };

  void print() 
  { 
    int  igrid, iside, iaxis;
    igrid = gridAxisSide.getGrid();
    iaxis = gridAxisSide.getAxis();
    iside = gridAxisSide.getSide();
    printf( "My name is %s (grid %d, axis %d, side %d).\n", myNameIs.c_str(), 
	    igrid, iaxis, iside);
  }

  void zero() { 
    gridAxisSide.setGrid(0);  gridAxisSide.setAxis(0); gridAxisSide.setSide(0);
  }

  //..data
  GridAxisSide gridAxisSide;
  std::string myNameIs;
};

class HugoArray {
public:
  HugoArray() : pHugo( NULL ) {};
  ~HugoArray() { this->freeMemory();  };

  void freeMemory();
  bool isAllocated();
  void resize( const GridAxisSide &graxsi );
  Hugo & get( const GridAxisSide &graxsiIndex );
  void set( const Hugo &hugo, GridAxisSide &graxsiIndex );

  Hugo *pHugo;
  GridAxisSide dimensions;
};

void HugoArray::freeMemory() {
  if( this->isAllocated() ) {
    delete [] pHugo;
  }
}

bool HugoArray::isAllocated() { return pHugo!=NULL; };

void HugoArray::resize( const GridAxisSide &graxsi )
{
  this->freeMemory();
  dimensions = graxsi;
  int size=dimensions.getSize();
  pHugo = new Hugo[ size ];
}

Hugo& HugoArray::get( const GridAxisSide &graxsi)
{
  int offset=graxsi.getOffset( this->dimensions );
  printf( " offset = %4d; ",offset);
  return pHugo[offset];
}

//
// ..we want a grid/edge array
//   .. for each grid, each axis, each side, we store an Interpolator
//   .. Hence, interpEdge[NGRIDS][NAXES][NSIDES]; //!!
//

int main(int argc, char **argv)
{
  const int nGrids=5;
  const int nAxis=3;
  const int nSides=2;
  
  //Hugo gridEdgeHugos[nGrids][nAxis][nSides];
  GridAxisSide dimensions;
  dimensions.setGrid( nGrids );
  dimensions.setAxis( nAxis );
  dimensions.setSide( nSides );

  HugoArray gridEdgeHugo;
  gridEdgeHugo.resize( dimensions );

  for( int ig=0; ig< nGrids; ++ig ) {
    for( int iaxis=0; iaxis< nAxis; iaxis++ ) { 
      for( int iside=0; iside< nSides; ++iside ) {
	GridAxisSide graxsi;
	graxsi.setGrid( ig );
	graxsi.setAxis( iaxis );
	graxsi.setSide( iside );
	//gridEdgeHugos[ig][iaxis][iside].setGridAxisSide( graxsi );
        gridEdgeHugo.get( graxsi ).setGridAxisSide( graxsi );
	printf("  set grid %d, axis %d, side %d\n",
	       ig, iaxis, iside );
      }
    }
  } // end for

  //..now print them
  for( int ig=0; ig< nGrids; ++ig ) {
    for( int iaxis=0; iaxis< nAxis; iaxis++ ) { 
      for( int iside=0; iside< nSides; ++iside ) {
	//gridEdgeHugos[ig][iaxis][iside].print();
	GridAxisSide graxsi;
	graxsi.setGrid( ig );
	graxsi.setAxis( iaxis );
	graxsi.setSide( iside );
	gridEdgeHugo.get( graxsi ).print();
      }
    }
  }
  printf("..done..\n");
}
