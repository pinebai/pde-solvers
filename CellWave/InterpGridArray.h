#ifndef INTERP_GRID_ARRAY_H
#define INTERP_GRID_ARRAY_H "InterpGridArray.h"

#include "Overture.h"
#include "InterpolatePoints.h"
#include <vector>
#include <stdio.h>
#include <iostream>

//#include "InterpolatePoints.h"

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

  //..data
  int  igrid, iside, iaxis;

};

class InterpGridArray {
 public:
  InterpGridArray() : pInterp( NULL ) {};
  ~InterpGridArray() { this->freeMemory();  };

 private:
  InterpGridArray( const InterpGridArray &X);
  InterpGridArray &operator=( const InterpGridArray &X);

 public:
  void freeMemory();
  bool isAllocated();
  void resize( const GridAxisSide &graxsi );
  const GridAxisSide &getDimensions() { return dimensions; };

  InterpolatePoints & get( const GridAxisSide &graxsiIndex );
  void set( const InterpolatePoints &interp, GridAxisSide &graxsiIndex );

  InterpolatePoints *pInterp;
  GridAxisSide dimensions;
};

#endif
