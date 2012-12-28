#ifndef FLUX_BC_H
#define FLUX_BC_H "FluxBC.h"

//
// Flux BC class for CellWave
//
#include <iostream>
#include <string>
#include <vector>

#include "Overture.h"
#include "CompositeGridOperators.h"
#include "InterpolatePoints.h"
#include "ArraySimple.h"

#include "FluxBCBase.h"

class FluxBC : public FluxBCBase {
 public:
 
  FluxBC( );
  ~FluxBC( );

 private:
  FluxBC( FluxBC &X);
  FluxBC & operator=( FluxBC &X);
 public:

  void updateToMatchGrid( CompositeGrid &cg, 
			  int bcID= standardFluxBCID);
  void setCompositeGrid( CompositeGrid &cg ); 

  void setupInterpolation(const int numberOfComponentFields=1);

  void applyBoundaryCondition( realCompositeGridFunction &qto,
			       int ic);

  //..data
  //InterpGridArray interpGridEdges;
  ArraySimple<int>        edgeStart,  edgeEnd;
  ArraySimple<bool>       isInterpolatedEdge;
  ArraySimple<realArray>  edgeCoords;

  realArray    xyzInterpolate;
  realArray    valuesInterpolate;
  IntegerArray ignoreGrid;
  
};

#endif

