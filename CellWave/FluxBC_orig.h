#ifndef FLUX_BC_H
#define FLUX_BC_H "FluxBC.h"

//
// Flux BC class for CellWave
//
#include <iostream>
#include <string>

#include "Overture.h"
#include "CompositeGridOperators.h"
#include "InterpolatePoints.h"

//#include "interpolateFluxBoundary.h"

#define FLUXBC_ForBoundary(side,axis) \
    for( axis=0; axis<mg.numberOfDimensions(); axis++ ) \
      for( side=0; side<=1; side++ )   

class OGFunction;
// extern intArray Overture::nullIntegerDistributedArray();

class FluxBC {
 public:

  enum FluxBCType {
    noFluxBC          = 0,
    interpolateFluxBC = 1
  }; 

  enum StandardBCID {
    standardFluxBCID=2
  };

  FluxBC( );
  ~FluxBC( );

 private:
  FluxBC( FluxBC &X);
  FluxBC & operator=( FluxBC &X);
 public:

  virtual void updateToMatchGrid( CompositeGrid &cg, 
			  int bcID= standardFluxBCID);
  virtual void setCompositeGrid( CompositeGrid &cg ); 
  
  virtual void  setInterpolant( Interpolant &interpolant );
  virtual void  setOperators( CompositeGridOperators &operators );
  
  virtual FluxBCType setFluxBCType( std::string &name );  
  virtual FluxBCType convertToFluxBCType( std::string &name );
  virtual FluxBCType getFluxBCType();

  virtual void setFluxBoundaryID( int id );
  virtual int  getFluxBoundaryID();

  virtual void   setFluxCoefficient( double f );
  virtual double getFluxCoefficient( );

  virtual void    setDebug( int debug_ );
  virtual int     getDebug();
  virtual double  getInterpTimer();

  virtual void applyBoundaryCondition( realCompositeGridFunction &qfrom,
				       realCompositeGridFunction &qto,
				       int ic);


  // define some interpolation functions
  
  int interpolatePoints(int iThisGrid,
			const realArray & positionToInterpolate,
			const realCompositeGridFunction & u,
			realArray & uInterpolated, 
			intArray & indexGuess=Overture::nullIntegerDistributedArray(),
			intArray & interpoleeGrid=Overture::nullIntegerDistributedArray(),
			intArray & wasInterpolated=Overture::nullIntegerDistributedArray());
  
  //..data
  InterpolatePoints          interpolator;

  CompositeGrid              *pCg;
  Interpolant                *pInterp;
  CompositeGridOperators     *pOp;

  int                        debug;         // general debug flag
  int                        interpDebug;   // debug flag for interpolation, bitflags 1+2+4+8;

  int                        idFluxBoundary; 
  FluxBCType                 fluxBCType;
  double                     fluxCoefficient; // F=fluxCoeff;  u_n = F[ u ], F should be a realArray

  double                     timerInterpCode;
};

#endif

