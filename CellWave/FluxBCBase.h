#ifndef FLUX_BC_BASE_H
#define FLUX_BC_BASE_H "FluxBCBase.h"

//
// Flux BC Base class for CellWave
//
#include <iostream>
#include <string>

#include "Overture.h"
#include "CompositeGridOperators.h"
#include "InterpolatePoints.h"

//#include "interpolateFluxBoundary.h"

class OGFunction;
// extern intArray Overture::nullIntegerDistributedArray();

class FluxBCBase {
 public:

  enum FluxBCType {
    noFluxBC          = 0,
    interpolateFluxBC = 1
  }; 

  enum StandardBCID {
    standardFluxBCID=2
  };

  FluxBCBase( );
  ~FluxBCBase( );

 private:
  FluxBCBase( FluxBCBase &X);
  FluxBCBase & operator=( FluxBCBase &X);
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

  virtual double  getInterpSetupTimer();
  virtual double  getInterpTimer();

  virtual void applyBoundaryCondition( realCompositeGridFunction &qfrom,
				       realCompositeGridFunction &qto,
				       int ic, double tcomp=0.);

  virtual void applyBoundaryCondition( realCompositeGridFunction &qto,
				       int ic) {};

  virtual void useTwilightZoneFlow( OGFunction &tzSolution );
  virtual void doNotUseTwilightZoneFlow();
  virtual bool isUsingTwilightZoneFlow();

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

  double                     timerInterpSetupCode;
  double                     timerInterpCode;

  //..Twilight zone flow
  OGFunction *pExactSolution;
  bool        useTZ;
};

#endif



