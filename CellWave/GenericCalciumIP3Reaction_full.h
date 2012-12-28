/// Brief description: CellWave::GenericReaction is the rxn base class 
///
///   CellWave::GenericCalciumIP3Reaction  -- base class for Ca2+/IP3 dynamics
///     part of CellWave
///
#ifndef CELLWAVE_GENERIC_CALCIUM_IP3_REACTION_H
#define CELLWAVE_GENERIC_CALCIUM_IP3_REACTION_H

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "GenericReaction.h"

namespace CellWave {

  class GenericCalciumIP3Reaction : public GenericReaction {
  public:

    GenericCalciumIP3Reaction() 
      { // do nothing
      };

    ~GenericCalciumIP3Reaction()
      { // do nothing
      };

  private: // do not use copy constructors for now
    GenericCalciumIP3Reaction( const GenericCalciumIP3Reaction &X) 
      { };
    GenericCalciumIP3Reaction & operator=( const GenericCalciumIP3Reaction &X)
      { };
  public:

    /* virtual */ double  getDiffusionCoefficient( const int component);
    /* virtual */ bool    isDiffusive( const int component);

    /* virtual */ std::string getTitle();
    /* virtual */ std::string getLongComponentName( const int component );
    /* virtual */ std::string getShortComponentName( const int component );

    /* virtual */ int     getNumberOfSpecies();
    /* virtual */ double  getLengthScale();

    //.. set parameters
    /* virtual */ bool
      readParameterFile( CellWave::ParameterReader &param);

    //.. looping calls to array interface, for initial data & rhs
    /* virtual */ void 
      callInitialDataLoop( const double &tcomp,
			   const int&nd,    const int &ncomp,
			   const int &nd1a, const int &nd1b,
			   const int &nd2a, const int &nd2b,
			   const int &nd3a, const int &nd3b,
			   const int &n1a,  const int &n1b,
			   const int &n2a,  const int &n2b,
			   const int &n3a,  const int &n3b,
			   const double &x, double &q );

    /* virtual */ void 
      callRHSLoop( const double &tcomp,
		   const int&nd,    const int &ncomp,
		   const int &nd1a, const int &nd1b,
		   const int &nd2a, const int &nd2b,
		   const int &nd3a, const int &nd3b,
		   const int &n1a,  const int &n1b,
		   const int &n2a,  const int &n2b,
		   const int &n3a,  const int &n3b,
		   const double &x, const double &q, 
		   double &rhs);
  
  };//end class GenericCalciumIP3Reaction

}; //end namespace CellWave
#endif
