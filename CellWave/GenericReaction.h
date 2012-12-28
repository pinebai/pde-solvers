/// Brief description: CellWave::GenericReaction is the rxn base class 
///
///   CellWave::GenericReaction  -- abstract base class
///     part of CellWave
///
#ifndef CELLWAVE_GENERICREACTION_H
#define CELLWAVE_GENERICREACTION_H

#include <math.h>
#include <stdio.h>
#include <iostream>
#include <string>

#include "ParameterReader.h"

namespace CellWave {

  class GenericReaction {
  private:
    std::string          reactionName;
  protected:
    void                 setReactionName( const std::string &name);
  public:

    GenericReaction() 
      { // do nothing
      };

    ~GenericReaction()
      { // do nothing
      };

  private: // do not use copy constructors for now
    GenericReaction( const GenericReaction &X) 
      { };
    GenericReaction & operator=( const GenericReaction &X)
      { };
  public:
    std::string     getReactionName();

    virtual double  getDiffusionCoefficient( const int component);
    virtual bool    isDiffusive( const int component);

    virtual double  getFluxBCCoefficient( const int component );
    virtual bool    hasFluxBC( const int component );

    virtual std::string getTitle();
    virtual std::string getLongComponentName( const int component );
    virtual std::string getShortComponentName( const int component );

    virtual int     getNumberOfSpecies();
    virtual double  getMaximumLengthScale() { return -1.; };
    virtual void    setMaximumLengthScale( double len ) { };

    double  getMinimumDiffusionLimit() { return 1e-16;} //visc<this is ==0.

    //.. set parameters
    virtual bool
      readParameterFile( CellWave::ParameterReader &param);

    virtual int
      getNumberOfDimensions() { return 0; };

    //.. print out parameters
    virtual void
      printParameters( int iPrintChannel );

    void printOneParameter(int ip,std::string name,double param,double def);
    void printOneParameter(int ip,std::string name,double param,int def);
    void printOneParameter(int ip,std::string name,int param, int def);
    void printOneParameter(int ip,std::string name,
			   std::string param,std::string def);

    //.. looping calls to array interface, for initial data & rhs
    virtual void 
      callInitialDataLoop( const double &tcomp,
			   const int&nd,    const int &ncomp,
			   const int &nd1a, const int &nd1b,
			   const int &nd2a, const int &nd2b,
			   const int &nd3a, const int &nd3b,
			   const int &n1a,  const int &n1b,
			   const int &n2a,  const int &n2b,
			   const int &n3a,  const int &n3b,
			   const double &x, double &q );

    virtual void 
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

    virtual void
      computeRapidBuffers( const double &tcomp,
		   const int&nd,    const int &ncomp,
		   const int &nd1a, const int &nd1b,
		   const int &nd2a, const int &nd2b,
		   const int &nd3a, const int &nd3b,
		   const int &n1a,  const int &n1b,
		   const int &n2a,  const int &n2b,
		   const int &n3a,  const int &n3b,
		   const double &x, const double &q);

  };//end class GenericReaction

}; //end namespace CellWave
#endif
