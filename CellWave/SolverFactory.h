/// Brief description: Factory that creates solver instances
/// 
///   The main program in CellWave doesn't know anything about the
///   specific solver class being used. The main program simply
///   uses a 'GenericSolver', which in fact is a pointer to
///   a derived specific solver class.
///
///   The specific solver is instantiated by this class, SolverFactory.
///   Hence SolverFactory knows about all the different solver types.
///
#ifndef CELLWAVE_SOLVER_FACTORY_H
#define CELLWAVE_SOLVER_FACTORY_H "SolverFactory.h"

#include <math.h>
#include <stdio.h>
#include <iostream.h>
#include "GenericSolver.h"

namespace CellWave {
  
  class SolverFactory {
  public:
    SolverFactory();
    ~SolverFactory();

  private:  // do not use copy constructors.
    SolverFactory( const SolverFactory &X);
    SolverFactory &operator=( const SolverFactory  &X);

  public:
    GenericSolver *getSolver( const std::string &name, CellWave::Info &data );

  };
};
#endif
