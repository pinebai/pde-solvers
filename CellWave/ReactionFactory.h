/// Brief description: Factory that creates reaction chem instances
/// 
///   The main program in CellWave doesn't know anything about the
///   specific rxn mechanism being used. The main program simply
///   uses a 'GenericReaction', which in fact is a pointer to
///   a derived specific reaction class.
///
///   The specific reaction is instantiated by this class, ReactionFactory.
///   Hence ReactionFactory knows about all the different reaction types.
///
#ifndef CELLWAVE_REACTION_FACTORY_H
#define CELLWAVE_REACTION_FACTORY_H "ReactionFactory.h"

#include <math.h>
#include <stdio.h>
#include <iostream.h>
#include "GenericReaction.h"

namespace CellWave {
  
  class ReactionFactory {
  public:
    ReactionFactory();
    ~ReactionFactory();

  private: // do not use copy constructors
    ReactionFactory( const ReactionFactory &X);
    ReactionFactory &operator=( const ReactionFactory &X);

  public:

    GenericReaction *getReaction( const std::string &name);
    
    // need: query functions
  }; // end class ReactionFactory
}; // end namespace CellWave
#endif
