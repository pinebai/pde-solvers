/// Brief description: Factory that creates reaction chem instances
/// 


#include <math.h>
#include <stdio.h>
#include <iostream.h>
#include "GenericReaction.h"

#include "ReactionFactory.h"

//..known reactions
#include "ReactionLiRinzelWagner.h"
#include "Reaction2Buffer.h"
#include "ReactionSlepchenko2Buffer.h"

using namespace CellWave;

//
// .. constructor
//
ReactionFactory::
ReactionFactory()
{
  //.. do nothing
};

ReactionFactory::
~ReactionFactory()
{
  //.. do nothing
};


GenericReaction* ReactionFactory::
getReaction( const std::string &name) 
{
  //..NEW REACTIONS MUST BE ADDED HERE.
  GenericReaction *pRxn = NULL;
  if      ( name == ReactionLiRinzelWagner::rxnType() ) {
    pRxn = new ReactionLiRinzelWagner();
  } 
  else if ( name == Reaction2Buffer::rxnType() ) {
    pRxn = new Reaction2Buffer();
  }
  else if ( name == ReactionSlepchenko2Buffer::rxnType() ) {
    pRxn = new ReactionSlepchenko2Buffer();
  }
  else {
    // throw error
    throw "Reaction Factory error";
  };
} 

