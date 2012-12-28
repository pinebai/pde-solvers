#include "Overture.h"
#include <stdio.h>
#include "InterpGridArray.h"

void InterpGridArray::freeMemory() {
  if( this->isAllocated() ) {
    delete [] pInterp;
  }
}

bool InterpGridArray::isAllocated() { return pInterp!=NULL; };

void InterpGridArray::resize( const GridAxisSide &graxsi )
{
  this->freeMemory();
  dimensions = graxsi;
  int size=dimensions.getSize();
  pInterp = new InterpolatePoints[ size ];
}

InterpolatePoints& InterpGridArray::get( const GridAxisSide &graxsi)
{
  int offset=graxsi.getOffset( this->dimensions );
  printf( " offset = %4d; ",offset);
  return pInterp[offset];
}

