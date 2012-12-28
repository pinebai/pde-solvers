//
// read & parse nucleus data files (.cwn)
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>
#include "NucleusGridFunction.h"

int main(int argc, char **argv)
{
  const std::string cn_filename="Grids/testNuclei.cwn";
  bool okFlag=false;

  CellWave::NucleusGridFunction nucleusFunc;
  nucleusFunc.readCellNucleusFile(cn_filename);
  nucleusFunc.printInfo();

  CellWave::NucleusGridFunction::IDVector nuclei;
  printf("Nuclei for each grid\n");
  for(int i=0; i<=6;  ++i ) {
    nucleusFunc.getGridNuclei( i, nuclei );
    printf(" nuclei.size() = %d  ", nuclei.size());
    if ( nuclei.size() > 0) {
      printf("..grid %d: ", i);
      for( int j=0; j< nuclei.size(); ++j ) {
	printf(" %d ", nuclei[j]);
      }
    }
    printf("\n");
  }
}
