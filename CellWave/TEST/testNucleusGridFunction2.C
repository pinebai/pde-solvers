//
// read & parse nucleus data files (.cwn)
//

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#include <iostream>

#include "Overture.h"
#include "PlotStuff.h"         

#include "NucleusGridFunction.h"

int main(int argc, char **argv)
{
  Overture::start(argc, argv);
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

  aString gridFileName="";
  if (argc>1) {
    gridFileName = argv[1];
  }
  else {
    printf("usage: %s <grid filename.hdf>\n", argv[0]);
    throw "error";
  }

  CompositeGrid cg;
  getFromADataBase( cg, gridFileName );
  cg.update();

  doubleCompositeGridFunction u( cg );
  u=1.;

  nucleusFunc.updateToMatchGrid( cg, u );
  nucleusFunc.evaluateGridFunction();
  //nucleusFunc.maskAGridFunction( u );

  bool openGraphicsWindow=TRUE;
  PlotStuff ps(openGraphicsWindow,"nucleus mask");  // create a PlotStuff object
  PlotStuffParameters psp;
  
  PlotIt::contour(ps,u,psp);  

  Overture::finish();
}
