#include <assert.h>
#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "CellWave.h"
//#include "getDiffusionDT.h"
#include "SolverLiRinzel.h"
#include "ReactionLiRinzelWagner.h"

int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  printf(" --------------------------------------------------------------------- \n");
  printf("   CellWave -- simulation of biochemical waves in cells\n");
  printf(" --------------------------------------------------------------------- \n");

  CellWave::Info data;

  if (argc>4) {
    assert( argv[1] != NULL );
    assert( argv[2] != NULL );
    assert( argv[3] != NULL );
    assert( argv[4] != NULL );
    data.nameOfOGFile   = argv[1];
    data.nameOfShowFile = argv[2];
    sscanf( argv[3], "%le", &data.timeStepSize);
    sscanf( argv[4], "%d",  &data.numberOfTimeSteps);
    //data.tmax=timeStepSize*totalNumberOfSteps;
    if (argc>5)  sscanf( argv[5], "%d", &data.saveEveryNthFrame);

    printf(".. got dt=%8.4e, num. steps=%d, save every %d frame, Tmax=%8.4e",
	   data.timeStepSize, data.numberOfTimeSteps, 
	   data.saveEveryNthFrame, data.maximumTime());
  } else {
    assert( argv[0] != NULL );
    printf("usage: %s <input grid> <output show file> <dt> <nsteps> <saveEvery, optional>\n", argv[0]);
    exit(-1);
  }

  //CellWave::ExplicitSolver solver( data );

  CellWave::SolverLiRinzel solver( data );
  solver.setup();
  solver.initialData();
  solver.solve();
  
  Overture::finish();          

  cout << "----------------------------------------------------------------\n";
  cout << "Computation completed after "<<  solver.getNumberOfTimeSteps ()
       << " steps at computational T="
       << solver.getTime() << endl;
  cout << "    Total CPU time ="<< solver.getTotalWallTime()
       << " sec, average CPU time/step =" 
       << solver.getTotalWallTime()/solver.getNumberOfTimeSteps()
       << " sec/step"<< endl;
  cout << "--> view the results by running:\n";
  cout << "     plotStuffx "<< data.nameOfShowFile <<endl;

  return 0;
    
}
