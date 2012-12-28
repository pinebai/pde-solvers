#include "Overture.h"
#include "CompositeGrid.h"
#include "getDiffusionDT.h"


int main(int argc, char **argv) 
{
  Overture::start( argc, argv);

  printf("testing getDiffusionDT.\n");

  aString nameOfOGFile;
  double  nu=0.;
  int     numberOfComponents=1;

  if (argc>2) {
    assert( argv[1] != NULL );
    assert( argv[2] != NULL );
    nameOfOGFile   = argv[1];
    sscanf( argv[2], "%le", &nu);
  } else {
    assert( argv[0] != NULL );
    printf("usage: %s <input grid> <diffusion coeff>\n", argv[0]);
    exit(-1);
  }

  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update();

  real dt=0.;
  real cfl=1.;

  numberOfComponents=2;
  realArray nuArray(numberOfComponents);
  nuArray(0) = 300.;
  nuArray(1) = nu;
  realArray dtArray, dtGridArray;

  getDiffusionDT( cfl, numberOfComponents, nuArray, dtArray, dtGridArray, cg);

  printf("Maximum timestep:\n=================\n");
  for (int ic=0; ic<numberOfComponents; ++ic) {
    nu = nuArray(ic);
    dt = dtArray(ic);
    printf("    component %3d, D=%8.4e, dt_max=%8.4e\n", ic, nu, dt);
    //printf("    component %3d, D=%8.4e, dt_max=%8.4e\n", ic, nuArray(ic), dtArray(ic));
  }
  printf("\n");
  printf("Maximum timestep per grid:\n=========================\n");
  for (int ig=0; ig<cg.numberOfGrids(); ++ig) {
    printf("**  Grid %3d:\n", ig);
    for (int ic=0; ic<numberOfComponents; ++ic) {
      printf("    component %3d, D=%8.4e, dt_max=%8.4e\n", ic, nuArray(ic), dtGridArray(ig,ic));
    }
  }
  printf("\n");
  printf("--> maximum available timestep globally: dt< %8.4e\n", min( dtArray )); 
  real dtmax = getDiffusionDT( cfl, numberOfComponents, nuArray, cg);
  printf("--> simple call, max dt globally:        dt< %8.4e\n",dtmax); 

}
