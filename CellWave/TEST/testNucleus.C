#include "Overture.h"
#include "Ogshow.h"
#include "CompositeGridOperators.h"
#include "PlotStuff.h"
#include "display.h"

#include "Nucleus.h"

// fortran names have an appended underscore:
#define mySolver mysolver_

extern "C"
{
  // fortran variables are passed by reference:
  void mySolver( const real &t, const real &dt,const real &a,const real&b,const real &nu, const int&nd,
     const int &nd1a,const int &nd1b,const int &nd2a,const int &nd2b,const int &nd3a,const int &nd3b,
     const int &n1a,const int &n1b,const int &n2a,const int &n2b,const int &n3a,const int &n3b,
               const real &x,const real &u, real &dudt );
}

int
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  printf("----------------------------------------------------------------------------\n");
  printf("   test 'class Nucleus'\n");
  printf("----------------------------------------------------------------------------\n");

  try {

  aString nameOfOGFile;
  if (argc>1) {
    nameOfOGFile = argv[1];
  }
  else {
    printf("usage: %s <overture grid file>\n",argv[0]);
    throw "error";
  }
  printf("..%s loading file %s...\n", argv[0], argv[1]);
 
  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update(MappedGrid::THEvertex | MappedGrid::THEmask);      // build vertices and mask
  Interpolant interpolant(cg);                                 // Make an interpolant
  CompositeGridOperators operators(cg);                        //operators for a CompositeGrid

  Range all;
  realCompositeGridFunction u(cg,all,all,all,1);               // create grid functions
  u.setOperators(operators);
  u.setName("Nucleus");                                        // name the grid function

  //..setup nucleus
  CellWave::Nucleus nucleus;
  nucleus.setID(0);
  CellWave::Nucleus::NucleusShape nucleusShape=CellWave::Nucleus::SphericalNucleus;
  //CellWave::Nucleus::NucleusShape nucleusShape=CellWave::Nucleus::BoxNucleus;
  nucleus.setShape( nucleusShape );
  nucleus.setBoundaryThickness( 20. );

  if( nucleusShape ==CellWave::Nucleus::SphericalNucleus) {
    nucleus.setCenter( -200,-100,0);
    //nucleus.setCenter( 0.,0.,0);
    nucleus.setRadius( 50.);
  }
  else if( nucleusShape ==CellWave::Nucleus::BoxNucleus) {
    nucleus.setCorners( -200,-100,0,
			-100, 0,  100);
  }

  bool openGraphicsWindow=TRUE;
  PlotStuff ps(openGraphicsWindow,"testNucleus");
  PlotStuffParameters psp;

  aString buff;                                              // buffer for sPrintF
  int numberOfTimeSteps=200;

  int grid;
  for( grid=0; grid<cg.numberOfComponentGrids(); grid++ ) // loop over grids
  {
    MappedGrid & mg = cg[grid];
    realArray & ug = u[grid];
    realArray & x = mg.vertex();  // array of vertices

    const IntegerArray & d = mg.dimension();
    const IntegerArray & gir= mg.gridIndexRange();
    const int nd=cg.numberOfDimensions();

    // call a fortran function to compute du/dt
    // (This function does not currently solve the convection diffusion equation)
    //mySolver( t,dt,a,b,nu,nd, d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2),
    //          gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
    //          *x.getDataPointer(),*ug.getDataPointer(), *dudtg.getDataPointer());
    const double tcomp=0.;
    const int ncomp=1;
    nucleus.getMaskArray(tcomp,nd, ncomp, d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2),
			 gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
			 *x.getDataPointer(),*ug.getDataPointer() );
    

    // display(ug,"ug","%6.3f");
    // display(dudtg,"dudtg","%6.3f");
  }

  psp.set(GI_PLOT_THE_OBJECT_AND_EXIT,false);
  psp.set(GI_TOP_LABEL,sPrintF(buff,"Nucleus"));  // set title
  ps.erase();
  PlotIt::contour(ps,u,psp);

  }
  catch( ... ) {
    printf("--error, exiting\n");
  }

  Overture::finish();
  return 0;

}



