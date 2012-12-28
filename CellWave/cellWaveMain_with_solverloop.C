#include <assert.h>
#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

#include "CellWave.h"
#include "getDiffusionDT.h"

int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  printf(" --------------------------------------------------------------------- \n");
  printf("   CellWave -- simulation of biochemical waves in cells\n");
  printf(" --------------------------------------------------------------------- \n");

  CellWave::Info data;

  //aString nameOfOGFile, nameOfShowFile;
  //cout << "example6>> Enter the name of the (old) overlapping grid file:" << endl;
  //cin >> nameOfOGFile;
  //cout << "example6>> Enter the name of the (new) show file (blank for none):" << endl;
  //cin >> nameOfShowFile;
  //real timeStepSize= 0.001;
  //int  totalNumberOfSteps=100;
  //int  saveEveryNthFrame =10;

  if (argc>4) {
    assert( argv[1] != NULL );
    assert( argv[2] != NULL );
    assert( argv[3] != NULL );
    assert( argv[4] != NULL );
    data.nameOfOGFile   = argv[1];
    data.nameOfShowFile = argv[2];
    sscanf( argv[3], "%le", &data.timeStepSize);
    sscanf( argv[4], "%d",  &data.totalNumberOfSteps);
    //data.tmax=timeStepSize*totalNumberOfSteps;
    if (argc>5)  sscanf( argv[5], "%d", data.saveEveryNthFrame);

    printf(".. got dt=%8.4e, num. steps=%d, save every %d frame, Tmax=%8.4e",
	   data.timeStepSize, data.totalNumberOfSteps, 
	   data.saveEveryNthFrame, data.tmax());
  } else {
    assert( argv[0] != NULL );
    printf("usage: %s <input grid> <output show file> <dt> <nsteps> <saveEvery, optional>\n", argv[0]);
    exit(-1);
  }

  // create and read in a CompositeGrid
  //CompositeGrid cg;
  getFromADataBase( data.cg, data.nameOfOGFile);
  data.cg.update();

  //Interpolant interpolant(cg);
  data.interpolant.updateToMatchGrid( data.cg );

  // ..create a show file
  // ..then save a general comment in the show file
  //Ogshow show( nameOfShowFile );                     
  data.show.open( data.nameOfShowFile );
  data.show.saveGeneralComment("alpha/CellWave for Calcium wave modeling");
  data.show.setFlushFrequency(100);   
  //CompositeGridOperators operators(cg);                        // operators for a CompositeGrid
  // operators.setOrderOfAccuracy(4);                          // for fourth order
  data.operators.updateToMatchGrid( data.cg );

  //
  // ..SETUP DATA
  //
  
  //....create grid functions
  CompositeGrid &cg = data.cg;

  int numberOfComponents=3; // p=[IP3], c=[Ca2+], h=inactivation, (b=buffer)
  Range all;
  realCompositeGridFunction q(cg,all,all,all,2*numberOfComponents);// create a grid function 
  q.setOperators( data.operators );                                 
  q.setName("all");              

  realCompositeGridFunction qn(cg,all,all,all,numberOfComponents);
  realCompositeGridFunction f(cg,all,all,all,numberOfComponents); //rhs in PDE
  realCompositeGridFunction fp(cg,all,all,all,numberOfComponents);//rhs prev step

  //CalciumChemistry chem(cg.numberOfDimensions());  // initializes parameters automatically (by constructor)
  CellWave::ReactionFactory rxnFactory;
  CellWave::GenericReaction *pChem  = rxnFactory.getReaction("LiRinzelWagner");
  CellWave::GenericReaction &chem = *pChem; // shorthand

  //..FIXME: won't work: need to set component numbers some otherway. 
  enum { cc=0, pc=1, hc=2,
	 bcfast=3, bcslow=4};
  //const int cc=chem.cc,                      
  //  pc=chem.pc, 
  //  hc=chem.hc, 
  //  bcfast=chem.bcfast, bcslow=  chem.bcslow;

  //..FIXME: main prog shouldn't know about the components

  //..FIXME --- make this work
#if 0
  int numComponents=chem.getNumberOfComponents();  //TODO
  q.setName( chem.getTitle() );
  for (int j=0; j<numComponents; ++j)  q.setName( chem.getComponentName( j ), j);
#endif 
 
  realCompositeGridFunction pconcentration, cconcentration, hconcentration;
  realCompositeGridFunction rhs_p, rhs_c, rhs_h; // components, and their RHS
  cconcentration.link(q,Range(cc,cc));  q.setName("Ca2+ (cytosol)",0); //FIXME
  pconcentration.link(q,Range(pc,pc));  q.setName("IP3",1);//FIXME
  hconcentration.link(q,Range(hc,hc));  q.setName("h Inhibitor",2);//FIXME

  rhs_c.link(f,Range(cc,cc));  f.setName("RHS for Ca2+",0);//FIXME
  rhs_p.link(f,Range(pc,pc));  f.setName("RHS for IP3",1);//FIXME
  rhs_h.link(f,Range(hc,hc));  f.setName("RHS for h Inhibitor",2);//FIXME

  //....set parameters
  //FIXME: parameters should be set from a parameter file.
  //FIXME: ALSO, hookup getDT to compute timestep
  //FIXME: what about the timestep from the reaction?? 

  //real tcomp=0, dt=.01;    
  //real tcomp=0, dt=.025; // works with D=300 um^2/sec 

  //real tcomp=0, dt=.00125;    // works with D=378 um^2/sec, dx_min= 5 um
  //int numberOfTimeSteps= 100000;
  //int saveEvery        = 200;

  //real tcomp=0, dt=.3;    // works with D=378 um^2/sec, dx_min= 5 um
  real tcomp=0, dt=.1;    // works with D=378 um^2/sec, dx_min= 5 um   //FIXME:
  int numberOfTimeSteps= 1000;                                         //FIXME
  int saveEvery        = 1;                                            //FIXME

  dt = timeStepSize;                      // from command line         //FIXME
  numberOfTimeSteps = totalNumberOfSteps; // from command line         //FIXME
  saveEvery         = saveEveryNthFrame;  // from command line         //FIXME

  realArray viscosity( chem.numberOfSpecies());
  viscosity( cc ) =  chem.getDiffusionCoefficient( int(cc) );
  viscosity( hc ) =  chem.getDiffusionCoefficient( int(hc) );
  viscosity( pc ) =  chem.getDiffusionCoefficient( int(pc) );

  real cfl =0.5;
  real dtDiffusion = getDiffusionDT(cfl, chem.numberOfSpecies(), viscosity, cg );
  //real dtDiffusion;
  dt = dtDiffusion; // or maybe dt = min(dt, dtDiffusion)
  printf("*** dt= %8.4e\n",  dt);

  //..set initial data!! --> allow hotspot in p,c,h (xc, yc, width)=exponential
  q=0.;
  
  for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
    //Index I1, I2,I3;
    //getIndex( cg[igrid].indexRange(), I1,I2,I3);

    MappedGrid & mg = cg[igrid];
    realArray & qArray = q[igrid];
    //realArray & dudtg = dudt[igrid];
    realArray & xArray = mg.vertex();  // array of vertices

    const IntegerArray & d = mg.dimension();
    const IntegerArray & gir= mg.gridIndexRange();
    const int nd   =cg.numberOfDimensions();
    const int ncomp=numberOfComponents;
    double tcomp=0.; //dummy time here
    
    //..explicit loop here: arguments are the info for a Fortran style loop
    chem.callInitialDataLoop(tcomp, 
			     nd, ncomp,
			     d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2), 
			     gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2), 
			     *xArray.getDataPointer(),*qArray.getDataPointer() );
    
  }; //end for igrid
  q.periodicUpdate();
  q.interpolate();

  //..start timestepping
  real wallTime=0.;
  int iWallTimeSteps=0;
  real totalWallTime=0.;
  char buffer[80];                                             // buffer for sprintf
  for( int istep=0; istep<=numberOfTimeSteps+1; istep++ )        // take some time steps
  {
    iWallTimeSteps++;
    real wallTimeThis=getCPU();

    if( istep % saveEvery == 0 )  // save solution every 'saveEvery' steps
    {
      double cmax=0., pmax=0., hmax=0.;
      for (int igrid=0; igrid < cg.numberOfComponentGrids(); ++igrid) {
	//where() ... add 'where mask>0' to get correct ip3 max 
	cmax = max(cmax, max(abs(cconcentration[igrid])) );
	pmax = max(pmax, max(abs(pconcentration[igrid])) );
	hmax = max(hmax, max(abs(hconcentration[igrid])) );
      }
      
      real avWallTime=wallTime/iWallTimeSteps;      
      
      printf("..step %5i (t=%6.2f): Max in c=%8.3e, p=%8.3e, h=%8.3e, CPU time/step %8.4e\n", 
	     istep, tcomp, cmax, pmax, hmax, avWallTime);
      wallTime =0.;     // reset
      iWallTimeSteps=0;

      show.startFrame();                                               // start a new frame
      show.saveComment(0,sPrintF(buffer,"Here is solution %i",istep)); // comment 0 (shown on plot)
      show.saveComment(1,sPrintF(buffer,"  t=%e ",tcomp));             // comment 1 (shown on plot)
      show.saveSolution( q );                                          // save the current grid function
    }

    //
    //..Perform one time step
    //
    //  (1) Chemistry: (Laplace terms imposed in step (2))
    //        dc/dt = beta*lam( (ip3 rel. )* (ca rel.) h^3 )(C_er - c) -nu_P *( )) + D*Laplace c
    //        dp/dt = Di*Laplace p - k_i *p 
    f = 0.;

    //..new style
    for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
      Index I1, I2,I3;
      getIndex( cg[igrid].gridIndexRange(), I1,I2,I3);
      
      MappedGrid & mg = cg[igrid];
      realArray & qArray   = q[igrid];
      realArray & rhsArray = f[igrid];
      //realArray & dudtg = dudt[igrid];
      realArray & xArray = mg.vertex();  // array of vertices
      
      const IntegerArray & d = mg.dimension();
      const IntegerArray & gir= mg.gridIndexRange();
      const int nd   =cg.numberOfDimensions();
      const int ncomp=numberOfComponents;

      //..explicit loop here: arguments are the info for a Fortran style loop
      chem.callRHSLoop(tcomp, 
		       nd, ncomp,
		       d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2), 
		       gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2), 
		       *xArray.getDataPointer(),*qArray.getDataPointer(),
		       *rhsArray.getDataPointer() );      

      //.. NEW TEST
      realArray pLaplace(I1,I2,I3);
      operators[igrid].derivative( MappedGridOperators::laplacianOperator, 
				   qArray, pLaplace, I1,I2,I3,pc);
      real dp= chem.getDiffusionCoefficient( int(pc) );
      rhsArray(I1,I2,I3,pc) += dp*pLaplace(I1,I2,I3); 

    }; //end for igrid

    //  (2) Diffusion
    real dp= chem.getDiffusionCoefficient( int(pc) );
    real dc= chem.getDiffusionCoefficient( int(cc) );

    //    rhs_p += dp *pconcentration.laplacian();  //didn't work
    rhs_c += dc *cconcentration.laplacian();  // move to chem loop
    rhs_h += 0.;

    //  (3) Update solution on next time (Adams-Bashforth, 2nd order, fully explicit) //FIXME -- adaptive timestepping
    if ( istep == 0 ) fp=f; // first step euler.
    //qn = q + .5*dt*( 3.*f - fp );
    //  (4) Update previous time levels
    //q=qn; 
    //fp=f;

    for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
      Index I1,I2,I3;
      const int lastcomp=hc;
      Range active(0, lastcomp);
      //getIndex( cg[igrid].gridIndexRange(), I1,I2,I3);
      
      qn[igrid](I1,I2,I3,active) = q[igrid](I1,I2,I3,active)
	+ .5*dt*( 3.*f[igrid](I1,I2,I3, active) - fp[igrid](I1,I2,I3, active) );

      //  (4) Update previous time levels      
      fp[igrid]=f[igrid];
      q[igrid](I1,I2,I3,active) = qn[igrid](I1,I2,I3,active);

      q[igrid](I1,I2,I3, cc+3) =  rhs_c[igrid](I1,I2,I3);
      q[igrid](I1,I2,I3, pc+3) =  rhs_p[igrid](I1,I2,I3);
      q[igrid](I1,I2,I3, hc+3) =  rhs_h[igrid](I1,I2,I3);

    }
    //  (5) Set bc's: internal (interpolate), and physical

    q.periodicUpdate();
    q.interpolate();                                           // interpolate

    // apply a neumann BC on all boundaries -- h does not get a bc
    //cconcentration.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    //pconcentration.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    q.applyBoundaryCondition(0,BCTypes::neumann,BCTypes::allBoundaries,0.);
    q.finishBoundaryConditions();

    // (6) update time to t_(n+1)
    tcomp +=dt;

    wallTimeThis = getCPU() - wallTimeThis;
    wallTime += wallTimeThis;
    totalWallTime += wallTimeThis;
  }
  Overture::finish();          

  cout << "----------------------------------------------------------------------\n";
  cout << "Computation completed after "<<  numberOfTimeSteps << " steps at computational T="
       << tcomp << endl;
  cout << "    Total CPU time ="<< totalWallTime 
       << " sec, average CPU time/step ="<< totalWallTime/numberOfTimeSteps<< " sec/step"<< endl;
  cout << "--> view the results by running:\n";
  cout << "     plotStuffx "<< nameOfShowFile <<endl;

  return 0;
    
}
