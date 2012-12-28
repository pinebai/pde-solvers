//
// --- TEST CODE for overset grid X. Laevis waves
//     --- 3D grid/ timing

#include "Overture.h"
#include "Ogshow.h"  
#include "CompositeGridOperators.h"

//
// ..CalciumChemistry data, in dimensional units
//
class CalciumChemistry {
public:             // in Wagner, Pearson, Keizer:

  //..parameters, as Table 1, in Wagner, Pearson, Keizer
  real nu_L;       
  real d_I;
  real k_p;        
  real nu_P;       
  real I_s;        
  real C_er;       

  real tau_0;       // inhibition relax time
  real d_act;      

  real diffCalcium; // =D, micro m^2/sec     diffusion coeff/calcium

  real lambda;     
  real d_inh;      
  real beta;       

  // .... params not used by Wagner et al
  real diffIP3;     // =0                  diffusion coeff/IP3
  real k_i;         // =0 or 0.2 1/sec,    IP3 degradation

  
  //..initial data values
  real calcium_0;   // initial calcium2+ level
  real ip3_0;       // initial IP3 level
  real h_0;         // initial inhibitor level

  //..initial blob data
  real xBlob_c;
  real yBlob_c;
  real blobWidth_c;
  real blobMin_c, blobMax_c;

  real xBlob_p;
  real yBlob_p;
  real blobWidth_p;
  real blobMin_p, blobMax_p;

  // for inhomog distribution of IP3, Wagner et al, figure 4.
  enum IP3DistributionType {HomogeneousIP3, BlobIP3, 
	RadialIP3, RadialWithGradientIP3} ;

  IP3DistributionType ip3distribution;
  real ip3dist_I_s, ip3dist_I_h, ip3dist_r_c, ip3dist_I_w; 
  real ip3dist_xstar, ip3dist_Iprime_h, ip3dist_Iprime_w;

  CalciumChemistry() {
    // ..set parameters from Table 1, in Wagner Pearson, Keizer
    nu_L       = 5e-4;  // nondimensional
    d_I        = 0.025; // micro M
    k_p        = 0.4;   // micro M
    nu_P       = 0.1;   // micro M
    I_s        = 0.12;  // micro M
    C_er       = 10.;   // micro M
    tau_0      = 4;     // sec
    //tau_0      = 25;     // sec  --- DIFFERENT from Wagner
    d_act      = 1.2;   // micro M
    diffCalcium= 300.;  // micro m^2/sec
    //diffCalcium =0.; //TRY THIS
    lambda     = 112.5; // 1/sec
    d_inh      = 1.5;   // micro M
    beta       = 0.053; // nondimensional

    // ..parameters NOT USED BY WAGNER et al
    diffIP3    = 0.; //300.;   // IP3 diffusion
    k_i        = 0.;     // 1/sec, IP3 degradation, 0=synthetic IP3S

    //..initial data
    calcium_0  = 0.1153;   // micro M
    ip3_0      = 0.12;     // micro M
    h_0        = 0.43;     // nondimensional
    //    h_0        = 0.93;     // nondimensional

    //..initial blob data
    // .... for Ca2+
    blobMin_c = calcium_0; 
    blobMax_c = 2.;
    xBlob_c = -500.;
    yBlob_c =0.;
    blobWidth_c = 60;

    // .... for IP3
    blobMin_p = ip3_0;
    blobMax_p = ip3_0;  //20.;
    xBlob_p = -500.;
    yBlob_p =0.;
    blobWidth_p = 200;

    // .... for inhomog. distribution of IP3, fig. 4 in Wagner et al
    ip3dist_I_s = 0.12;  // u Molar
    ip3dist_I_h = 1.0;  
    ip3dist_I_w = 0.015;
    ip3dist_r_c = 500.; // u meter

    ip3dist_xstar    = 167;   // u meter
    ip3dist_Iprime_h = 0.84;  // u Molar
    ip3dist_Iprime_w = 0.8; 

    //ip3distribution = HomogeneousIP3;
    //ip3distribution = RadialIP3;
    ip3distribution = RadialWithGradientIP3;

    if (ip3distribution == RadialWithGradientIP3) {
      // CONCAVE WAVE PARAMS -- FIGURE 6 in Wagner et al -- redefine some params
      blobMax_c = calcium_0;  //2.;
      diffCalcium = 377.36; // Derr = 20um^2/sec
      lambda      = 250.;   // 1/sec 
    }
  }; 
};


int 
main(int argc, char *argv[])
{
  Overture::start(argc,argv);  // initialize Overture

  printf(" ---------------------------------------------------------------------------- \n");
  printf("Solve: Intracellular Ca2+/ IP3 wave dynamics on an overset grid\n");
  printf("       modeling Xenopus Laevis fertilization waves.\n");
  printf("Save results in a show file, use plotStuff to view this file                  \n");
  printf(" ---------------------------------------------------------------------------- \n");

  aString nameOfOGFile, nameOfShowFile;
  cout << "example6>> Enter the name of the (old) overlapping grid file:" << endl;
  cin >> nameOfOGFile;
  cout << "example6>> Enter the name of the (new) show file (blank for none):" << endl;
  cin >> nameOfShowFile;

  // create and read in a CompositeGrid
  CompositeGrid cg;
  getFromADataBase(cg,nameOfOGFile);
  cg.update();

  Interpolant interpolant(cg);                                 // Make an interpolant

  Ogshow show( nameOfShowFile );                               // create a show file
  show.saveGeneralComment("Overset grid X. Laevis Calcium waves");// save a general comment in the show file
  show.setFlushFrequency(1);                                  // flush file every 10 frames
    
  CompositeGridOperators operators(cg);                        // operators for a CompositeGrid
  // operators.setOrderOfAccuracy(4);                          // for fourth order

  //
  // ..SETUP DATA
  //
  
  //....create grid functions
  int numberOfDimensions= cg.numberOfDimensions();
  int numberOfComponents=3; // p=[IP3], c=[Ca2+], h=inactivation, (b=buffer)
  Range all;
  realCompositeGridFunction q(cg,all,all,all,numberOfComponents);// create a grid function 
  q.setOperators(operators);                                 
  q.setName("all");                                              // name the grid function

  realCompositeGridFunction qn(cg,all,all,all,numberOfComponents);
  realCompositeGridFunction f(cg,all,all,all,numberOfComponents); //rhs in PDE
  realCompositeGridFunction fp(cg,all,all,all,numberOfComponents);//rhs prev step

  realCompositeGridFunction pconcentration, cconcentration, hconcentration;
  realCompositeGridFunction rhs_p, rhs_c, rhs_h; // components, and their RHS
  cconcentration.link(q,Range(0,0));  q.setName("Ca2+ (cytosol)",0);
  pconcentration.link(q,Range(1,1));  q.setName("IP3",1);
  hconcentration.link(q,Range(2,2));  q.setName("h Inhibitor",2);

  rhs_c.link(f,Range(0,0));  f.setName("RHS for Ca2+",0);
  rhs_p.link(f,Range(1,1));  f.setName("RHS for IP3",1);
  rhs_h.link(f,Range(2,2));  f.setName("RHS for h Inhibitor",2);

  //....set parameters
  CalciumChemistry chem;    // initializes parameters automatically (by constructor)
  //real tcomp=0, dt=.01;    
  //real tcomp=0, dt=.025; // works with D=300 um^2/sec 

  //real tcomp=0, dt=.00125;    // works with D=378 um^2/sec, dx_min= 5 um
  //int numberOfTimeSteps= 100000;
  //int saveEvery        = 200;

  // for spherical-explicit-1mm-v2.hdf
  //real tcomp=0., dt=0.01;
  real tcomp=0., dt=0.1;
  int numberOfTimeSteps  =  1000;
  int saveEvery        =  10;
  //int saveEvery        =  1;

  // for spherical-explicit-1mm.hdf
  //real tcomp=0, dt=.1;    // works with D=378 um^2/sec, dx_min= 5 um
  //int numberOfTimeSteps= 3000;
  //int saveEvery        = 5;

  //..set initial data!! --> allow hotspot in p,c,h (xc, yc, width)=exponential
  q=1.;
  pconcentration=chem.calcium_0;
  cconcentration=chem.ip3_0;
  hconcentration=chem.h_0;

 
  for (int igrid=0; igrid< cg.numberOfComponentGrids(); ++igrid) {
    Index I1, I2,I3;
    getIndex( cg[igrid].indexRange(), I1,I2,I3);

    cconcentration[igrid](I1,I2,I3) = chem.calcium_0;
      
    if ( chem.ip3distribution == CalciumChemistry::HomogeneousIP3 ) {
      pconcentration[igrid](I1,I2,I3) = chem.ip3_0;
    }
    else if ( chem.ip3distribution == CalciumChemistry::BlobIP3 ) {
      if (2 == numberOfDimensions ) {
	pconcentration[igrid](I1,I2,I3) 
	  = chem.blobMin_p + ( chem.blobMax_p - chem.blobMin_p)
	  * exp( -( pow( cg[igrid].vertex()(I1,I2,I3,axis1)-chem.xBlob_p,2)
		    +  pow( cg[igrid].vertex()(I1,I2,I3,axis2)-chem.yBlob_p,2))/pow( chem.blobWidth_p,2.));
      } 
      else  if (3 == numberOfDimensions ) {
	pconcentration[igrid](I1,I2,I3) 
	  = chem.blobMin_p + ( chem.blobMax_p - chem.blobMin_p)
	  * exp( -( pow( cg[igrid].vertex()(I1,I2,I3,axis1)-chem.xBlob_p,2)
		    +  pow( cg[igrid].vertex()(I1,I2,I3,axis2)-chem.yBlob_p,2)
		    +  pow( cg[igrid].vertex()(I1,I2,I3,axis3)-chem.yBlob_p,2) )/pow( chem.blobWidth_p,2.));
      }
    }
    else if ( (chem.ip3distribution == CalciumChemistry::RadialIP3)
   ||  (chem.ip3distribution == CalciumChemistry::RadialWithGradientIP3)) {
      if (2 == numberOfDimensions ) {
	pconcentration[igrid](I1,I2,I3) //fig 4, Wagner et al
	  = chem.ip3dist_I_s *( 1. + chem.ip3dist_I_h
			* exp( (sqrt( pow( cg[igrid].vertex()(I1,I2,I3,axis1),2)
				    + pow( cg[igrid].vertex()(I1,I2,I3,axis2),2))/chem.ip3dist_r_c
				    -1.)/chem.ip3dist_I_w));  
      }
      else if (3 == numberOfDimensions ) {
	pconcentration[igrid](I1,I2,I3) //fig 4, Wagner et al
	  = chem.ip3dist_I_s *( 1. + chem.ip3dist_I_h
			* exp( (sqrt( pow( cg[igrid].vertex()(I1,I2,I3,axis1),2)
				    + pow( cg[igrid].vertex()(I1,I2,I3,axis2),2)
				    + pow( cg[igrid].vertex()(I1,I2,I3,axis3),2))/chem.ip3dist_r_c
				    -1.)/chem.ip3dist_I_w));  
      }
    } 
    else throw "error: unknown IP3 distribution type\n";

    if (chem.ip3distribution == CalciumChemistry::RadialWithGradientIP3) {
      where( cg[igrid].vertex()(I1,I2,I3,axis1) < -chem.ip3dist_xstar ) { // for x+x* <0
	if (numberOfDimensions == 2 ){
          pconcentration[igrid](I1,I2,I3) //fig 4, Wagner et al, second equation
            += chem.ip3dist_Iprime_h *(
               exp( (sqrt( pow( cg[igrid].vertex()(I1,I2,I3,axis1),2)
                    + pow( cg[igrid].vertex()(I1,I2,I3,axis2),2))/chem.ip3dist_r_c
                      -1.)/chem.ip3dist_I_w)
              * exp( - pow( cg[igrid].vertex()(I1,I2,I3,axis2)
                          / (chem.ip3dist_r_c*chem.ip3dist_Iprime_w), 4.)));

	}
	else if (numberOfDimensions ==3 ) {
        pconcentration[igrid](I1,I2,I3) //fig 4, Wagner et al, second equation
          += chem.ip3dist_Iprime_h *(
             exp( (sqrt( pow( cg[igrid].vertex()(I1,I2,I3,axis1),2)
                    + pow( cg[igrid].vertex()(I1,I2,I3,axis2),2)
		    + pow( cg[igrid].vertex()(I1,I2,I3,axis3),2))/chem.ip3dist_r_c
                      -1.)/chem.ip3dist_I_w)
            * exp( - pow( cg[igrid].vertex()(I1,I2,I3,axis2)
                        / (chem.ip3dist_r_c*chem.ip3dist_Iprime_w), 4.)));
	  
	}
      }//end where
    }//end if
  }; //end for
  q.periodicUpdate();
  q.interpolate();

  //..start timestepping
  char buffer[80];                                             // buffer for sprintf
  real wallTime=0.;
  int  iWallTimeSteps=0;
  for( int istep=0; istep<=numberOfTimeSteps+1; istep++ )        // take some time steps
  {

    double cmax=0., pmax=0., hmax=0., cmin=1e100, pmin=1e100, hmin=1e100;
    for (int igrid=0; igrid < cg.numberOfComponentGrids(); ++igrid) {
      cmax = max(cmax, max(abs(cconcentration[igrid])) );
      pmax = max(pmax, max(abs(pconcentration[igrid])) );
      hmax = max(hmax, max(abs(hconcentration[igrid])) );

      cmin = min(cmin, min(abs(cconcentration[igrid])) );
      pmin = min(pmin, min(abs(pconcentration[igrid])) );
      hmin = min(hmin, min(abs(hconcentration[igrid])) );
    }
    iWallTimeSteps++;
    double wallTimeThis=getCPU();
    real avWallTime=wallTime/iWallTimeSteps;
    printf("..step %5i (t=%6.2f): Max in c=%8.4e, p=%8.4e, h=%8.4e, CPU time/step=%8.4e sec\n", 
	   istep, tcomp, cmax, pmax, hmax,avWallTime);
    wallTime=0;     //reset
    iWallTimeSteps=0;

    if( istep % saveEvery == 0 )  // save solution every 'saveEvery' steps
    {
      printf(" -- saving frame --\n");
      show.startFrame();
      show.saveComment(0,sPrintF(buffer,"X. Laevis waves: solution %i",istep));
      show.saveComment(1,sPrintF(buffer,"  t=%10.6f ",tcomp));
      show.saveSolution( q ); 
    }

    //
    //..Perform one time step
    //
    //  (1) Chemistry: (Laplace terms imposed in step (2))
    //        dc/dt = beta*lam( (ip3 rel. )* (ca rel.) h^3 )(C_er - c) -nu_P *( )) + D*Laplace c
    //        dp/dt = Di*Laplace p - k_i *p 
    f = 0.;

    for( int igrid=0; igrid< cg.numberOfComponentGrids(); igrid++ ) {
      const CalciumChemistry &a=chem; // short hand
      const realArray &c      = cconcentration[igrid];
      const realArray &p      = pconcentration[igrid];
      const realArray &h      = hconcentration[igrid];
      realArray & rhsArray_c    = rhs_c[igrid]; // use array's to do the chemistry
      realArray & rhsArray_p    = rhs_p[igrid]; // --> these are symbolic links, and
      realArray & rhsArray_h    = rhs_h[igrid]; //     update the gridfcn automatically

      rhsArray_c = a.beta*a.lambda*( 
		      (a.nu_L + pow( p/( p+ a.d_I),3)*pow(c/(c+a.d_act),3.)*pow(h,3))
	              *( a.C_er - c ) 
		      - a.nu_P* pow(c,2)/( pow(c,2) + (a.k_p)*(a.k_p))); 

      rhsArray_p = - a.k_i *p;

      rhsArray_h = (a.d_inh - (c+ a.d_inh)* h)/a.tau_0;//Wagner etal
      //rhsArray_h = (a.d_inh/(c+ a.d_inh) - h)/a.tau_0; //causes pulses
    }      
    

    //  (2) Diffusion
    real dp= chem.diffIP3, dc= (chem.beta)*(chem.diffCalcium);

    rhs_p += dp *pconcentration.laplacian(); 
    rhs_c += dc *cconcentration.laplacian();
    rhs_h += 0.;

    //  (3) Update solution on next time (Adams-Bashforth, 2nd order, fully explicit)
    if ( istep == 0 ) fp=f; // first step euler.
    qn = q + .5*dt*( 3.*f - fp );
    //  (4) Update previous time levels
    q=qn; 
    fp=f;

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

  }

  Overture::finish();          
  return 0;
    
}
