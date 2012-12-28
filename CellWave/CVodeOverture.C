//
// CVodeOverture -- interfacing CVode in serial mode to Overture
//

#include "CVodeOverture.h"
#define Ith(v,i)    NV_Ith_S(v,i-1)       /* Ith numbers components 1..NEQ */ 

//--CVodeOvertureGrid -- a single mapped grid
CVodeOvertureGrid::
CVodeOvertureGrid() 
{ 
  cvodeData.absoluteTolerance=1e-6; 
  cvodeData.relativeTolerance=1e-4;
}

CVodeOvertureGrid::
~CVodeOvertureGrid()
{
  //default
}


void CVodeOvertureGrid::
initialize(int numberOfComponents, MappedGrid &mg, realArray &uSolution )
{    
  setMappedGrid( mg );
  setGridFunction( numberOfComponents, uSolution );
  //compute number of equations here
  
  CVodeData &cvd = cvodeData;
  const int neq = getNumberOfGridPoints();
  cvd.setNumberOfEquations( neq );
  cvd.machEnv   = M_EnvInit_Serial( neq );
  cvd.abstol    = N_VNew(neq, cvd.machEnv); assert( cvd.abstol != NULL );
  
  for(int i=0; i<neq; ++i ) {
    Ith( cvd.abstol, i) = cvd.absoluteTolerance;
  }

  const realtype time0=0.;
  cvd.cvode_mem = CVodeMalloc( neq, rhsFunction, time0, cvd.solutionVector,
			       cvd.getSolver(), cvd.getNonlinearSolver(),
			       SV, &cvd.relativeTolerance, cvd.abstol,
			       NULL, NULL, FALSE,
			       cvd.iopt, cvd.ropt, cvd.machEnv );
  cvd.solverFlag= CVSpgmr( cvd.cvode_mem, NONE, MODIFIED_GS, 0, 0.0, NULL, NULL,
			   NULL, NULL, NULL); // no precond/jac at this time
  if( cvd.solverFlag != SUCCESS ) {
    printf("CVSpgmr failed for grid %d.\n", getGridNumber());
  }
  
  //NOTES:
  //(1) CVode documentation lists N_VMAKE which doesn't exist. The calling seq. changed.
  //    Now we should use <new vec>= N_VMake( length, dataptr, machEnv)
}

void CVodeOvertureGrid::
setRHSFunction(   void (*rhs)(integertype N, 
			      realtype t, 
			      N_Vector y, 
			      N_Vector ydot,
			      void *f_data))
{
  rhsFunction = rhs; 
}

void CVodeOvertureGrid::
setMappedGrid(MappedGrid &mg )
{
  pMg = &mg;
  computeOvertureDimensions();
}

void CVodeOvertureGrid::
setGridFunction( int numberOfComponents, realArray &uSolution  )
{
  double *udata = getArrayDataPointer( uSolution );
  ncomponents   = numberOfComponents;
  const int neq = getNumberOfGridPoints();
  cvodeData.solutionVector= N_VMake( neq, udata,cvodeData.machEnv );
}


double *CVodeOvertureGrid::
getArrayDataPointer( const realArray &u ) // return Fortran array ptr to A++ realArray
{
  return( u.getDataPointer() );
}

void CVodeOvertureGrid::
computeOvertureDimensions( )
{
  assert( pMg != NULL );
  MappedGrid &mg = *pMg;
  const IntegerArray & d = mg.dimension();
  const IntegerArray & gir= mg.gridIndexRange();
  
  ndimensions=mg.numberOfDimensions();
  //nd1a,nd1b,nd2a,nd2b,nd3a,nd3b = d(0,0),d(1,0),d(0,1),d(1,1),d(0,2),d(1,2)    
  nd1a = d(0,0);
  nd1b = d(1,0);
  nd2a = d(0,1);
  nd2b = d(1,1);
  nd3a = d(0,2);
  nd3b = d(1,2);
  // n1a,n1b,n2a,n2b,n3a,n3b = gir(0,0),gir(1,0),gir(0,1),gir(1,1),gir(0,2),gir(1,2),
  n1a = gir(0,0);
  n1b = gir(1,0);
  n2a = gir(0,1);
  n2b = gir(1,1);
  n3a = gir(0,2);
  n3b = gir(1,2);
}

// ------------- CVodeOverture

CVodeOverture::
CVodeOverture( )
{
  numberOfGrids = 0;
  cvodeGrid = NULL;
}

CVodeOverture::
~CVodeOverture()
{
  delete cvodeGrid;
}

void CVodeOverture::
initialize( CompositeGrid &cg, realCompositeGridFunction &q)
{
  numberOfGrids= cg.numberOfComponentGrids();
  cvodeGrid = new CVodeOvertureGrid[ numberOfGrids ];
  
  numberOfComponents = q[0].numberOfComponents(); 
  
  for( int igrid=0; igrid< numberOfGrids; ++igrid ) {
    MappedGrid &mg = cg[igrid];
    realArray &u   = q[igrid];
    CVodeOvertureGrid &cvg= cvodeGrid[igrid];
    int ncomp      = q[igrid].numberOfComponents();

    if( ncomp != numberOfComponents ) {
      printf("CVodeOverture -- numberOfComponents[%d]=%d, ",igrid, ncomp);
      printf("differs from numberOfComponents[0]=%d\n",numberOfComponents);
    }
    
    cvg.initialize( ncomp, mg, u );

  }
}
