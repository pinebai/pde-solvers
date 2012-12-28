//
// Who to blame:  Geoff Chesshire
//

#include "CompositeGrid.h"
#include "display.h"

#ifdef USE_STL
RCVector_STATIC_MEMBER_DATA(CompositeGrid)
#endif // USE_STL

//  Define a triple for loop.  The macro is defined only within this file.
#define COMPOSITE_GRID_FOR_3(range,i,j,k)                              \
        for (k=((Integer*)(range))[4]; k<=((Integer*)(range))[5]; k++) \
        for (j=((Integer*)(range))[2]; j<=((Integer*)(range))[3]; j++) \
        for (i=((Integer*)(range))[0]; i<=((Integer*)(range))[1]; i++)

//
// class CompositeGrid:
//
// Public member functions:
//
// Default constructor.
//
// If numberOfDimensions_==0 (e.g., by default) then create a null
// CompositeGrid.  Otherwise create a CompositeGrid
// with the given number of dimensions and number of component grids.
//
CompositeGrid::CompositeGrid(
  const Integer numberOfDimensions_,
  const Integer numberOfComponentGrids_):
  GridCollection() {
    className = "CompositeGrid";
    master=this;
    rcData = new
      CompositeGridData(numberOfDimensions_, numberOfComponentGrids_);
    isCounted = LogicalTrue;
    rcData->incrementReferenceCount();
    updateReferences();
}
//
// Copy constructor.  (Does a deep copy by default.)
//
CompositeGrid::CompositeGrid(
  const CompositeGrid& x,
  const CopyType       ct):
  GridCollection() {
    className = "CompositeGrid";
    master=this;
    switch (ct) {
      case DEEP:
      case NOCOPY:
        rcData = (CompositeGridData*)
          ((ReferenceCounting*)x.rcData)->virtualConstructor(ct);
        isCounted = LogicalTrue;
        rcData->incrementReferenceCount();
      break;
      case SHALLOW:
        rcData = x.rcData;
        isCounted = x.isCounted;
        if (isCounted) rcData->incrementReferenceCount();
      break;
    } // end switch
    updateReferences();
}
//
// Destructor.
//
CompositeGrid::~CompositeGrid()
  { if (isCounted && rcData->decrementReferenceCount() == 0) delete rcData; }
//
// Assignment operator.  (Does a deep copy.)
//
CompositeGrid& CompositeGrid::
operator=(const CompositeGrid& x) 
{
//  GridCollection::operator=(x);
  if (rcData != x.rcData) {
    if (rcData->getClassName() == x.rcData->getClassName())
    {
      (ReferenceCounting&)*rcData = (ReferenceCounting&)*x.rcData;
      updateReferences();

      // *wdh* 000612 : to get refinementLevels etc.
      if( rcData->computedGeometry & THEbaseGrid )
        rcData->update(THEbaseGrid); 
      if( rcData->computedGeometry & THEcomponentGrid )
        rcData->update(THEcomponentGrid); 
      if( rcData->computedGeometry & THErefinementLevel )
        rcData->update(THErefinementLevel); 
      if( rcData->computedGeometry & THEmultigridLevel)
        rcData->update(THEmultigridLevel);
    } 
    else
    {
      CompositeGrid& y =
	*(CompositeGrid*)x.virtualConstructor();
      reference(y); delete &y;
    } // end if
    master=x.master;
  } // end if
  return *this;
}

//\begin{>>CompositeGridInclude.tex}{\subsubsection{reference(CompositeGrid)}} 
void CompositeGrid::
reference(const CompositeGrid& x) 
// ===========================================================
// /Description:
//    Make a reference.  (Does a shallow copy.)
//\end{CompositeGridInclude.tex}
// ===========================================================
{
  GridCollection::reference(x);
  if (rcData != x.rcData) 
  {
    if (isCounted && rcData->decrementReferenceCount() == 0)
      delete rcData;
    rcData = x.rcData;
    isCounted = x.isCounted;
    if (isCounted) rcData->incrementReferenceCount();
    // *wdh* updateReferences();
    master=x.master;
  } // end if
  updateReferences();   // *wdh* 000322 -- we must always do this since the number of grids etc. may have changed.

}

void CompositeGrid::
reference(CompositeGridData& x) 
{
  GridCollection::reference(x);
  if (rcData != &x) 
  {
    if (isCounted && rcData->decrementReferenceCount() == 0)
      delete rcData;
    rcData = &x;
    isCounted = !x.uncountedReferencesMayExist();
    if (isCounted) rcData->incrementReferenceCount();
    // *wdh* updateReferences();
  } // end if
  updateReferences();   // *wdh* 000322 -- we must always do this since the number of grids etc. may have changed.
}

//
// Break a reference.  (Replaces with a deep copy.)
//
void CompositeGrid::breakReference() {
//  GridCollection::breakReference();
    if (!isCounted || rcData->getReferenceCount() != 1) {
        CompositeGrid x = *this; // Uses the (deep) copy constructor.
        reference(x);
    } // end if
}
//
// Change the grid to be all vertex-centered.
//
void CompositeGrid::changeToAllVertexCentered() {
    Logical isAllVertexCentered = LogicalTrue;
    for (Integer i=0; i<numberOfGrids(); i++)
      isAllVertexCentered =
      isAllVertexCentered && grid[i].isAllVertexCentered();
    if (!isAllVertexCentered)
      destroy(EVERYTHING & ~GridCollection::EVERYTHING);
    GridCollection::changeToAllVertexCentered();
}
//
// Change the grid to be all cell-centered.
//
void CompositeGrid::changeToAllCellCentered() {
    Logical isAllCellCentered = LogicalTrue;
    for (Integer i=0; i<numberOfGrids(); i++)
      isAllCellCentered =
      isAllCellCentered && grid[i].isAllCellCentered();
    if (!isAllCellCentered)
      destroy(EVERYTHING & ~GridCollection::EVERYTHING);
    GridCollection::changeToAllCellCentered();
}

//\begin{>>CompositeGridInclude.tex}{\subsubsection{changeInterpolationWidth}} 
int CompositeGrid::
changeInterpolationWidth( int width )
// ===========================================================
// /Description:
//   Reduce interpolation width of an already computed overlapping grid.
// For example you may want to use 2-point interpolation instead of 3-point
// interpolation. This routine will adjust the interpoleeLocation array
// and interpolationWidth array. It is currently not possible to reset the
// grid to it's original width (although this could be supported)
//.
// /width (input) : a positive integer. The interpolation with can only be
// decreased.
//\end{CompositeGridInclude.tex}
// ===========================================================
{
  assert( width>=1 );
  
  const IntegerArray & iw0 = interpolationWidth();
  IntegerArray & iw = (IntegerArray &) iw0;
  
  int g;
  for( g=0; g<numberOfComponentGrids(); g++ )
  {
    if( numberOfInterpolationPoints(g)>0 )
    {
      intArray & vWidth = variableInterpolationWidth[g];
      intArray & il = interpoleeLocation[g];
      intArray & ig = interpoleeGrid[g];
      const realArray & ic = interpolationCoordinates[g];

      int indexPosition;
      real relativeOffset,px;
      for( int i=0; i<numberOfInterpolationPoints(g); i++ )
      {
	int oldWidth=vWidth(i);
	if( width < oldWidth )
	{
	  int gridi = ig(i);   // **** could vectorize this loop since list is sorted by interpolee
	  MappedGrid & cgridi = (*this)[gridi];

	  for( int axis=axis1; axis<numberOfDimensions(); axis++ ) 
	  {
	    indexPosition=il(i,axis);
	    relativeOffset=ic(i,axis)/cgridi.gridSpacing(axis)+cgridi.indexRange(Start,axis);
	    // for 3-pt interpolation : 0<= px <=3 and normally .5<= px <= 2.5 if centred.
	    px= cgridi.isCellCentered(axis)  ? relativeOffset-indexPosition-.5  : relativeOffset-indexPosition;

	    //......interpolation width less than maximum allowed
	    if( px > width/2. )
	    {
	      // we need to increase the interpoleeLocation
	      int ipx=min(int(px-(width-2)/2.),oldWidth-width);
	      il(i,axis)+=ipx;

	      // printf("grid=%i, i=%i gridi=%i, axis=%i px=%8.2e shift=%i\n",g,i,gridi,axis,px,ipx);
	    }
	  }
	}
    
      }
      vWidth=min(vWidth,width);
    }
  } // end for g
  iw=min(iw,width);
  
  return 0;
}



//
// Check that the data structure is self-consistent.
//
void CompositeGrid::consistencyCheck() const {
    GridCollection::consistencyCheck();
    if (rcData != GridCollection::rcData) {
        cerr << className << "::consistencyCheck():  "
             << "rcData != GridCollection::rcData for "
             << getClassName() << " " << getGlobalID() << "." << endl;
        assert(rcData == GridCollection::rcData);
    } // end if
    numberOfInterpolationPoints      .Test_Consistency();
    numberOfImplicitInterpolationPoints.Test_Consistency();
    interpolationStartEndIndex.Test_Consistency();

//    numberOfInterpoleePoints         .Test_Consistency();
    interpolationIsImplicit          .Test_Consistency();
//    backupInterpolationIsImplicit    .Test_Consistency();
    interpolationWidth               .Test_Consistency();
//    backupInterpolationWidth         .Test_Consistency();
    interpolationOverlap             .Test_Consistency();
    maximumHoleCuttingDistance       .Test_Consistency();
//    backupInterpolationOverlap       .Test_Consistency();
//    interpolationConditionLimit      .Test_Consistency();
//    backupInterpolationConditionLimit.Test_Consistency();
    interpolationPreference          .Test_Consistency();
    mayInterpolate                   .Test_Consistency();
//    mayBackupInterpolate             .Test_Consistency();
    mayCutHoles                      .Test_Consistency();
    sharedSidesMayCutHoles           .Test_Consistency();
    multigridCoarseningRatio         .Test_Consistency();
    multigridProlongationWidth       .Test_Consistency();
    multigridRestrictionWidth        .Test_Consistency();
//    interpoleeGridRange              .Test_Consistency();
    interpolationCoordinates         .consistencyCheck();
    interpoleeGrid                   .consistencyCheck();
    variableInterpolationWidth       .consistencyCheck();
//    interpoleePoint                  .consistencyCheck();
    interpoleeLocation               .consistencyCheck();
    interpolationPoint               .consistencyCheck();
//    interpolationCondition           .consistencyCheck();
    multigridLevel                   .consistencyCheck();
//    inverseCondition                 .consistencyCheck();
    inverseCoordinates               .consistencyCheck();
    inverseGrid                      .consistencyCheck();
}
// Here is the master grid.
CompositeGrid & CompositeGrid::
masterGridCollection()
{
  assert( master!=0 );
  return *master;
}


//
// "Get" and "put" database operations.
//
Integer CompositeGrid::get(
  const GenericDataBase& db,
  const aString&         name) {
    Integer returnValue = rcData->get(db, name);
    updateReferences();
    return returnValue;
}
Integer CompositeGrid::put(
  GenericDataBase& db,
  const aString&   name) const
  { return rcData->put(db, name); }
//
// Set references to reference-counted data.
//
void CompositeGrid::updateReferences(const Integer what) {
    GridCollection::reference(*rcData);
#define REFERENCE(x) x.reference(rcData->x)
#define REF_ARRAY(x) \
    if (x.getDataPointer() != rcData->x.getDataPointer()) REFERENCE(x)
    REF_ARRAY(numberOfInterpolationPoints);
    REF_ARRAY(numberOfImplicitInterpolationPoints);
    REF_ARRAY(interpolationStartEndIndex);

//    REF_ARRAY(numberOfInterpoleePoints);
    REF_ARRAY(interpolationIsImplicit);
//    REF_ARRAY(backupInterpolationIsImplicit);
    REF_ARRAY(interpolationWidth);
//    REF_ARRAY(backupInterpolationWidth);
    REF_ARRAY(interpolationOverlap);
    REF_ARRAY(maximumHoleCuttingDistance);
   
//    REF_ARRAY(backupInterpolationOverlap);
//    REF_ARRAY(interpolationConditionLimit);
//    REF_ARRAY(backupInterpolationConditionLimit);
    REF_ARRAY(interpolationPreference);
    REF_ARRAY(mayInterpolate);
//    REF_ARRAY(mayBackupInterpolate);
    REF_ARRAY(mayCutHoles);
    REF_ARRAY(sharedSidesMayCutHoles);
    REF_ARRAY(multigridCoarseningRatio);
    REF_ARRAY(multigridProlongationWidth);
    REF_ARRAY(multigridRestrictionWidth);
//    REF_ARRAY(interpoleeGridRange);
    REFERENCE(interpolationCoordinates);
    REFERENCE(interpoleeGrid);
    REFERENCE(variableInterpolationWidth);
//    REFERENCE(interpoleePoint);
    REFERENCE(interpoleeLocation);
    REFERENCE(interpolationPoint);
//    REFERENCE(interpolationCondition);
    REFERENCE(multigridLevel);
//    REFERENCE(inverseCondition);
    REFERENCE(inverseCoordinates);
    REFERENCE(inverseGrid);
#undef REFERENCE
#undef REF_ARRAY
#define SET_GRID(x) if (x.gridCollectionData == rcData) x.gridCollection = this
//    SET_GRID(inverseCondition);
    SET_GRID(inverseCoordinates);
    SET_GRID(inverseGrid);
#undef SET_GRID
    GridCollection::updateReferences(what);

    int i;
#ifdef USE_STL
/* is this correct? */
#define FOR_COLLECTION(X) \
  for( i=list.begin(); i<=list.begin(); i++ ) \
     X[i].master=this;
#else
#define FOR_COLLECTION(X) \
  for( i=0; i<X.getLength(); i++ ) \
    X[i].master=this;
#endif
  FOR_COLLECTION(multigridLevel);
#undef FOR_COLLECTION
}
//
// Update the grid, sharing the data of another grid.
//
Integer CompositeGrid::update(
  GenericGridCollection& x,
  const Integer          what,
  const Integer          how)
{
  Integer upd = rcData->update(*((CompositeGrid&)x).rcData, what, how);
  updateReferences(what);

  // We need to assign the mappedGrid in each geometry array *wdh*
  CompositeGrid & cg = (CompositeGrid&)x;
  for( int k=0; k<cg.numberOfGrids(); k++ )
    cg[k].updateMappedGridPointers(what);
    
  return upd;
}
//
// Destroy optional grid data.
//
void CompositeGrid::destroy(const Integer what) {
    rcData->destroy(what);
    updateReferences();
}


//\begin{>>CompositeGridInclude.tex}{\subsubsection{add(Mapping)}}
int CompositeGrid::
add(Mapping & map)
// ==========================================================================
// /Description: 
//   Add a new grid, built from a Mapping
//\end{CompositeGridInclude.tex}
//==========================================================================
{
  MappedGrid g(map);
  if( numberOfGrids()>0 && g.getGridType()!=GenericGrid::unstructuredGrid )
  {
    // set number of ghost points equal to that from the previous grid.
    MappedGrid & mg = (*this)[numberOfGrids()-1];
    for( int axis=0; axis<mg.numberOfDimensions(); axis++ )
      for( int side=Start; side<=End; side++ )
      {
	g.setNumberOfGhostPoints(side,axis,mg.numberOfGhostPoints(side,axis));
        // printf("CompositeGrid::add: set ghsot points to %i\n",mg.numberOfGhostPoints(side,axis));
      }
    
  }
  return add(g);
}
    
//\begin{>>CompositeGridInclude.tex}{\subsubsection{add(MappedGrid)}}
int CompositeGrid::
add(MappedGrid & g)
// ==========================================================================
// /Description: 
//    Add a new grid. The grid collection will keep a reference to g.
//\end{CompositeGridInclude.tex}
//==========================================================================
{
  int returnValue;
  const int numberOfOldGrids=numberOfGrids();
  
  returnValue = GridCollection::add(g);
  updateReferences();

  // update some of the parameter arrays:
  // interpolationWidth

  const int l = 0; // mg level
  const int n=numberOfOldGrids;
  Range G=numberOfOldGrids, all;
  
  numberOfInterpolationPoints(n) = 0;
  if( numberOfOldGrids>0 )
  {
    const int n0=n-1; // assign new values from this grid, choose the last one
    rcData->interpolationIsImplicit(G,n,l)= rcData->interpolationIsImplicit(G,n0,l);
    rcData->interpolationIsImplicit(n,G,l)= rcData->interpolationIsImplicit(n0,G,l);
    rcData->interpolationIsImplicit(n,n,l)= rcData->interpolationIsImplicit(n0,n0,l);

    for (int axis=0; axis<3; axis++) 
    {
      rcData->interpolationWidth(axis,G,n,l) = rcData->interpolationWidth(axis,G,n0,l);
      rcData->interpolationWidth(axis,n,G,l) = rcData->interpolationWidth(axis,n0,G,l);
      rcData->interpolationWidth(axis,n,n,l) = rcData->interpolationWidth(axis,n0,n0,l);


      rcData->interpolationOverlap(axis,G,n,l) = rcData->interpolationOverlap(axis,G,n0,l);
      rcData->interpolationOverlap(axis,n,G,l) = rcData->interpolationOverlap(axis,n0,G,l);
      rcData->interpolationOverlap(axis,n,n,l) = rcData->interpolationOverlap(axis,n0,n0,l);
    
      rcData->multigridCoarseningRatio(axis,n,l) = rcData->multigridCoarseningRatio(axis,n0,l);
      rcData->multigridProlongationWidth(axis,n,l) = rcData->multigridProlongationWidth(axis,n0,l);
      rcData->multigridRestrictionWidth(axis,n,l) = rcData->multigridRestrictionWidth(axis,n0,l);

    } 
    rcData->maximumHoleCuttingDistance(all,all,n)=rcData->maximumHoleCuttingDistance(all,all,n0);
  }
  else
  {
    // this is the first grid in the collection, use default values:
    rcData->interpolationIsImplicit(G,n,l)= true;
    rcData->interpolationIsImplicit(n,G,l)= true;
    rcData->interpolationIsImplicit(n,n,l)= true;

    for (int axis=0; axis<3; axis++) 
    {
      rcData->interpolationWidth(axis,G,n,l) = 3;
      rcData->interpolationWidth(axis,n,G,l) = 3;
      rcData->interpolationWidth(axis,n,n,l) = 3;


      rcData->interpolationOverlap(axis,G,n,l) = 0.;
      rcData->interpolationOverlap(axis,n,G,l) = 0.;
      rcData->interpolationOverlap(axis,n,n,l) = 0.;
    
      rcData->multigridCoarseningRatio(axis,n,l) = 2;
      rcData->multigridProlongationWidth(axis,n,l) = 3;
      rcData->multigridRestrictionWidth(axis,n,l) = 3;

    } 
    rcData->maximumHoleCuttingDistance(all,all,n)=SQRT(.1*REAL_MAX);
  }

  rcData->mayInterpolate(n,G,l) = true;
  rcData->mayInterpolate(G,n,l) = true;
  rcData->mayInterpolate(n,n,l) = true;

  rcData->interpolationPreference(n,G,l) = -1; // what should this be ?
  rcData->interpolationPreference(G,n,l) = -1;
  rcData->interpolationPreference(n,n,l) = -1;

  rcData->mayCutHoles(n,G)=true;
  rcData->mayCutHoles(G,n)=true;

  rcData->sharedSidesMayCutHoles(n,G)=false;
  rcData->sharedSidesMayCutHoles(G,n)=false;
  

  return returnValue;
}


    // delete a grid:
int CompositeGrid::
deleteGrid(Integer k)
{
  IntegerArray gridsToDelete(1);
  gridsToDelete=k;

  return deleteGrid(gridsToDelete);

}

// delete a list of grids:
int CompositeGrid::
deleteGrid(const IntegerArray & gridsToDelete )
{
  // see GenericGridCollection::deleteMultigridLevels
  int numberToDelete=gridsToDelete.getLength(0);
  if( numberToDelete==0 )
  {
    printf("CompositeGrid::deleteGrid:WARNING: no grids were specified to be deleted\n");
    return -1;
  }
  int newNumberOfGrids=numberOfGrids()-numberToDelete;
  if( newNumberOfGrids<0 )
  {
    printf("CompositeGrid::deleteGrid:ERROR:trying to delete %i grids, but the collection only has %i grids\n",
	   numberToDelete,numberOfGrids());
    return 1;
  }
  if( min(gridsToDelete)<0 || max(gridsToDelete)>=numberOfGrids() )
  {
    printf("CompositeGrid::deleteGrid:ERROR:trying to delete an invalid grid. There are %i grids\n",numberOfGrids());
    gridsToDelete.display("gridsToDelete");
    Overture::abort("error");
  }

  display(gridsToDelete,"CompositeGrid::deleteGrid: gridsToDelete");

  // make a list of grids that remain:
  IntegerArray save(numberOfGrids()),ia,gridsToSave(numberOfGrids()-numberToDelete);
  save.seqAdd(0,1);        // all grids: 0,1,2,3,...
  save(gridsToDelete)=-1;  // mark deleted grids
  ia=(save>=0).indexMap(); // list of remaining grids
  gridsToSave=save(ia);    // compressed list of remaining grids

  display(gridsToSave,"gridsToSave");

  int returnValue=deleteGrid(gridsToDelete,gridsToSave);

  return returnValue;
}

// delete a list of grids (use this when you also know the list of grids to save).
int CompositeGrid::
deleteGrid(const IntegerArray & gridsToDelete, const IntegerArray & gridsToSave )
// should this function be protected?
{
  int returnValue=0;
  // shift values before deleting the grids.

  // we need to re-assign values in the arrays found in setNumberOfDimensionsAndGrids
  // before they are reshaped to be a smaller size. 
  // Note that a(3,4) and b= a.reshape(2,3) will satisfy a(0:1,0:2)==b(0:1,0:2)
  //
  int numberToDelete=gridsToDelete.getLength(0);
  int newNumberOfGrids=gridsToSave.getLength(0);

  if( min(gridsToDelete) > max(gridsToSave) )
  {
    // not need to shift parameters values since grids to be deleted are on
    // the end of the list.
  }
  else
  { 
    // shift parameters around in preparation for grids to be deleted.
    // These arrays will be resized below.

    Range R=newNumberOfGrids;

//     // do this until A++ is fixed.
//     IntegerArray temp(R,R);
//     for( int k=0; k<newNumberOfGrids; k++ )
//       temp(R,k)=mayCutHoles(gridsToSave,gridsToSave(k));
//     mayCutHoles(R,R)=temp;


    numberOfInterpolationPoints(R)=numberOfInterpolationPoints(gridsToSave)*1;
    numberOfImplicitInterpolationPoints(R)=numberOfImplicitInterpolationPoints(gridsToSave)*1;
    int i,j,k;
    for( i=0; i<4; i++ )
      for( j=0; j<newNumberOfGrids; j++ )
      for( k=0; k<newNumberOfGrids; k++ )
        interpolationStartEndIndex(i,j,k)=interpolationStartEndIndex(i,gridsToSave(j),gridsToSave(k));

    for( i=0; i<3; i++ )
      for( j=0; j<2; j++ )
      for( k=0; k<newNumberOfGrids; k++ )
	maximumHoleCuttingDistance(i,j,k)=maximumHoleCuttingDistance(i,j,gridsToSave(k));       


    for( j=0; j<newNumberOfGrids; j++ )
      for( k=0; k<newNumberOfGrids; k++ )
      {
	mayCutHoles(j,k)=mayCutHoles(gridsToSave(j),gridsToSave(k));
	sharedSidesMayCutHoles(j,k)=sharedSidesMayCutHoles(gridsToSave(j),gridsToSave(k));
      }
    
    for( int l=0; l<numberOfMultigridLevels(); l++ )
    {
      for( j=0; j<newNumberOfGrids; j++ )
      for( k=0; k<newNumberOfGrids; k++ )
        interpolationIsImplicit(j,k,l)=interpolationIsImplicit(gridsToSave(j),gridsToSave(k),l);
      for( i=0; i<3; i++ )
      {
        for( j=0; j<newNumberOfGrids; j++ )
        for( k=0; k<newNumberOfGrids; k++ )
	{
	  interpolationWidth(i,j,k,l)=interpolationWidth(i,gridsToSave(j),gridsToSave(k),l);
	  interpolationOverlap(i,j,k,l)=interpolationOverlap(i,gridsToSave(j),gridsToSave(k),l);
	}
	for( j=0; j<newNumberOfGrids; j++ )
	{
	  multigridCoarseningRatio(i,j,l)=multigridCoarseningRatio(i,gridsToSave(j),l);
	  multigridProlongationWidth(i,j,l)=multigridProlongationWidth(i,gridsToSave(j),l);
	  multigridRestrictionWidth (i,j,l)=multigridRestrictionWidth(i,gridsToSave(j),l);
	}
      }
      for( j=0; j<newNumberOfGrids; j++ )
      for( k=0; k<newNumberOfGrids; k++ )
      {
	interpolationPreference(j,k,l)=interpolationPreference(gridsToSave(j),gridsToSave(k),l);
	mayInterpolate(j,k,l)=mayInterpolate(gridsToSave(j),gridsToSave(k),l);
      }
      
    }

    // printf("CompositeGrid::deleteGrid:WARNING: This case is not finished\n");
  }

  returnValue=GridCollection::deleteGrid(gridsToDelete,gridsToSave); // this will resize parameter arrays.
  updateReferences();

  printf("CompositeGrid::deleteGrids: newNumberOfGrids=%i, numberOfGrids=%i, numberOfComponentGrids=%i \n",
        newNumberOfGrids,numberOfGrids(),numberOfComponentGrids());
  

  return returnValue;
}




//
// Add a refinement grid to the collection.
//
Integer CompositeGrid::
addRefinement(const IntegerArray& range,  // The indexRange of the refinement grid.
	      const IntegerArray& factor, // The refinement factor w.r.t. level-1.
	      const Integer&      level,  // The refinement level number.
	      const Integer       k)     // The index of an ancestor of the refinement.
{
  Range G=numberOfComponentGrids();  // original number of grids
  
  Integer n = rcData->addRefinement(range, factor, level, k);

  updateReferences();
  return n;
}


//! Replace refinement levels "level0" and higher
/*!
  This function is used by the AMR Regrid class in order to efficiently replace a collection
  of refinement grids. This function avoids the overhead of calling addRefinement and deleteRefinement
  many times.

 \param level0,numberOfRefinementLevels0 : replace and/or add levels level0,..,numberOfRefinementLevels0-1
 \param gridInfo[bg][lev](0:ni-1,0:ng-1) : info defining a new refinement grid on base grid bg 
       and refinement level=level0+lev, lev=0,1,.... If we let 
             IntegerArray & info = gridInfo[bg][lev]
       then the number of new refinement grids is given by info.getLength(1).
       The first 6 entries in info define the range(0:1,0:2) of the refinement grid and the
       next three entries define the refinement factors along each axis,
          info(0,g) = range(0,0)
          info(1,g) = range(1,0)
          info(2,g) = range(0,1)
          info(3,g) = range(1,1)
          info(4,g) = range(0,2)
          info(5,g) = range(1,2)
          info(6,g) = factor(0)    
          info(7,g) = factor(1) 
          info(8,g) = factor(2) 

 */
Integer CompositeGrid::
replaceRefinementLevels(int level0, int numberOfRefinementLevels0, IntegerArray **gridInfo )
{
  int returnValue=rcData->replaceRefinementLevels(level0,numberOfRefinementLevels0,gridInfo );
  updateReferences();
  return returnValue;
}


//
// Delete all multigrid levels of refinement grid k.
//
void CompositeGrid::
deleteRefinement(const Integer& k)
{
   rcData->deleteRefinement(k);
   updateReferences();
}

//
// Delete all grids with refinement level greater than the given level.
//
void CompositeGrid::deleteRefinementLevels(const Integer level) {
    rcData->deleteRefinementLevels(level);
    updateReferences();
}
//
// Reference x[i] for refinementLevelNumber(i) <= level.
// Delete all other grids.
//
void CompositeGrid::referenceRefinementLevels(
  GenericGridCollection& x,
  const Integer          level) {
    rcData->referenceRefinementLevels(*((CompositeGrid&)x).rcData, level);
    updateReferences();
}
//
// Add a multigrid coarsening of grid k.
//
Integer CompositeGrid::addMultigridCoarsening(
  const IntegerArray& factor, // The coarsening factor w.r.t level-1
  const Integer&      level,  // The multigrid level number.
  const Integer       k) {    // The index of the corresponding grid
                              // at any finer multigrid level.
    Integer n = rcData->addMultigridCoarsening(factor, level, k);
    updateReferences();
    return n;
}
//
// Add multigrid coarsenings of grids in order to complete the multigrid levels.
//
void CompositeGrid::makeCompleteMultigridLevels() {
    rcData->makeCompleteMultigridLevels();
    updateReferences();
}
//
// Delete grid k, a multigrid coarsening, and all of its multigrid coarsenings.
//
void CompositeGrid::deleteMultigridCoarsening(const Integer& k) {
    rcData->deleteMultigridCoarsening(k);
    updateReferences();
}
//
// Delete all of the grids with multigrid level greater than the given level.
//
void CompositeGrid::deleteMultigridLevels(const Integer level) {
    rcData->deleteMultigridLevels(level);
    updateReferences();
}
//
// Set the number of grids.
//
void CompositeGrid::setNumberOfGrids(const Integer& numberOfGrids_) {
    rcData->setNumberOfGrids(numberOfGrids_);
    updateReferences();
}
//
// Set the number of dimensions.
//
void CompositeGrid::setNumberOfDimensions(
  const Integer& numberOfDimensions_) {
    rcData->setNumberOfGrids(numberOfDimensions_);
    updateReferences();
}
//
// Set the number of dimensions and grids.
//
void CompositeGrid::setNumberOfDimensionsAndGrids(
  const Integer& numberOfDimensions_,
  const Integer& numberOfGrids_) {
    rcData->setNumberOfDimensionsAndGrids(numberOfDimensions_, numberOfGrids_);
    updateReferences();
}

//\begin{>>CompositeGridInclude.tex}{\subsubsection{sizeOf}}
real CompositeGrid::
sizeOf(FILE *file /* = NULL */, bool returnSizeOfReference /* = false */ ) const
// ==========================================================================
// /Description: 
//   Return number of bytes allocated by this grid; optionally print detailed info to a file
//
// /file (input) : optinally supply a file to write detailed info to. Choose file=stdout to
// write to standard output.
// /returnSizeOfReference (input): if true only count the items that would not be referenced if this
//   CompositeGrid were referenced to another.
// /Return value: the number of bytes.
//\end{CompositeGridInclude.tex}
//==========================================================================
{
  real size=sizeof(*this);
  size+=GridCollection::sizeOf(file,returnSizeOfReference)-sizeof(GridCollection);
  
  size+=numberOfInterpolationPoints.elementCount()*sizeof(int);
  size+=numberOfImplicitInterpolationPoints.elementCount()*sizeof(int);
  size+=interpolationStartEndIndex.elementCount()*sizeof(int);
  size+=interpolationIsImplicit.elementCount()*sizeof(int);
  size+=interpolationWidth.elementCount()*sizeof(int);
  size+=interpolationPreference.elementCount()*sizeof(int);
  size+=mayInterpolate.elementCount()*sizeof(int);
  size+=mayCutHoles.elementCount()*sizeof(int);
  size+=sharedSidesMayCutHoles.elementCount()*sizeof(int);
  size+=multigridCoarseningRatio.elementCount()*sizeof(int);
  size+=multigridProlongationWidth.elementCount()*sizeof(int);
  size+=multigridRestrictionWidth.elementCount()*sizeof(int);

  size+=interpolationOverlap.elementCount()*sizeof(real);
  size+=maximumHoleCuttingDistance.elementCount()*sizeof(real);

  if( numberOfGrids()>0 )
  {
    for( int g=0; g<numberOfGrids(); g++ )
    {
      if( g<interpolationCoordinates.getLength() )
        size+=interpolationCoordinates[g].elementCount()*sizeof(real);
      if( g<interpoleeGrid.getLength() )
        size+=interpoleeGrid[g].elementCount()*sizeof(int);
      if( g<variableInterpolationWidth.getLength() )
        size+=variableInterpolationWidth[g].elementCount()*sizeof(int);
//      if( g<interpoleePoint.getLength() )
//	size+=interpoleePoint[g].elementCount()*sizeof(int);
      if( g<interpoleeLocation.getLength() )
	size+=interpoleeLocation[g].elementCount()*sizeof(int);
      if( g<interpolationPoint.getLength() )
	size+=interpolationPoint[g].elementCount()*sizeof(int);
//       if( g<interpolationCondition.getLength() )
//         size+=interpolationCondition[g].elementCount()*sizeof(real);
    }
  }

  // what about multigridLevel ?? what arrays are not shared?

  size+=inverseCoordinates.sizeOf();
  size+=inverseGrid.sizeOf();
  if( rcData!=NULL )
  {
    size+=rcData->hybridConnectivity.sizeOf();

    typedef TrivialArray<BoundaryAdjustment,Range>     BoundaryAdjustmentArray;
    typedef TrivialArray<BoundaryAdjustmentArray,Range>BoundaryAdjustmentArray2;
    BoundaryAdjustmentArray2 & boundaryAdjustment = rcData->boundaryAdjustment;
    if( boundaryAdjustment.getNumberOfElements() ) 
    {
      for( int k1=0; k1<numberOfComponentGrids(); k1++ )
	for( int k2=0; k2<numberOfComponentGrids(); k2++ )
	{
	  BoundaryAdjustmentArray& bA12 = boundaryAdjustment(k1,k2);
	  if( bA12.getNumberOfElements()>0 ) 
	  {
	    for( int axis=0; axis<numberOfDimensions(); axis++ )
	    for( int side=0; side<=1; side++ )
	    {
	      BoundaryAdjustment& bA = bA12(side,axis);
              size+=bA.sizeOf();
	    }
	  }
	}
    }
  }
  
  if( multigridLevel.getLength()>0 )
  {
    for( int l=0; l<numberOfMultigridLevels(); l++ )
    {
      size+=multigridLevel[l].sizeOf(file,true);  // count up data that is not referenced.
    }
  }
  

  return size;
}


//\begin{>>CompositeGridInclude.tex}{\subsubsection{setOverlapParameters}}
int CompositeGrid::
setOverlapParameters()
// ==========================================================================
// /Description: 
//    Assign default values to the overlap parameters such as interpolationWidth etc.
//
//\end{CompositeGridInclude.tex}
//==========================================================================
{
  Range all;
  int k1;
  for( k1=0; k1<numberOfGrids(); k1++) 
  {
    maximumHoleCuttingDistance(all,all,k1)=SQRT(.1*REAL_MAX); // this will be squared
  }
    

  int l;
  for( l=0; l<numberOfMultigridLevels(); l++ )
  {
    for( k1=0; k1<numberOfGrids(); k1++) 
    {
      for (int k2=0; k2<numberOfGrids(); k2++) 
      {
	interpolationIsImplicit(k1,k2,l) = true;

        int axis;
	for( axis=0; axis<numberOfDimensions(); axis++ )
	{
	  interpolationWidth(axis,k1,k2,l)              =3;
	  interpolationOverlap(axis,k1,k2,l)            = .5;
	  multigridCoarseningRatio(axis,k1,l)           = 2;
	  multigridProlongationWidth(axis,k1,l)         = 3;
	  multigridRestrictionWidth(axis,k1,l)          = 3;
	}
	// note: may cut holes does not have multigrid levels.
	mayCutHoles(k1,k2)=TRUE;   
	sharedSidesMayCutHoles(k1,k2)=FALSE;
      
	for( axis=numberOfDimensions(); axis<3; axis++ )
	{
	  interpolationWidth(axis,k1,k2,l)              =1;
	  interpolationOverlap(axis,k1,k2,l)            = -.5;
	  multigridCoarseningRatio(axis,k1,l)           = 1;
	  multigridProlongationWidth(axis,k1,l)         = 1;
	  multigridRestrictionWidth(axis,k1,l)          = 1;
	}
      
	interpolationPreference(k1,k2,l)           = k1;
	mayInterpolate(k1,k2,l)                    = true;
	  
      } // for k2
    }
  }
  return 0;
}

//\begin{>>CompositeGridInclude.tex}{\subsubsection{setOverlapParameters}}
int CompositeGrid::
setOverlapParameters(CompositeGrid & cg)
// ==========================================================================
// /Description: 
//    Assign values to the overlap parameters such as interpolationWidth etc.
//  based on the values in the grid cg.
//
//\end{CompositeGridInclude.tex}
//==========================================================================
{
  if( cg.numberOfGrids()!=numberOfGrids() )
    setOverlapParameters();
  
  int numberOfGridsToSet=min(cg.numberOfGrids(),numberOfGrids());

  Range all;
  int k1;
  for( k1=0; k1<numberOfGridsToSet; k1++) 
  {
    maximumHoleCuttingDistance(all,all,k1)=cg.maximumHoleCuttingDistance(all,all,k1);
  }
    

  int l;
  for( l=0; l<numberOfMultigridLevels(); l++ )
  {
    for( k1=0; k1<numberOfGridsToSet; k1++) 
    {
      for (int k2=0; k2<numberOfGridsToSet; k2++) 
      {
	interpolationIsImplicit(k1,k2,l) = cg.interpolationIsImplicit(k1,k2,l);

        int axis;
	for( axis=0; axis<numberOfDimensions(); axis++ )
	{
	  interpolationWidth(axis,k1,k2,l)              =cg.interpolationWidth(axis,k1,k2,l);
	  interpolationOverlap(axis,k1,k2,l)            =cg.interpolationOverlap(axis,k1,k2,l);
	  multigridCoarseningRatio(axis,k1,l)           =cg.multigridCoarseningRatio(axis,k1,l);
	  multigridProlongationWidth(axis,k1,l)         =cg.multigridProlongationWidth(axis,k1,l);
	  multigridRestrictionWidth(axis,k1,l)          =cg.multigridRestrictionWidth(axis,k1,l);
	}
	// note: may cut holes does not have multigrid levels.
	mayCutHoles(k1,k2)=cg.mayCutHoles(k1,k2);
	sharedSidesMayCutHoles(k1,k2)=cg.sharedSidesMayCutHoles(k1,k2);
      
	for( axis=numberOfDimensions(); axis<3; axis++ )
	{
	  interpolationWidth(axis,k1,k2,l)              =1;
	  interpolationOverlap(axis,k1,k2,l)            = -.5;
	  multigridCoarseningRatio(axis,k1,l)           = 1;
	  multigridProlongationWidth(axis,k1,l)         = 1;
	  multigridRestrictionWidth(axis,k1,l)          = 1;
	}
      
	interpolationPreference(k1,k2,l)           = cg.interpolationPreference(k1,k2,l);
	mayInterpolate(k1,k2,l)                    = cg.mayInterpolate(k1,k2,l);
	  
      } // for k2
    }
  }
  return 0;
}



// Initialize the CompositeGrid with the given number of dimensions and grids.
// These grids have their gridNumbers, baseGridNumbers and componentGridNumbers
// set to [0, ..., numberOfGrids_-1], and their refinementLevelNumbers and
// multigridLevelNumbers set to zero.
//
void CompositeGrid::initialize(
  const Integer& numberOfDimensions_,
  const Integer& numberOfGrids_) {
    GridCollection::initialize(numberOfDimensions_, numberOfGrids_);
    rcData->initialize(numberOfDimensions_, numberOfGrids_);
}
//
// Stream output operator.
//
ostream& operator<<(ostream& s, const CompositeGrid& g) {
    Integer i, k, k1, k2, l;
    s << (GridCollection&)g << endl
      << "  numberOfCompleteMultigridLevels() =  "
      <<  g.numberOfCompleteMultigridLevels() << endl
      << "  epsilon()                         =  "
      <<  g.epsilon() << endl
      << "  numberOfInterpolationPoints       = [";
    for (i=0; i<g.numberOfGrids(); i++)
      s << (i ? "," : "")
        << g.numberOfInterpolationPoints(i);
    s << "]" << endl
//      << "  numberOfInterpoleePoints          = [";
//    for (i=0; i<g.numberOfGrids(); i++)
//      s << (i ? "," : "")
//        << g.numberOfInterpoleePoints(i);
//    s << "]" << endl
      << "  interpolationIsAllExplicit()      =  "
      <<  (g.interpolationIsAllExplicit() ? 'T' : 'F') << endl
      << "  interpolationIsAllImplicit()      =  "
      <<  (g.interpolationIsAllImplicit() ? 'T' : 'F') << endl
      << "  interpolationIsImplicit           = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
            if (k2) s << ";";
            for (k1=0; k1<g.numberOfComponentGrids(); k1++)
              s << (k1 ? "," : "")
                << (g.interpolationIsImplicit(k1,k2,l) ? 'T' : 'F');
        } // end for
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
    s << "]" << endl
//       << "  backupInterpolationIsImplicit     = [";
//     for (l=0; l<g.numberOfMultigridLevels(); l++) {
//         if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
//         for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
//             if (k2) s << ";";
//             for (k1=0; k1<g.numberOfComponentGrids(); k1++)
//               s << (k1 ? "," : "")
//                 << (g.backupInterpolationIsImplicit(k1,k2,l) ? 'T' : 'F');
//         } // end for
//         if (g.numberOfMultigridLevels() > 1) s << ")";
//     } // end for
//    s << "]" << endl
      << "  interpolationWidth                = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
            if (k2) s << ";";
            for (k1=0; k1<g.numberOfComponentGrids(); k1++)
              s << (k1 ? "," : "")
                << g.interpolationWidth(0,k1,k2,l) << ":"
                << g.interpolationWidth(1,k1,k2,l) << ":"
                << g.interpolationWidth(2,k1,k2,l);
        } // end for
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
    s << "]" << endl
//       << "  backupInterpolationWidth          = [";
//     for (l=0; l<g.numberOfMultigridLevels(); l++) {
//         if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
//         for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
//             if (k2) s << ";";
//             for (k1=0; k1<g.numberOfComponentGrids(); k1++)
//               s << (k1 ? "," : "")
//                 << g.backupInterpolationWidth(0,k1,k2,l) << ":"
//                 << g.backupInterpolationWidth(1,k1,k2,l) << ":"
//                 << g.backupInterpolationWidth(2,k1,k2,l);
//         } // end for
//         if (g.numberOfMultigridLevels() > 1) s << ")";
//     } // end for
//    s << "]" << endl
      << "  interpolationOverlap              = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
            if (k2) s << ";";
            for (k1=0; k1<g.numberOfComponentGrids(); k1++)
              s << (k1 ? "," : "")
                << g.interpolationOverlap(0,k1,k2,l) << ":"
                << g.interpolationOverlap(1,k1,k2,l) << ":"
                << g.interpolationOverlap(2,k1,k2,l);
        } // end for
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
    s << "]" << endl
//       << "  backupInterpolationOverlap        = [";
//     for (l=0; l<g.numberOfMultigridLevels(); l++) {
//         if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
//         for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
//             if (k2) s << ";";
//             for (k1=0; k1<g.numberOfComponentGrids(); k1++)
//               s << (k1 ? "," : "")
//                 << g.backupInterpolationOverlap(0,k1,k2,l) << ":"
//                 << g.backupInterpolationOverlap(1,k1,k2,l) << ":"
//                 << g.backupInterpolationOverlap(2,k1,k2,l);
//         } // end for
//         if (g.numberOfMultigridLevels() > 1) s << ")";
//     } // end for
//     s << "]" << endl
//       << "  interpolationConditionLimit       = [";
//     for (l=0; l<g.numberOfMultigridLevels(); l++) {
//         if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
//         for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
//             if (k2) s << ";";
//             for (k1=0; k1<g.numberOfComponentGrids(); k1++)
//               s << (k1 ? "," : "")
//                 << g.interpolationConditionLimit(k1,k2,l);
//         } // end for
//         if (g.numberOfMultigridLevels() > 1) s << ")";
//     } // end for
//     s << "]" << endl
//       << "  backupInterpolationConditionLimit = [";
//     for (l=0; l<g.numberOfMultigridLevels(); l++) {
//         if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
//         for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
//             if (k2) s << ";";
//             for (k1=0; k1<g.numberOfComponentGrids(); k1++)
//               s << (k1 ? "," : "")
//                 << g.backupInterpolationConditionLimit(k1,k2,l);
//         } // end for
//         if (g.numberOfMultigridLevels() > 1) s << ")";
//     } // end for
//     s << "]" << endl
      << "  interpolationPreference           = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
            if (k2) s << ";";
            for (k1=0; k1<g.numberOfComponentGrids(); k1++)
              s << (k1 ? "," : "")
                << g.interpolationPreference(k1,k2,l);
        } // end for
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
    s << "]" << endl
      << "  mayInterpolate                    = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
            if (k2) s << ";";
            for (k1=0; k1<g.numberOfComponentGrids(); k1++)
              s << (k1 ? "," : "")
                << (g.mayInterpolate(k1,k2,l) ? 'T' : 'F');
        } // end for
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
    s << "]" << endl
//       << "  mayBackupInterpolate              = [";
//     for (l=0; l<g.numberOfMultigridLevels(); l++) {
//         if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
//         for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
//             if (k2) s << ";";
//             for (k1=0; k1<g.numberOfComponentGrids(); k1++)
//               s << (k1 ? "," : "")
//                 << (g.mayBackupInterpolate(k1,k2,l) ? 'T' : 'F');
//         } // end for
//         if (g.numberOfMultigridLevels() > 1) s << ")";
//     } // end for
//     s << "]" << endl
      << "  mayCutHoles                       = [";
    for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
        if (k2) s << ";";
        for (k1=0; k1<g.numberOfComponentGrids(); k1++)
          s << (k1 ? "," : "")
            << (g.mayCutHoles(k1,k2) ? 'T' : 'F');
    } // end for
    s << "]" << endl
      << "  multigridCoarseningRatio          = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k=0; k<g.numberOfComponentGrids(); k++)
          s << (k ? "," : "")
            << g.multigridCoarseningRatio(0,k,l) << ":"
            << g.multigridCoarseningRatio(1,k,l) << ":"
            << g.multigridCoarseningRatio(2,k,l);
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
    s << "]" << endl
      << "  multigridProlongationWidth        = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k=0; k<g.numberOfComponentGrids(); k++)
          s << (k ? "," : "")
            << g.multigridProlongationWidth(0,k,l) << ":"
            << g.multigridProlongationWidth(1,k,l) << ":"
            << g.multigridProlongationWidth(2,k,l);
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
    s << "]" << endl
      << "  multigridRestrictionWidth         = [";
    for (l=0; l<g.numberOfMultigridLevels(); l++) {
        if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
        for (k=0; k<g.numberOfComponentGrids(); k++)
          s << (k ? "," : "")
            << g.multigridRestrictionWidth(0,k,l) << ":"
            << g.multigridRestrictionWidth(1,k,l) << ":"
            << g.multigridRestrictionWidth(2,k,l);
        if (g.numberOfMultigridLevels() > 1) s << ")";
    } // end for
//    s << "]" << endl
//       << "  interpoleeGridRange               = [";
//     for (l=0; l<g.numberOfMultigridLevels(); l++) {
//         if (g.numberOfMultigridLevels() > 1) s << (l ? ",(" : "(");
//         for (k2=0; k2<g.numberOfComponentGrids(); k2++) {
//             if (k2) s << ";";
//             for (k1=0; k1<=g.numberOfComponentGrids(); k1++)
//               s << (k1 ? "," : "")
//                 << g.interpoleeGridRange(k1,k2,l);
//         } // end for
//         if (g.numberOfMultigridLevels() > 1) s << ")";
//     } // end for
   return s
      << "]";
}

//
// class CompositeGridData:
//
CompositeGridData::CompositeGridData(
  const Integer numberOfDimensions_,
  const Integer numberOfComponentGrids_):
  GridCollectionData(numberOfDimensions_, numberOfComponentGrids_) {
    className = "CompositeGridData";
    initialize(numberOfDimensions_, numberOfComponentGrids_);
}
CompositeGridData::CompositeGridData(
  const CompositeGridData& x,
  const CopyType           ct):
  GridCollectionData() {
    className = "CompositeGridData";
    initialize(numberOfDimensions, numberOfGrids);
    if (ct != NOCOPY) *this = x;
}
CompositeGridData::~CompositeGridData() { }
CompositeGridData& CompositeGridData::
operator=(const CompositeGridData& x) 
{
  GridCollectionData::operator=(x);

  numberOfCompleteMultigridLevels   = x.numberOfCompleteMultigridLevels;
  epsilon                           = x.epsilon;
  numberOfInterpolationPoints.redim(0);
  numberOfInterpolationPoints       = x.numberOfInterpolationPoints;

  numberOfImplicitInterpolationPoints.redim(0); 
  numberOfImplicitInterpolationPoints=x.numberOfImplicitInterpolationPoints;
  
  interpolationStartEndIndex.redim(0);
  interpolationStartEndIndex=x.interpolationStartEndIndex;
  
//  numberOfInterpoleePoints.redim(0);
//  numberOfInterpoleePoints          = x.numberOfInterpoleePoints;
  interpolationIsAllExplicit        = x.interpolationIsAllExplicit;
  interpolationIsAllImplicit        = x.interpolationIsAllImplicit;
  interpolationIsImplicit.redim(0);
  interpolationIsImplicit           = x.interpolationIsImplicit;
//  backupInterpolationIsImplicit.redim(0);
//  backupInterpolationIsImplicit     = x.backupInterpolationIsImplicit;
  interpolationWidth.redim(0);
  interpolationWidth                = x.interpolationWidth;
//  backupInterpolationWidth.redim(0);
//  backupInterpolationWidth          = x.backupInterpolationWidth;
  interpolationOverlap.redim(0);
  interpolationOverlap              = x.interpolationOverlap;
  maximumHoleCuttingDistance.redim(0);
  maximumHoleCuttingDistance=x.maximumHoleCuttingDistance;
//  backupInterpolationOverlap.redim(0);
//  backupInterpolationOverlap        = x.backupInterpolationOverlap;
  interpolationPreference.redim(0);
  interpolationPreference           = x.interpolationPreference;
  mayInterpolate.redim(0);
  mayInterpolate                    = x.mayInterpolate;
//  mayBackupInterpolate.redim(0);
//  mayBackupInterpolate              = x.mayBackupInterpolate;
  mayCutHoles.redim(0);
  mayCutHoles                       = x.mayCutHoles;
  sharedSidesMayCutHoles.redim(0);
  sharedSidesMayCutHoles            = x.sharedSidesMayCutHoles;
  multigridCoarseningRatio.redim(0);
  multigridCoarseningRatio          = x.multigridCoarseningRatio;
  multigridProlongationWidth.redim(0);
  multigridProlongationWidth        = x.multigridProlongationWidth;
  multigridRestrictionWidth.redim(0);
  multigridRestrictionWidth         = x.multigridRestrictionWidth;
//  interpoleeGridRange.redim(0);
//  interpoleeGridRange               = x.interpoleeGridRange;
//  interpolationConditionLimit.redim(0);
//  interpolationConditionLimit       = x.interpolationConditionLimit;
//  backupInterpolationConditionLimit.redim(0);
//  backupInterpolationConditionLimit = x.backupInterpolationConditionLimit;

  Integer upd = NOTHING, des = NOTHING;
#ifdef USE_STL
  if (x.interpolationCoordinates.size())
    upd |= THEinterpolationCoordinates;
  else des |= THEinterpolationCoordinates;
  if (x.interpoleeGrid          .size())
    upd |= THEinterpoleeGrid;
  else des |= THEinterpoleeGrid;
//   if (x.variableInterpolationWidth.size())
//     upd |= THEvariableInterpolationWidth;
//   else des |= THEvariableInterpolationWidth;
  if (x.interpoleeLocation      .size())
    upd |= THEinterpoleeLocation;
  else des |= THEinterpoleeLocation;
  if (x.interpolationPoint      .size())
    upd |= THEinterpolationPoint;
  else des |= THEinterpolationPoint;
//   if (x.interpolationCondition  .size())
//     upd |= THEinterpolationCondition;
//  else des |= THEinterpolationCondition;
#else
  if (x.interpolationCoordinates.getLength())
    upd |= THEinterpolationCoordinates;
  else des |= THEinterpolationCoordinates;
  if (x.interpoleeGrid          .getLength())
    upd |= THEinterpoleeGrid;
  else des |= THEinterpoleeGrid;
//   if (x.variableInterpolationWidth.getLength())
//     upd |= THEvariableInterpolationWidth;
//   else des |= THEvariableInterpolationWidth;
  if (x.interpoleeLocation      .getLength())
    upd |= THEinterpoleeLocation;
  else des |= THEinterpoleeLocation;
  if (x.interpolationPoint      .getLength())
    upd |= THEinterpolationPoint;
  else des |= THEinterpolationPoint;
//   if (x.interpolationCondition  .getLength())
//     upd |= THEinterpolationCondition;
//   else des |= THEinterpolationCondition;
#endif // USE_STL
  if ( // x.inverseCondition  .gridCollectionData == (GridCollectionData *)(&x) &&
      x.inverseCoordinates.gridCollectionData == (GridCollectionData *)(&x) &&
      x.inverseGrid       .gridCollectionData == (GridCollectionData *)(&x))
    upd |= THEinverseMap; else des |= THEinverseMap;
  if (upd &= ~des) update(upd, COMPUTEnothing); if (des) destroy(des);
//
//****************************************************************
//***** Assume that the x data have the expected dimensions. *****
//*****          This is probably a bad assumption.          *****
//****************************************************************
//
#ifdef USE_STL
  if (x.interpolationCoordinates.size())
    interpolationCoordinates = x.interpolationCoordinates;
  if (x.interpoleeGrid.size()) {
    interpoleeGrid         = x.interpoleeGrid;
//    interpoleePoint        = x.interpoleePoint;
    variableInterpolationWidth = x.variableInterpolationWidth;
  } // end if
  if (x.interpoleeLocation.size())
    interpoleeLocation       = x.interpoleeLocation;
  if (x.interpolationPoint.size())
    interpolationPoint       = x.interpolationPoint;
//   if (x.interpolationCondition.size())
//     interpolationCondition   = x.interpolationCondition;
#else
  if (x.interpolationCoordinates.getLength())
    interpolationCoordinates = x.interpolationCoordinates;
  if (x.interpoleeGrid.getLength()) {
    interpoleeGrid         = x.interpoleeGrid;
//    interpoleePoint        = x.interpoleePoint;
    variableInterpolationWidth=x.variableInterpolationWidth;
  } // end if
  if (x.interpoleeLocation.getLength())
    interpoleeLocation       = x.interpoleeLocation;
  if (x.interpolationPoint.getLength())
    interpolationPoint       = x.interpolationPoint;
//   if (x.interpolationCondition.getLength())
//     interpolationCondition   = x.interpolationCondition;
#endif // USE_STL
//   if ( x.inverseCondition.gridCollectionData == this) {
//     inverseCondition = x.inverseCondition;
//     inverseCondition.updateToMatchGrid(*this);
//   }
  if (x.inverseCoordinates.gridCollectionData == this) {
    inverseCoordinates = x.inverseCoordinates;
    if( !inverseCoordinates.isNull() ) // *wdh* 011126
      inverseCoordinates.updateToMatchGrid(*this);
  }
  if (x.inverseGrid.gridCollectionData == this) {
    inverseGrid = x.inverseGrid;
    if( !inverseGrid.isNull() ) // *wdh* 011126
      inverseGrid.updateToMatchGrid(*this);
  }
  computedGeometry |= upd & x.computedGeometry;

  hybridConnectivity = x.hybridConnectivity;

  return *this;
}
void CompositeGridData::reference(const CompositeGridData& x) {
    cerr << "CompositeGridData::reference(const CompositeGridData&) "
         << "was called!" << endl;
    GridCollectionData::reference(x);
}
void CompositeGridData::breakReference() {
    cerr << "CompositeGridData::breakReference() was called!" << endl;
    GridCollectionData::breakReference();
}
void CompositeGridData::consistencyCheck() const {
    GridCollectionData::              consistencyCheck();
    numberOfInterpolationPoints      .Test_Consistency();
    numberOfImplicitInterpolationPoints.Test_Consistency(); 
    interpolationStartEndIndex.Test_Consistency();          

//    numberOfInterpoleePoints         .Test_Consistency();
    interpolationIsImplicit          .Test_Consistency();
//    backupInterpolationIsImplicit    .Test_Consistency();
    interpolationWidth               .Test_Consistency();
//    backupInterpolationWidth         .Test_Consistency();
    interpolationOverlap             .Test_Consistency();
    maximumHoleCuttingDistance.Test_Consistency();
//    backupInterpolationOverlap       .Test_Consistency();
//    interpolationConditionLimit      .Test_Consistency();
//    backupInterpolationConditionLimit.Test_Consistency();
    interpolationPreference          .Test_Consistency();
    mayInterpolate                   .Test_Consistency();
//    mayBackupInterpolate             .Test_Consistency();
    mayCutHoles                      .Test_Consistency();
    sharedSidesMayCutHoles           .Test_Consistency();
    multigridCoarseningRatio         .Test_Consistency();
    multigridProlongationWidth       .Test_Consistency();
    multigridRestrictionWidth        .Test_Consistency();
//    interpoleeGridRange              .Test_Consistency();
    interpolationCoordinates         .consistencyCheck();
    interpoleeGrid                   .consistencyCheck();
    variableInterpolationWidth       .consistencyCheck();
//    interpoleePoint                  .consistencyCheck();
    interpoleeLocation               .consistencyCheck();
    interpolationPoint               .consistencyCheck();
//    interpolationCondition           .consistencyCheck();
    multigridLevel                   .consistencyCheck();
//    inverseCondition                 .consistencyCheck();
    inverseCoordinates               .consistencyCheck();
    inverseGrid                      .consistencyCheck();
}
Integer CompositeGridData::get(
  const GenericDataBase& db,
  const aString&         name) {
    Integer returnValue = 0;
    GenericDataBase& dir = *db.virtualConstructor();
    db.find(dir, name, getClassName());
    dir.setMode(GenericDataBase::streamInputMode);

    returnValue |= GridCollectionData::get(dir, "GridCollectionData");

    const Integer computedGeometry0 = computedGeometry;
    initialize(numberOfDimensions, numberOfGrids);

    returnValue |= dir.get(numberOfCompleteMultigridLevels,
                          "numberOfCompleteMultigridLevels");
    returnValue |= dir.get(epsilon,
                          "epsilon");
#if defined GNU || defined __PHOTON || defined __DECCXX
    {
    Integer foo;
    returnValue |= dir.get(foo,
                          "interpolationIsAllExplicit");
    interpolationIsAllExplicit = foo;
    returnValue |= dir.get(foo,
                          "interpolationIsAllImplicit");
    interpolationIsAllImplicit = foo;
    }
#else
    returnValue |= dir.get(interpolationIsAllExplicit,
                          "interpolationIsAllExplicit");
    returnValue |= dir.get(interpolationIsAllImplicit,
                          "interpolationIsAllImplicit");
#endif // defined GNU || defined __PHOTON || defined __DECCXX

    if (numberOfGrids > 0) {
        returnValue |= dir.get(numberOfInterpolationPoints,
                              "numberOfInterpolationPoints");
        returnValue |= dir.get(numberOfImplicitInterpolationPoints,
                              "numberOfImplicitInterpolationPoints");
        returnValue |= dir.get(interpolationStartEndIndex,
                              "interpolationStartEndIndex");
//        returnValue |= dir.get(numberOfInterpoleePoints,
//                              "numberOfInterpoleePoints");
    } // end if

    if (numberOfComponentGrids > 0) {
        returnValue |= dir.get(interpolationIsImplicit,
                              "interpolationIsImplicit");
//        returnValue |= dir.get(backupInterpolationIsImplicit,
//                              "backupInterpolationIsImplicit");
        returnValue |= dir.get(interpolationWidth,
                              "interpolationWidth");
//        returnValue |= dir.get(backupInterpolationWidth,
//                              "backupInterpolationWidth");
        returnValue |= dir.get(interpolationOverlap,
                              "interpolationOverlap");
        returnValue |= dir.get(maximumHoleCuttingDistance,
                              "maximumHoleCuttingDistance");
//        returnValue |= dir.get(backupInterpolationOverlap,
//                              "backupInterpolationOverlap");
//        returnValue |= dir.get(interpolationConditionLimit,
//                              "interpolationConditionLimit");
//        returnValue |= dir.get(backupInterpolationConditionLimit,
//                              "backupInterpolationConditionLimit");
        returnValue |= dir.get(interpolationPreference,
                              "interpolationPreference");
        returnValue |= dir.get(mayInterpolate,
                              "mayInterpolate");
//        returnValue |= dir.get(mayBackupInterpolate,
//                              "mayBackupInterpolate");
        returnValue |= dir.get(mayCutHoles,
                              "mayCutHoles");
        returnValue |= dir.get(sharedSidesMayCutHoles,
                              "sharedSidesMayCutHoles");
        returnValue |= dir.get(multigridCoarseningRatio,
                              "multigridCoarseningRatio");
        returnValue |= dir.get(multigridProlongationWidth,
                              "multigridProlongationWidth");
        returnValue |= dir.get(multigridRestrictionWidth,
                              "multigridRestrictionWidth");
//        returnValue |= dir.get(interpoleeGridRange,
//                              "interpoleeGridRange");
    } // end if

    CompositeGridData::update((GenericGridCollectionData&)*this,
      computedGeometry0 & (EVERYTHING & ~THElists), COMPUTEnothing);
    computedGeometry = computedGeometry0 & ~THElists;

    for (Integer i=0; i<numberOfGrids; i++)
      if (numberOfInterpolationPoints(i)) {
        char thing_i[32];
        if (computedGeometry & THEinterpolationCoordinates) {
            sprintf(thing_i,      "interpolationCoordinates[%d]", i);
            returnValue |=
              dir.get(interpolationCoordinates[i], thing_i);
        } // end if
        if (computedGeometry & THEinterpoleeGrid)
        {
            sprintf(thing_i,      "interpoleeGrid[%d]", i);
            returnValue |= dir.get(interpoleeGrid[i], thing_i);
            sprintf(thing_i,      "variableInterpolationWidth[%d]", i);
            int rt = dir.get(variableInterpolationWidth[i], thing_i);
            returnValue |= rt;
            if( rt!=0 )
	    {
              printf("Giving default values for variableInterpolationWidth : %i\n",max(interpolationWidth));
              variableInterpolationWidth[i].redim(numberOfInterpolationPoints(i));
              variableInterpolationWidth[i]=max(interpolationWidth);
	    }
        } // end if
        if (computedGeometry & THEinterpoleeLocation) {
            sprintf(thing_i,      "interpoleeLocation[%d]", i);
            returnValue |= dir.get(interpoleeLocation[i], thing_i);
        } // end if
        if (computedGeometry & THEinterpolationPoint) {
            sprintf(thing_i,      "interpolationPoint[%d]", i);
            returnValue |= dir.get(interpolationPoint[i], thing_i);
        } // end if
//         if (computedGeometry & THEinterpolationCondition) {
//             sprintf(thing_i,      "interpolationCondition[%d]", i);
//             returnValue |= dir.get(interpolationCondition[i], thing_i);
//         } // end if
    } // end if, end for

//     if (computedGeometry & THEinterpoleeGrid)
//      for (Integer i=0; i<numberOfGrids; i++)
//       if (numberOfInterpoleePoints(i))
//       {
//         char thing_i[32];
//         sprintf(thing_i, "interpoleePoint[%d]", i);             // *wdh* what is this ?
//         returnValue |= dir.get(interpoleePoint[i], thing_i);
//       } // end if, end for, end if

    if (computedGeometry & THEinverseMap) returnValue |=
//      inverseCondition  .get(dir, "inverseCondition")   |
      inverseCoordinates.get(dir, "inverseCoordinates") |
      inverseGrid       .get(dir, "inverseGrid");

    CompositeGridData::update((GenericGridCollectionData&)*this,
      computedGeometry0 & THElists, COMPUTEnothing);
    computedGeometry = computedGeometry0;

    // // //
    // kkc 5/23/01 added io of hybrid connectivity
    //
    int hybridConnectivitySaved=0;
    dir.get(hybridConnectivitySaved,"hybridConnectivitySaved");   // *wdh* do this way for streaming mode
    if ( hybridConnectivitySaved ) 
    {
      int unstructuredGridIndex=-1;                    
      dir.get(unstructuredGridIndex,"unstructuredGridIndex");
      intArray ugi;
      dir.get(ugi,"hybridUVertex2GridIndex");
      intArray bfacem;
      dir.get(bfacem,"hybridBoundaryFaceMapping");
	
      intArray *gi2uvptr = new intArray[numberOfGrids-1];
      intArray *gv2uvptr = new intArray[numberOfGrids-1];
      aString buff;
      for ( int g=0; g<numberOfGrids-1; g++ )
      {
	buff ="";
	sPrintF(buff,"hybridGridIndex2UVertex[%d]",g);
	dir.get(gi2uvptr[g],buff);
	buff ="";
	sPrintF(buff,"hybridGridVertex2UVertex[%d]",g);
	dir.get(gv2uvptr[g],buff);
      }

      hybridConnectivity.setCompositeGridHybridConnectivity(unstructuredGridIndex,
							    gi2uvptr,
							    ugi,
							    gv2uvptr,
							    bfacem);
    }

    // // //

    delete &dir;
    return returnValue;
}
Integer CompositeGridData::put(
  GenericDataBase& db,
  const aString&   name) const {
    Integer returnValue = 0;
    GenericDataBase& dir = *db.virtualConstructor();
    db.create(dir, name, getClassName());
    dir.setMode(GenericDataBase::streamOutputMode);

    returnValue |= GridCollectionData::put(dir, "GridCollectionData");

    returnValue |= dir.put(numberOfCompleteMultigridLevels,
                          "numberOfCompleteMultigridLevels");
    returnValue |= dir.put(epsilon, "epsilon");
    returnValue |= dir.put(interpolationIsAllExplicit,
                          "interpolationIsAllExplicit");
    returnValue |= dir.put(interpolationIsAllImplicit,
                          "interpolationIsAllImplicit");

    if (numberOfGrids > 0) {
        returnValue |= dir.put(numberOfInterpolationPoints,
                              "numberOfInterpolationPoints");
        returnValue |= dir.put(numberOfImplicitInterpolationPoints,
                              "numberOfImplicitInterpolationPoints");
        returnValue |= dir.put(interpolationStartEndIndex,
                              "interpolationStartEndIndex");
//        returnValue |= dir.put(numberOfInterpoleePoints,
//                              "numberOfInterpoleePoints");
    } // end if

    if (numberOfComponentGrids > 0) {
        returnValue |= dir.put(interpolationIsImplicit,
                              "interpolationIsImplicit");
//        returnValue |= dir.put(backupInterpolationIsImplicit,
//                              "backupInterpolationIsImplicit");
        returnValue |= dir.put(interpolationWidth,
                              "interpolationWidth");
//        returnValue |= dir.put(backupInterpolationWidth,
//                              "backupInterpolationWidth");
        returnValue |= dir.put(interpolationOverlap,
                              "interpolationOverlap");
        returnValue |= dir.put(maximumHoleCuttingDistance,
                              "maximumHoleCuttingDistance");
//        returnValue |= dir.put(backupInterpolationOverlap,
//                              "backupInterpolationOverlap");
//        returnValue |= dir.put(interpolationConditionLimit,
//                              "interpolationConditionLimit");
//        returnValue |= dir.put(backupInterpolationConditionLimit,
//                              "backupInterpolationConditionLimit");
        returnValue |= dir.put(interpolationPreference,
                              "interpolationPreference");
        returnValue |= dir.put(mayInterpolate,
                              "mayInterpolate");
//        returnValue |= dir.put(mayBackupInterpolate,
//                              "mayBackupInterpolate");
        returnValue |= dir.put(mayCutHoles,
                              "mayCutHoles");
        returnValue |= dir.put(sharedSidesMayCutHoles,
                              "sharedSidesMayCutHoles");
        returnValue |= dir.put(multigridCoarseningRatio,
                              "multigridCoarseningRatio");
        returnValue |= dir.put(multigridProlongationWidth,
                              "multigridProlongationWidth");
        returnValue |= dir.put(multigridRestrictionWidth,
                              "multigridRestrictionWidth");
//        returnValue |= dir.put(interpoleeGridRange,
//                              "interpoleeGridRange");
    } // end if

    for (Integer i=0; i<numberOfGrids; i++)
      if (numberOfInterpolationPoints(i)) {
        char thing_i[32];
        if (computedGeometry &  THEinterpolationCoordinates) {
            sprintf(thing_i,      "interpolationCoordinates[%d]", i);
            returnValue |= dir.put(interpolationCoordinates[i], thing_i);
        } // end if
        if (computedGeometry &  THEinterpoleeGrid) {
            sprintf(thing_i,      "interpoleeGrid[%d]", i);
            returnValue |= dir.put(interpoleeGrid[i], thing_i);
            sprintf(thing_i,      "variableInterpolationWidth[%d]", i);
            returnValue |= dir.put(variableInterpolationWidth[i], thing_i);
        } // end if
        if (computedGeometry &  THEinterpoleeLocation) {
            sprintf(thing_i,      "interpoleeLocation[%d]", i);
            returnValue |= dir.put(interpoleeLocation[i], thing_i);
        } // end if
        if (computedGeometry &  THEinterpolationPoint) {
            sprintf(thing_i,      "interpolationPoint[%d]", i);
            returnValue |= dir.put(interpolationPoint[i], thing_i);
        } // end if
//         if (computedGeometry &  THEinterpolationCondition) {
//             sprintf(thing_i,      "interpolationCondition[%d]", i);
//             returnValue |= dir.put(interpolationCondition[i], thing_i);
//         } // end if
    } // end if, end for

//     if (computedGeometry &  THEinterpoleeGrid)
//      for (Integer i=0; i<numberOfGrids; i++)
//       if (numberOfInterpoleePoints(i)) {
//         char thing_i[32];
//         sprintf(thing_i, "interpoleePoint[%d]", i);
//         returnValue |= dir.put(interpoleePoint[i], thing_i);
//     } // end if, end for, end if

    if (computedGeometry & THEinverseMap) returnValue |=
//      inverseCondition  .put(dir, "inverseCondition")   |
      inverseCoordinates.put(dir, "inverseCoordinates") |
      inverseGrid       .put(dir, "inverseGrid");

    // // //
    // kkc 5/23/01 added output of hybrid connectivity
    //
    int hybridConnectivitySaved=hybridConnectivity.getUnstructuredGridIndex()>-1;
    returnValue |= dir.put(hybridConnectivitySaved,"hybridConnectivitySaved");  // *wdh* do this way for streaming mode
    if ( hybridConnectivity.getUnstructuredGridIndex()>-1 ) 
    {
      returnValue |= dir.put(hybridConnectivity.getUnstructuredGridIndex(),"unstructuredGridIndex");

      returnValue |= dir.put(hybridConnectivity.getBoundaryFaceMapping(),"hybridBoundaryFaceMapping");
      returnValue |= dir.put(hybridConnectivity.getUVertex2GridIndex(),"hybridUVertex2GridIndex");
	
      aString buff="";
      for ( int g=0; g<numberOfGrids-1; g++ )
      {
	buff ="";
	const intArray & gi2uv = hybridConnectivity.getGridIndex2UVertex(g);
	sPrintF(buff,"hybridGridIndex2UVertex[%d]",g);
	dir.put(gi2uv,buff);
	buff ="";
	const intArray & gv2uv = hybridConnectivity.getGridVertex2UVertex(g);
	sPrintF(buff,"hybridGridVertex2UVertex[%d]",g);
	dir.put(gv2uv,buff);
      }
    }

    // // //

    delete &dir;
    return returnValue;
}
Integer CompositeGridData::update(
  GenericGridCollectionData& x,
  const Integer              what,
  const Integer              how) {
    Integer upd = GridCollectionData::update(x, what & ~THEmultigridLevel, how);
    CompositeGridData& y = (CompositeGridData&)x;
    Integer computeNeeded =
      how & COMPUTEgeometry         ? what :
      how & COMPUTEgeometryAsNeeded ? what & ~computedGeometry :
                                      NOTHING;
//
//  Compute interpolationIsAllExplicit and interpolationIsAllImplicit from
//  values of interpolationIsImplicit, backupInterpolationIsImplicit,
//  mayInterpolate and mayBackupInterpolate.
//
    interpolationIsAllExplicit =
    interpolationIsAllImplicit = LogicalTrue;
    for (Integer l=0; l<numberOfMultigridLevels; l++)
      for (Integer k1=0; k1<numberOfComponentGrids; k1++)
        for (Integer k2=0; k2<numberOfComponentGrids; k2++)
          if (k1 != k2) {
            if ((mayInterpolate(k1,k2,l) &&
                 interpolationIsImplicit(k1,k2,l)) 
		//  || ( mayBackupInterpolate(k1,k2,l) && backupInterpolationIsImplicit(k1,k2,l))
               )
              interpolationIsAllExplicit = LogicalFalse;
            if ((mayInterpolate(k1,k2,l) &&
                 !interpolationIsImplicit(k1,k2,l)) 
		//  || ( mayBackupInterpolate(k1,k2,l) && !backupInterpolationIsImplicit(k1,k2,l))
               )
              interpolationIsAllImplicit = LogicalFalse;
    } // end if, end for, end for, end for

#ifdef USE_STL
    if (what & THEinterpolationCoordinates) {
        Integer i = numberOfGrids - interpolationCoordinates.size();
        if (i < 0) interpolationCoordinates.erase(
          interpolationCoordinates.begin() + numberOfGrids,
          interpolationCoordinates.end());
        else for (Integer j=0; j<i; j++) interpolationCoordinates.push_back();
    } // end if
    if (what & THEinterpoleeGrid) 
    {
      Integer i = numberOfGrids - interpoleeGrid.size();
      if (i < 0) 
        interpoleeGrid.erase(interpoleeGrid.begin() + numberOfGrids, interpoleeGrid.end());
      else 
        for (Integer j=0; j<i; j++) interpoleeGrid.push_back();
      i = numberOfGrids - variableInterpolationWidth.size();
      if (i < 0)
         variableInterpolationWidth.erase(variableInterpolationWidth.begin() + numberOfGrids,
                                          variableInterpolationWidth.end());
      else
         for (Integer j=0; j<i; j++) variableInterpolationWidth.push_back();
//       i = numberOfGrids - interpoleePoint.size();
//       if (i < 0)
//          interpoleePoint.erase(	interpoleePoint.begin() + numberOfGrids,interpoleePoint.end());
//       else
//          for (Integer j=0; j<i; j++) interpoleePoint.push_back();
    } // end if
    if (what & THEinterpoleeLocation) {
        Integer i = numberOfGrids - interpoleeLocation.size();
        if (i < 0) interpoleeLocation.erase(
          interpoleeLocation.begin() + numberOfGrids,
          interpoleeLocation.end());
        else for (Integer j=0; j<i; j++) interpoleeLocation.push_back();
    } // end if
    if (what & THEinterpolationPoint) {
        Integer i = numberOfGrids - interpolationPoint.size();
        if (i < 0) interpolationPoint.erase(
          interpolationPoint.begin() + numberOfGrids,
          interpolationPoint.end());
        else for (Integer j=0; j<i; j++) interpolationPoint.push_back();
    } // end if
//     if (what & THEinterpolationCondition) {
//         Integer i = numberOfGrids - interpolationCondition.size();
//         if (i < 0) interpolationCondition.erase(
//           interpolationCondition.begin() + numberOfGrids,
//           interpolationCondition.end());
//         else for (Integer j=0; j<i; j++) interpolationCondition.push_back();
//     } // end if
#else
    if (what & THEinterpolationCoordinates) {
        while (interpolationCoordinates.getLength() < numberOfGrids)
          interpolationCoordinates     .addElement();
        while (interpolationCoordinates.getLength() > numberOfGrids)
          interpolationCoordinates     .deleteElement();
    } // end if
    if (what & THEinterpoleeGrid) 
    {
        while (interpoleeGrid          .getLength() < numberOfGrids)
          interpoleeGrid               .addElement();
        while (interpoleeGrid          .getLength() > numberOfGrids)
          interpoleeGrid               .deleteElement();

        while (variableInterpolationWidth.getLength() < numberOfGrids)
          variableInterpolationWidth     .addElement();
        while (variableInterpolationWidth.getLength() > numberOfGrids)
          variableInterpolationWidth     .deleteElement();

//         while (interpoleePoint         .getLength() < numberOfGrids)
//           interpoleePoint              .addElement();
//         while (interpoleePoint         .getLength() > numberOfGrids)
//           interpoleePoint              .deleteElement();
    } // end if
    if (what & THEinterpoleeLocation) {
        while (interpoleeLocation      .getLength() < numberOfGrids)
          interpoleeLocation           .addElement();
        while (interpoleeLocation      .getLength() > numberOfGrids)
          interpoleeLocation           .deleteElement();
    } // end if
    if (what & THEinterpolationPoint) {
        while (interpolationPoint      .getLength() < numberOfGrids)
          interpolationPoint           .addElement();
        while (interpolationPoint      .getLength() > numberOfGrids)
          interpolationPoint           .deleteElement();
    } // end if
//     if (what & THEinterpolationCondition) {
//         while (interpolationCondition  .getLength() < numberOfGrids)
//           interpolationCondition       .addElement();
//         while (interpolationCondition  .getLength() > numberOfGrids)
//           interpolationCondition       .deleteElement();
//     } // end if
#endif // USE_STL

    for (Integer i=0; i<numberOfGrids; i++) {
        if (numberOfInterpolationPoints(i)) {
            if (what & THEinterpolationCoordinates) {
                if (&y != this && i <
#ifdef USE_STL
                  y.interpolationCoordinates.size() &&
#else
                  y.interpolationCoordinates.getLength() &&
#endif // USE_STL
                  y.interpolationCoordinates[i].elementCount() ==
                    numberOfInterpolationPoints(i) * numberOfDimensions &&
                  y.interpolationCoordinates[i].getBase(0) == 0 &&
                  y.interpolationCoordinates[i].getBound(0) ==
                    numberOfInterpolationPoints(i) - 1 &&
                  y.interpolationCoordinates[i].getBase(1) == 0 &&
                  y.interpolationCoordinates[i].getBound(1) ==
                    numberOfDimensions - 1)
                  interpolationCoordinates[i].
                    reference(y.interpolationCoordinates[i]);
                if (interpolationCoordinates[i].elementCount() !=
                    numberOfInterpolationPoints(i) * numberOfDimensions ||
                    interpolationCoordinates[i].getBase(0) != 0 ||
                    interpolationCoordinates[i].getBound(0) !=
                    numberOfInterpolationPoints(i) - 1 ||
                    interpolationCoordinates[i].getBase(1) != 0 ||
                    interpolationCoordinates[i].getBound(1) !=
                    numberOfDimensions - 1) {
                    interpolationCoordinates[i].redim(
                      numberOfInterpolationPoints(i), numberOfDimensions);
                    interpolationCoordinates[i] = (Real)0.;
                    if (how & COMPUTEgeometryAsNeeded)
                      computeNeeded  |=  THEinterpolationCoordinates;
                    computedGeometry &= ~THEinterpolationCoordinates;
                    upd              |=  THEinterpolationCoordinates;
                } // end if
            } // end if
            if (what & THEinterpoleeGrid)
            {
                if (&y != this && i <
#ifdef USE_STL
                  y.interpoleeGrid.size() &&
#else
                  y.interpoleeGrid.getLength() &&
#endif // USE_STL
                  y.interpoleeGrid[i].elementCount() == numberOfInterpolationPoints(i) &&
                  y.interpoleeGrid[i].getBase(0) == 0 && 
		    y.interpoleeGrid[i].getBound(0) == numberOfInterpolationPoints(i) - 1)
		{
                  interpoleeGrid[i].reference(y.interpoleeGrid[i]);
                  variableInterpolationWidth[i].reference(y.variableInterpolationWidth[i]);
		}
                if (interpoleeGrid[i].elementCount() !=numberOfInterpolationPoints(i) ||
                    interpoleeGrid[i].getBase(0) != 0 || 
                    interpoleeGrid[i].getBound(0) != numberOfInterpolationPoints(i) - 1)
                {
                    interpoleeGrid[i].redim(numberOfInterpolationPoints(i));
                    interpoleeGrid[i] = 0;
                    variableInterpolationWidth[i].redim(numberOfInterpolationPoints(i));
                    variableInterpolationWidth[i] = 0;
                    if (how & COMPUTEgeometryAsNeeded)
                      computeNeeded  |=  THEinterpoleeGrid;
                    computedGeometry &= ~THEinterpoleeGrid;
                    upd              |=  THEinterpoleeGrid;
                } // end if
            } // end if
            if (what & THEinterpoleeLocation) {
                if (&y != this && i <
#ifdef USE_STL
                  y.interpoleeLocation.size() &&
#else
                  y.interpoleeLocation.getLength() &&
#endif // USE_STL
                  y.interpoleeLocation[i].elementCount() !=
                    numberOfInterpolationPoints(i) * numberOfDimensions &&
                  y.interpoleeLocation[i].getBase(0) == 0 &&
                  y.interpoleeLocation[i].getBound(0) ==
                    numberOfInterpolationPoints(i) - 1 &&
                  y.interpoleeLocation[i].getBase(1) == 0 &&
                  y.interpoleeLocation[i].getBound(1) ==
                    numberOfDimensions - 1)
                  interpoleeLocation[i].
                    reference(y.interpoleeLocation[i]);
                if (interpoleeLocation[i].elementCount() !=
                    numberOfInterpolationPoints(i) * numberOfDimensions ||
                    interpoleeLocation[i].getBase(0) != 0 ||
                    interpoleeLocation[i].getBound(0) !=
                    numberOfInterpolationPoints(i) - 1 ||
                    interpoleeLocation[i].getBase(1) != 0 ||
                    interpoleeLocation[i].getBound(1) !=
                    numberOfDimensions - 1) {
                    interpoleeLocation[i].redim(
                      numberOfInterpolationPoints(i), numberOfDimensions);
                    interpoleeLocation[i] = 0;
                    if (how & COMPUTEgeometryAsNeeded)
                      computeNeeded  |=  THEinterpoleeLocation;
                    computedGeometry &= ~THEinterpoleeLocation;
                    upd              |=  THEinterpoleeLocation;
                } // end if
            } // end if
            if (what & THEinterpolationPoint) {
                if (&y != this && i <
#ifdef USE_STL
                  y.interpolationPoint.size() &&
#else
                  y.interpolationPoint.getLength() &&
#endif // USE_STL
                  y.interpolationPoint[i].elementCount() ==
                    numberOfInterpolationPoints(i) * numberOfDimensions &&
                  y.interpolationPoint[i].getBase(0) == 0 &&
                  y.interpolationPoint[i].getBound(0) ==
                    numberOfInterpolationPoints(i) - 1 &&
                  y.interpolationPoint[i].getBase(1) == 0 &&
                  y.interpolationPoint[i].getBound(1) ==
                    numberOfDimensions - 1)
                  interpolationPoint[i].
                    reference(y.interpolationPoint[i]);
                if (interpolationPoint[i].elementCount() !=
                    numberOfInterpolationPoints(i) * numberOfDimensions ||
                    interpolationPoint[i].getBase(0) != 0 ||
                    interpolationPoint[i].getBound(0) !=
                    numberOfInterpolationPoints(i) - 1 ||
                    interpolationPoint[i].getBase(1) != 0 ||
                    interpolationPoint[i].getBound(1) !=
                    numberOfDimensions - 1) {
                    interpolationPoint[i].redim(
                      numberOfInterpolationPoints(i), numberOfDimensions);
                    interpolationPoint[i] = 0;
                    if (how & COMPUTEgeometryAsNeeded)
                      computeNeeded  |=  THEinterpolationPoint;
                    computedGeometry &= ~THEinterpolationPoint;
                    upd              |=  THEinterpolationPoint;
                } // end if
            } // end if
//             if (what & THEinterpolationCondition) {
//                 if (&y != this && i <
// #ifdef USE_STL
//                   y.interpolationCondition.size() &&
// #else
//                   y.interpolationCondition.getLength() &&
// #endif // USE_STL
//                   y.interpolationCondition[i].elementCount() ==
//                     numberOfInterpolationPoints(i) &&
//                   y.interpolationCondition[i].getBase(0) == 0 &&
//                   y.interpolationCondition[i].getBound(0) ==
//                     numberOfInterpolationPoints(i) - 1)
//                   interpolationCondition[i].
//                     reference(y.interpolationCondition[i]);
//                 if (interpolationCondition[i].elementCount() !=
//                     numberOfInterpolationPoints(i) ||
//                     interpolationCondition[i].getBase(0) != 0 ||
//                     interpolationCondition[i].getBound(0) !=
//                     numberOfInterpolationPoints(i) - 1) {
//                     interpolationCondition[i].redim(
//                       numberOfInterpolationPoints(i));
//                     interpolationCondition[i] = (Real)0.;
//                     if (how & COMPUTEgeometryAsNeeded)
//                       computeNeeded  |=  THEinterpolationCondition;
//                     computedGeometry &= ~THEinterpolationCondition;
//                     upd              |=  THEinterpolationCondition;
//                 } // end if
//             } // end if
        } else { // (numberOfInterpolationPoints == 0)
            if (what & THEinterpolationCoordinates)
              interpolationCoordinates[i].redim(0);
            if (what & THEinterpoleeGrid)
	    {
              interpoleeGrid[i]            .redim(0);
	      variableInterpolationWidth[i].redim(0);
	    }
            if (what & THEinterpoleeLocation)
              interpoleeLocation[i]      .redim(0);
            if (what & THEinterpolationPoint)
              interpolationPoint[i]      .redim(0);
//             if (what & THEinterpolationCondition)
//               interpolationCondition[i]  .redim(0);
        } // end if
//         if (what & THEinterpoleeGrid) {
//             if (numberOfInterpoleePoints(i)) {
//                 if (&y != this && i <
// #ifdef USE_STL
//                   y.interpoleePoint.size() &&
// #else
//                   y.interpoleePoint.getLength() &&
// #endif // USE_STL
//                   y.interpoleePoint[i].elementCount() ==
//                     numberOfInterpoleePoints(i) &&
//                   y.interpoleePoint[i].getBase(0) == 0 &&
//                   y.interpoleePoint[i].getBound(0) ==
//                     numberOfInterpoleePoints(i) - 1)
//                   interpoleePoint[i].
//                     reference(y.interpoleePoint[i]);
//                 if (interpoleePoint[i].elementCount() !=
//                     numberOfInterpoleePoints(i) ||
//                     interpoleePoint[i].getBase(0) != 0 ||
//                     interpoleePoint[i].getBound(0) !=
//                     numberOfInterpoleePoints(i) - 1) {
//                     interpoleePoint[i].redim(
//                       numberOfInterpoleePoints(i));
//                     interpoleePoint[i] = 0;
//                     if (how & COMPUTEgeometryAsNeeded)
//                       computeNeeded  |=  THEinterpoleeGrid;
//                     computedGeometry &= ~THEinterpoleeGrid;
//                     upd              |=  THEinterpoleeGrid;
//                 } // end if
//             } else { // (numberOfInterpoleePoints == 0)
//                 interpoleePoint[i].redim(0);
//             } // end if
//         } // end if
    } // end for

    const Range all, nd1 = numberOfDimensions;
    if (what & THEinverseMap) {
        if (&y != this) {
//            inverseCondition  .reference(y.inverseCondition);
            inverseCoordinates.reference(y.inverseCoordinates);
            inverseGrid       .reference(y.inverseGrid);
            if (y.computedGeometry &   THEinverseMap) {
                computedGeometry   |=  THEinverseMap;
            } else if (how         &   COMPUTEgeometryAsNeeded) {
                computedGeometry   &= ~THEinverseMap;
                computeNeeded      |=  THEinverseMap;
            } // end if
        } // end if
        if((  // inverseCondition  .updateToMatchGrid(*this, all, all, all) |
	  // *wdh*     inverseCoordinates.updateToMatchGrid(*this, nd1, all, all, all) |
	  inverseCoordinates.updateToMatchGrid(*this, all, all, all, nd1) |
	  inverseGrid       .updateToMatchGrid(*this, all, all, all)) &
           RealCompositeGridFunction::updateResized) {
	  // inverseCondition       =   (Real)0.;
	  inverseCoordinates     =   (Real)0.;
	  inverseGrid            =   -1;
	  if (how                &   COMPUTEgeometryAsNeeded)
	    computeNeeded        |=  THEinverseMap;
	  computedGeometry       &= ~THEinverseMap;
	  upd                    |=  THEinverseMap;
	} // end if

// *wdh* 000202
//         boundaryAdjustment.redim(numberOfGrids,numberOfBaseGrids);
//         for (Integer k1=0; k1<numberOfGrids; k1++) {
//             MappedGrid& g1 = grid[k1];
//             for (Integer k2=0; k2<numberOfBaseGrids; k2++) {
//                 TrivialArray<BoundaryAdjustment,Range>&
//                   bA12 = boundaryAdjustment(k1,k2);
//                 assert(baseGridNumber(k2) == k2);
//                 if (k2 == baseGridNumber(k1)) bA12.redim(0); else {
//                     MappedGrid& g2 = grid[k2];
//                     bA12.redim(2,numberOfDimensions);
//                     Logical noSharedBoundaries = LogicalTrue;
//                     for (Integer kd1=0; kd1<numberOfDimensions; kd1++)
//                       for (Integer ks1=0; ks1<2; ks1++) {
//                         BoundaryAdjustment& bA = bA12(ks1,kd1);
//                         Logical needAdjustment = LogicalFalse;
//                         for (Integer kd2=0; kd2<numberOfDimensions; kd2++)
//                           for (Integer ks2=0; ks2<2; ks2++)
//                             if (g1.boundaryCondition(ks1,kd1) > 0 &&
//                                 g2.boundaryCondition(ks2,kd2) > 0 &&
//                                 g1.sharedBoundaryFlag(ks1,kd1) &&
//                                 g2.sharedBoundaryFlag(ks2,kd2) ==
//                                 g1.sharedBoundaryFlag(ks1,kd1))
//                               needAdjustment = LogicalTrue;
//                         if (needAdjustment) {
//                             noSharedBoundaries = LogicalFalse;
//                             const Integer side = ks1 ?
//                               RealMappedGridFunction::endingGridIndex :
//                               RealMappedGridFunction::startingGridIndex;
//                             const Range d0 = numberOfDimensions,
//                               d1 = kd1==0 ? Range(side,side) : Range(),
//                               d2 = kd1==1 ? Range(side,side) : Range(),
//                               d3 = kd1==2 ? Range(side,side) : Range();
//                             if (( bA.boundaryAdjustment
//                                     .updateToMatchGrid(g1, d1, d2, d3, d0) |
//                                   bA.acrossGrid
//                                     .updateToMatchGrid(g1, d1, d2, d3, d0) |
//                                   bA.oppositeBoundary
//                                     .updateToMatchGrid(g1, d1, d2, d3, d0) )
//                               & RealMappedGridFunction::updateResized) {
//                                 bA.computedGeometry &= ~THEinverseMap;
//                                 if (how          &  COMPUTEgeometryAsNeeded)
//                                   computeNeeded  |=  THEinverseMap;
//                                 computedGeometry &= ~THEinverseMap;
//                                 upd              |=  THEinverseMap;
//                             } // end if
//                         } else {
//                             bA.computedGeometry &= ~THEinverseMap;
//                             bA.boundaryAdjustment.destroy();
//                             bA.acrossGrid        .destroy();
//                             bA.oppositeBoundary  .destroy();
//                         } // end if
//                     } // end for, end for
//                     if (noSharedBoundaries) bA12.redim(0);
//                 } // end if
//             } // end for
//         } // end for
    } // end if what

    if (what &                THEmultigridLevel)
      upd |= updateCollection(THEmultigridLevel | (what & ~THElists),
        numberOfMultigridLevels, multigridLevel,
        GridCollectionData::multigridLevel,
        GenericGridCollectionData::multigridLevel, multigridLevelNumber);

    upd |= computeGeometry(computeNeeded, how);

    return upd;
}
void CompositeGridData::destroy(const Integer what) {
#ifdef USE_STL
    if (what & THEinterpolationCoordinates) interpolationCoordinates.erase(
      interpolationCoordinates.begin(),     interpolationCoordinates.end());
    if (what & THEinterpoleeGrid)
    {
      interpoleeGrid.erase(interpoleeGrid.begin(),interpoleeGrid.end());
      variableInterpolationWidth.erase(variableInterpolationWidth.begin(),variableInterpolationWidth.end());
    }
//     if (what & THEinterpoleeGrid)           interpoleePoint.erase(
//       interpoleePoint.begin(),              interpoleePoint.end());
    if (what & THEinterpoleeLocation)       interpoleeLocation.erase(
      interpoleeLocation.begin(),           interpoleeLocation.end());
    if (what & THEinterpolationPoint)       interpolationPoint.erase(
      interpolationPoint.begin(),           interpolationPoint.end());
//     if (what & THEinterpolationCondition)   interpolationCondition.erase(
//       interpolationCondition.begin(),       interpolationCondition.end());
    if (what & THEmultigridLevel)           multigridLevel.erase(
      multigridLevel.begin(),               multigridLevel.end());
#else
    if (what & THEinterpolationCoordinates)
      while (interpolationCoordinates.getLength())
        interpolationCoordinates.deleteElement();
    if (what & THEinterpoleeGrid)
    {
      while (interpoleeGrid.getLength())
        interpoleeGrid          .deleteElement();
      while (variableInterpolationWidth.getLength())
        variableInterpolationWidth.deleteElement();
    }
//     if (what & THEinterpoleeGrid)
//       while (interpoleePoint.getLength())
//         interpoleePoint         .deleteElement();
    if (what & THEinterpoleeLocation)
      while (interpoleeLocation.getLength())
        interpoleeLocation      .deleteElement();
    if (what & THEinterpolationPoint)
      while (interpolationPoint.getLength())
        interpolationPoint      .deleteElement();
//     if (what & THEinterpolationCondition)
//       while (interpolationCondition.getLength())
//         interpolationCondition  .deleteElement();
    if (what & THEmultigridLevel)
      multigridLevel.reference(ListOfCompositeGrid());


#endif // USE_STL
    if (what & THEinverseMap) {
//         inverseCondition  .destroy();
         inverseCoordinates.destroy();
         inverseGrid       .destroy();
         boundaryAdjustment.redim(0);
    } // end if

    GridCollectionData::destroy(what);
}


//! Replace refinement level "level0" and higher
/*!

 \param level0,numberOfRefinementLevels0 : replace and/or add levels level0,..,numberOfRefinementLevels0-1
 \param gridInfo[bg][l](0:ni-1,0:ng-1) : info defining a new refinement grid on base grid bg 
       and refinement level=level0+l
 */
Integer CompositeGridData::
replaceRefinementLevels(int level0, int numberOfRefinementLevels0, IntegerArray **gridInfo )
{
  int returnValue=GridCollectionData::replaceRefinementLevels(level0,numberOfRefinementLevels0,gridInfo );
  if( returnValue!=0 ) return returnValue;
  
  // redimension arrays in the CompositeGrid.
  setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids);

  // assign values to the composite grid arrays.


  for( Integer k10=0; k10<numberOfGrids; k10++)
  {
    if( refinementLevelNumber(k10) >= level0 ) 
    {
      // this must be a new grid

      Integer k1 = componentGridNumber(k10), k3 = baseGridNumber(k10);
      assert(k3 == componentGridNumber(k3));

      const Integer l = multigridLevelNumber(k1);
      numberOfInterpolationPoints(k1) = 0;

      for (Integer k20=0; k20<numberOfGrids; k20++)
      {
        // *wdh* if( multigridLevelNumber(k20)  == l &&
        // 	    (refinementLevelNumber(k20) == refinementLevelNumber(n) ||
        // 	     refinementLevelNumber(k20) == refinementLevelNumber(n) - 1)) 
	if( multigridLevelNumber(k20)  == l )
	{
	  MappedGrid& g_20 = grid[k20];
	  Integer k2 = componentGridNumber(k20), k4 = baseGridNumber(k20);
	  assert(k4 == componentGridNumber(k4));
	  if (k20 == k10 || k3 != k4) 
	  {
	    // Interpolation from self or from an unrelated grid.
	    interpolationIsImplicit(k1,k2,l) = interpolationIsImplicit(k3,k4,l);
	    for (Integer kd=0; kd<3; kd++) 
	    {
	      interpolationWidth(kd,k1,k2,l) = interpolationWidth(kd,k3,k4,l);
	      interpolationOverlap(kd,k1,k2,l) = interpolationOverlap(kd,k3,k4,l);
	    } // end for

	  } 
	  else //  if (k3 == k4) 
	  {
	    //   Interpolation from a related grid.
	    interpolationIsImplicit(k1,k2,l) = LogicalFalse;
	    Integer kd;
	    for (kd=0; kd<numberOfDimensions; kd++) 
	    {
	      interpolationWidth(kd,k1,k2,l) = 1;
	      if (g_20.mapping().getGridDimensions(kd) == 1)
	      {
		// Interpolation on a surface grid.
		interpolationOverlap(kd,k1,k2,l)       = (Real)-.5;
	      } 
	      else if (refinementLevelNumber(k20) == refinementLevelNumber(k10))
	      {
		// Interpolation from a grid at the same refinement level.
		interpolationOverlap(kd,k1,k2,l) = amax1(epsilon,
							 (Real).5 * interpolationWidth(kd,k1,k2,l) - (Real)1. +
							 (Real).5 * (g_20.discretizationWidth(kd) - 1));
	      } 
	      else 
	      {
		//  Implicit interpolation from a parent (coarser) grid.
		// *wdh* 000821 interpolationIsImplicit(k1,k2,l) = LogicalTrue;
		interpolationIsImplicit(k1,k2,l) = false;
		interpolationWidth(kd,k1,k2,l) = interpolationWidth(kd,k3,k4,l);
		interpolationOverlap(kd,k1,k2,l) = amax1(epsilon,
							 (Real).5 * interpolationWidth(kd,k1,k2,l) - (Real)1.);
	      } // end if

	      if (refinementLevelNumber(k20) == refinementLevelNumber(k10) - 1) 
	      {
		//  Coarsen for multigrid in the same way as for any parent.
		multigridCoarseningRatio(kd,k1,l) = multigridCoarseningRatio(kd,k3,l);
		multigridProlongationWidth(kd,k1,l) = multigridProlongationWidth(kd,k3,l);
		multigridRestrictionWidth(kd,k1,l) = multigridRestrictionWidth(kd,k3,l);
	      } // end if
	    } // end for
	    for (kd=numberOfDimensions; kd<3; kd++) 
	    {
	      interpolationWidth(kd,k1,k2,l) = 1;
	      interpolationWidth(kd,k1,k2,l) = interpolationWidth(kd,k3,k4,l);
	      interpolationOverlap(kd,k1,k2,l) = interpolationOverlap(kd,k3,k4,l);
	    } // end for
	  } // end if
	  mayInterpolate(k1,k2,l) = k1 == k2 ? LogicalFalse : k3 == k4 ? LogicalTrue  :   mayInterpolate(k3,k4,l);
	  if (refinementLevelNumber(k20) == refinementLevelNumber(k10)) 
	  {
	    mayInterpolate(k2,k1,l)= k1 == k2 ? LogicalFalse : k3 == k4 ? LogicalTrue  : mayInterpolate(k4,k3,l);
	  } // end if
	} // end if
      } // end for k20

      //  Initially disallow interpolation to or from the new grid.
      //
      interpolationPreference(k1,Range(0,numberOfComponentGrids-1),l) = -1;
      interpolationPreference(Range(0,numberOfComponentGrids-1),k1,l) = -1;

    } // end if
  }

  // wdh: whay are these here?
//   Integer k1 = componentGridNumber(n), k3 = baseGridNumber(n);
//   assert(k3 == componentGridNumber(k3));
//   for (Integer kd=0; kd<3; kd++)
//   {
//     multigridCoarseningRatio(kd,k1,l) = multigridCoarseningRatio(kd,k3,l);
//     multigridProlongationWidth(kd,k1,l) = multigridProlongationWidth(kd,k3,l);
//     multigridRestrictionWidth(kd,k1,l) = multigridRestrictionWidth(kd,k3,l);
//   } // end for
// //

  if( computedGeometry & GridCollection::THEmultigridLevel )
    printf("***** CompositeGrid: replaceRefinementLevels START THEmultigridLevel!\n");
  return 0;
}


Integer CompositeGridData::
addRefinement(
  const IntegerArray& range,
  const IntegerArray& factor,
  const Integer&      level,
  const Integer       k) 
// ======================================================================================================
//   /Description:
//     Add a refinement.
//
// /Return value: grid number of the new grid added.
// ======================================================================================================
{
  Integer n = GridCollectionData::addRefinement(range, factor, level,  k);
  setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids);
  const Integer l = multigridLevelNumber(n);
  numberOfInterpolationPoints(n) = 0;

  for( Integer k10=0; k10<numberOfGrids; k10++)
  {
    if( multigridLevelNumber(k10)  == l ) // *wdh* && (refinementLevelNumber(k10) == refinementLevelNumber(n))) 
    {
      Integer k1 = componentGridNumber(k10), k3 = baseGridNumber(k10);
      assert(k3 == componentGridNumber(k3));
      for (Integer k20=0; k20<numberOfGrids; k20++)
      {
        // *wdh* if( multigridLevelNumber(k20)  == l &&
        // 	    (refinementLevelNumber(k20) == refinementLevelNumber(n) ||
        // 	     refinementLevelNumber(k20) == refinementLevelNumber(n) - 1)) 
	if( multigridLevelNumber(k20)  == l )
	{
	  MappedGrid& g_20 = grid[k20];
	  Integer k2 = componentGridNumber(k20), k4 = baseGridNumber(k20);
	  assert(k4 == componentGridNumber(k4));
	  if (k20 == k10 || k3 != k4) 
	  {
	    // Interpolation from self or from an unrelated grid.
	    interpolationIsImplicit(k1,k2,l) = interpolationIsImplicit(k3,k4,l);
	    for (Integer kd=0; kd<3; kd++) 
	    {
	      interpolationWidth(kd,k1,k2,l) = interpolationWidth(kd,k3,k4,l);
	      interpolationOverlap(kd,k1,k2,l) = interpolationOverlap(kd,k3,k4,l);
	    } // end for

	  } 
	  else //  if (k3 == k4) 
	  {
	    //   Interpolation from a related grid.
	    interpolationIsImplicit(k1,k2,l) = LogicalFalse;
	    Integer kd;
	    for (kd=0; kd<numberOfDimensions; kd++) 
	    {
	      interpolationWidth(kd,k1,k2,l) = 1;
	      if (g_20.mapping().getGridDimensions(kd) == 1)
	      {
		// Interpolation on a surface grid.
		interpolationOverlap(kd,k1,k2,l)       = (Real)-.5;
	      } 
	      else if (refinementLevelNumber(k20) == refinementLevelNumber(k10))
	      {
		// Interpolation from a grid at the same refinement level.
		interpolationOverlap(kd,k1,k2,l) = amax1(epsilon,
							 (Real).5 * interpolationWidth(kd,k1,k2,l) - (Real)1. +
							 (Real).5 * (g_20.discretizationWidth(kd) - 1));
	      } 
	      else 
	      {
		//  Implicit interpolation from a parent (coarser) grid.
		// *wdh* 000821 interpolationIsImplicit(k1,k2,l) = LogicalTrue;
		interpolationIsImplicit(k1,k2,l) = false;
		interpolationWidth(kd,k1,k2,l) = interpolationWidth(kd,k3,k4,l);
		interpolationOverlap(kd,k1,k2,l) = amax1(epsilon,
							 (Real).5 * interpolationWidth(kd,k1,k2,l) - (Real)1.);
	      } // end if

	      if (refinementLevelNumber(k20) == refinementLevelNumber(k10) - 1) 
	      {
		//  Coarsen for multigrid in the same way as for any parent.
		multigridCoarseningRatio(kd,k1,l) = multigridCoarseningRatio(kd,k3,l);
		multigridProlongationWidth(kd,k1,l) = multigridProlongationWidth(kd,k3,l);
		multigridRestrictionWidth(kd,k1,l) = multigridRestrictionWidth(kd,k3,l);
	      } // end if
	    } // end for
	    for (kd=numberOfDimensions; kd<3; kd++) 
	    {
	      interpolationWidth(kd,k1,k2,l) = 1;
	      interpolationWidth(kd,k1,k2,l) = interpolationWidth(kd,k3,k4,l);
	      interpolationOverlap(kd,k1,k2,l) = interpolationOverlap(kd,k3,k4,l);
	    } // end for
	  } // end if
	  mayInterpolate(k1,k2,l) = k1 == k2 ? LogicalFalse : k3 == k4 ? LogicalTrue  :   mayInterpolate(k3,k4,l);
	  if (refinementLevelNumber(k20) == refinementLevelNumber(n)) 
	  {
	    mayInterpolate(k2,k1,l)= k1 == k2 ? LogicalFalse : k3 == k4 ? LogicalTrue  : mayInterpolate(k4,k3,l);
	  } // end if
	} // end if
      }
    } // end if
  }

  Integer k1 = componentGridNumber(n), k3 = baseGridNumber(n);
  assert(k3 == componentGridNumber(k3));
  for (Integer kd=0; kd<3; kd++)
  {
    multigridCoarseningRatio(kd,k1,l) = multigridCoarseningRatio(kd,k3,l);
    multigridProlongationWidth(kd,k1,l) = multigridProlongationWidth(kd,k3,l);
    multigridRestrictionWidth(kd,k1,l) = multigridRestrictionWidth(kd,k3,l);
  } // end for
//
//  Initially disallow interpolation to or from the new grid.
//
  interpolationPreference(k1,Range(0,numberOfComponentGrids-1),l) = -1;
  interpolationPreference(Range(0,numberOfComponentGrids-1),k1,l) = -1;
  return n;
}

void CompositeGridData::deleteRefinement(const Integer& k) {
    if (k < 0 || k >= numberOfGrids) {
        cout << "CompositeGridData::deleteRefinement(k = "
             << k << "):  Grid " << k << " does not exist." << endl;
        assert(k >= 0); assert(k < numberOfGrids);
    } else if (refinementLevelNumber(k) == 0) {
        cout << "CompositeGridData::deleteRefinement(k = "
             << k << "):  Grid k = " << k << " is not a refinement." << endl;
        assert(refinementLevelNumber(k) != 0);
    } // end if
    CompositeGridData::deleteMultigridCoarsening(k);
}
void CompositeGridData::deleteRefinementLevels(const Integer level) {
    Integer i = numberOfGrids, j = i - 1;
    while (i--) if (refinementLevelNumber(i) > level && i < j--) {
        Range r1(i, j), r2 = r1 + 1;
        numberOfInterpolationPoints(r1) = numberOfInterpolationPoints(r2);
//        numberOfInterpoleePoints(r1)    = numberOfInterpoleePoints(r2);
    } // end if, end while
    GridCollectionData::deleteRefinementLevels(level);
    setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids);
}
void CompositeGridData::referenceRefinementLevels(
  GenericGridCollectionData& x,
  const Integer              level) {
    GridCollectionData::referenceRefinementLevels(x, level);
    setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids);
    CompositeGridData& y = (CompositeGridData&)x;
    for (Integer i=0, j=0; i<y.numberOfGrids; i++)
      if (y.refinementLevelNumber(i) <= level) {
        numberOfInterpolationPoints(j) = y.numberOfInterpolationPoints(i);
//        numberOfInterpoleePoints(j)    = y.numberOfInterpoleePoints(i);
        j++;
    } // end if, end for
}
Integer CompositeGridData::
addMultigridCoarsening( const IntegerArray& factor,
			const Integer&      level,
			const Integer       k)
{
  Integer n = GridCollectionData::addMultigridCoarsening(factor, level,  k);
  setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids);
  MappedGrid& g_n = grid[n];
  Integer k1 = componentGridNumber(n), k2, kd;

  assert( k1>=0 && k1<numberOfGrids);  // *wdh* 981118
    
  for (k2=0; k2<numberOfComponentGrids; k2++) 
  {
    interpolationIsImplicit(k1,k2,level) = interpolationIsImplicit(k1,k2,level-1);
    mayInterpolate(k1,k2,level) =  mayInterpolate(k1,k2,level-1);
    for (kd=0; kd<3; kd++) 
    {
      interpolationWidth(kd,k1,k2,level) =interpolationWidth(kd,k1,k2,level-1);
      interpolationOverlap(kd,k1,k2,level) = 	interpolationOverlap(kd,k1,k2,level-1);
    } // end for
  } // end for
  //
  //  Find the first component grid at the same multigrid level
  //  and refinement level.  If it is a different grid, we use its
  //  stencil widths.  Otherwise if level > 1, we use the stencil
  //  widths of the same component grid at the next-finer multigrid
  //  level.  Otherwise we use default stencil widths of two.
  //
  for (k2=0; k2<numberOfGrids; k2++) 
    if (k2 == n) 
      continue;
  else if (multigridLevelNumber(k2) == multigridLevelNumber(n) &&
	   baseGridNumber(k2)       == baseGridNumber(n)) 
    break;
  if (k2 != numberOfGrids) k2 = componentGridNumber(k2);

  for (kd=0; kd<numberOfDimensions; kd++)
  {
    multigridCoarseningRatio(kd,k1,level)   = factor(kd);
    multigridProlongationWidth(kd,k1,level) =
      k2 != numberOfGrids ? multigridProlongationWidth(kd,k2,level) :
      level > 1 ? multigridProlongationWidth(kd,k1,level-1) : 2;
    assert(multigridProlongationWidth(kd,k1,level) > 0);
    multigridRestrictionWidth(kd,k1,level)  =
      k2 != numberOfGrids ? multigridRestrictionWidth(kd,k2,level) :
      level > 1 ? multigridRestrictionWidth(kd,k1,level-1) : 2;
    assert(multigridRestrictionWidth(kd,k1,level) > 0);
    if (g_n.isCellCentered(kd)) 
    {
      //          multigridCoarseningRatio and multigridRestrictionWidth
      //          must be both odd or both even.
      if ((multigridRestrictionWidth(kd,k1,level)-multigridCoarseningRatio(kd,k1,level)) % 2)
	multigridRestrictionWidth(kd,k1,level)++;
    } 
    else 
    {
      //          multigridRestrictionWidth must be odd.
      if (multigridRestrictionWidth(kd,k1,level) % 2 == 0)
	multigridRestrictionWidth(kd,k1,level)++;
      //          multigridProlongationWidth must be even.
      if (multigridProlongationWidth(kd,k1,level) % 2)
	multigridProlongationWidth(kd,k1,level)++;
    } // end if
  } // end for
  for (kd=numberOfDimensions; kd<3; kd++) 
  {
    multigridCoarseningRatio(kd,k1,level)   = 1;
    multigridProlongationWidth(kd,k1,level) = 1;
    multigridRestrictionWidth(kd,k1,level)  = 1;
  } // end for
  //
  //  Initially disallow interpolation to or from the new grid.
  //
  interpolationPreference(k1,Range(0,numberOfComponentGrids-1),level) = -1;
  interpolationPreference(Range(0,numberOfComponentGrids-1),k1,level) = -1;
  return n;
}

void CompositeGridData::makeCompleteMultigridLevels() {
//
//  Find the coarsest multigrid level of each component grid.
//
    IntegerArray coarseLevel(numberOfComponentGrids); coarseLevel = -1;
    IntegerArray coarseGrid(numberOfComponentGrids);  coarseGrid  = -1;
    Integer g;
    for (g=0; g<numberOfGrids; g++) {
        const Integer k = componentGridNumber(g);
        if (coarseLevel(k) < multigridLevelNumber(g))
          { coarseLevel(k) = multigridLevelNumber(g); coarseGrid(k)  = g; }
    } // end for
//
//  Add multigrid coarsenings to the component grids that need more of them.
//  Use the same coarsening ratio as that of the corresponding base grid.
//  Note that if there is initially only one multigrid level of a base grid,
//  then all of the new multigrid levels will have coarsening factors of one.
//
    for (g=0; g<numberOfGrids; g++) {
        const Integer k = componentGridNumber(g), l = coarseLevel(k) + 1;
        if (l < numberOfCompleteMultigridLevels &&
          multigridLevelNumber(g) == l - 1) {
            const Integer b = baseGridNumber(g);
            assert(b == componentGridNumber(b));
            coarseLevel(k) = l; coarseGrid(k) = addMultigridCoarsening(
              multigridCoarseningRatio(Range(0,2),b,l), l, coarseGrid(k));
            assert(coarseGrid(k) > g);
        } // end if
    } // end for
}

void CompositeGridData::
deleteMultigridCoarsening(const Integer& k) 
{
  // printf("** CompositeGridData::deleteMultigridCoarsening k=%i\n",k);
  
  if (k < 0 || k >= numberOfGrids)
  {
    cout << "CompositeGridData::deleteMultigridCoarsening(k = "
	 << k << "):  Grid " << k << " does not exist." << endl;
    assert(k >= 0); assert(k < numberOfGrids);
  } 
  else if (multigridLevelNumber(k) == 0 && refinementLevelNumber(k) == 0)
  {
    cout << "CompositeGridData::deleteMultigridCoarsening(k = "
	 << k << "):  Grid k = " << k << " is not a multigrid coarsening."
	 << endl;
    assert(multigridLevelNumber(k) != 0 || refinementLevelNumber(k) != 0);
  } // end if
  Integer i, lastGrid=numberOfGrids-1;
  Range allGrids,allLevels;
  for( i=numberOfGrids-1; i>=0; i-- )
  {
    if (componentGridNumber(i) == componentGridNumber(k) && multigridLevelNumber(i) >= multigridLevelNumber(k) )
    {
      Range r1(i, lastGrid-1), r2=r1+1;
      numberOfInterpolationPoints(r1) = numberOfInterpolationPoints(r2);

      // *wdh* now update arrays
      interpolationIsImplicit(r1,allGrids,allLevels) = interpolationIsImplicit(r2,allGrids,allLevels);
      interpolationIsImplicit(allGrids,r1,allLevels) = interpolationIsImplicit(allGrids,r1,allLevels);
      for( int kd=0; kd<3; kd++ ) 
      {
        interpolationWidth(kd,r1,allGrids,allLevels)=interpolationWidth(kd,r2,allGrids,allLevels);
        interpolationWidth(kd,allGrids,r1,allLevels)=interpolationWidth(kd,allGrids,r2,allLevels);

        interpolationOverlap(kd,r1,allGrids,allLevels)=interpolationOverlap(kd,r2,allGrids,allLevels);
        interpolationOverlap(kd,allGrids,r1,allLevels)=interpolationOverlap(kd,allGrids,r2,allLevels);
      }
      
      lastGrid--;
    } // end if
  }
  
  GridCollectionData::deleteMultigridCoarsening(k);
  setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids);

}

void CompositeGridData::deleteMultigridLevels(const Integer level) {
    if (level < 0) {
        cout << "CompositeGridData::deleteMultigridLevel(level = "
             << level << "):  Multigrid level " << level << " does not exist."
             << endl;
        assert(level >= 0);
    } else if (level < numberOfMultigridLevels-1) {
        Integer i = numberOfGrids, j = i - 1;
        while (i--) if (multigridLevelNumber(i) > level && i < j--) {
            Range r1(i, j), r2 = r1 + 1;
            numberOfInterpolationPoints(r1) = numberOfInterpolationPoints(r2);
//            numberOfInterpoleePoints(r1)    = numberOfInterpoleePoints(r2);
        } // end if, end while
    } // end if
    GridCollectionData::deleteMultigridLevels(level);
    setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids);
}
void CompositeGridData::setNumberOfGrids(const Integer& numberOfGrids_)
  { setNumberOfDimensionsAndGrids(numberOfDimensions, numberOfGrids_); }
void CompositeGridData::setNumberOfDimensions(
  const Integer& numberOfDimensions_)
  { setNumberOfDimensionsAndGrids(numberOfDimensions_, numberOfGrids); }

void CompositeGridData::
setNumberOfDimensionsAndGrids(
  const Integer& numberOfDimensions_,
  const Integer& numberOfGrids_) 
{
  GridCollectionData::setNumberOfDimensionsAndGrids(numberOfDimensions_, numberOfGrids_);

  const Integer n = numberOfGrids_-numberOfInterpolationPoints.elementCount();
  if (n) 
  {
    numberOfInterpolationPoints.resize(numberOfGrids);
    numberOfImplicitInterpolationPoints.resize(numberOfGrids);
    interpolationStartEndIndex.resize(4,numberOfGrids,numberOfGrids);  
    if (n > 0) 
    {
      const Range newGrids(numberOfGrids_ - n, numberOfGrids_ - 1);
      numberOfInterpolationPoints(newGrids) = 0;
      Range all;
      for( int i=0; i<4; i++ )
      {
        interpolationStartEndIndex(i,all     ,newGrids)=-1;
        interpolationStartEndIndex(i,newGrids,all     )=-1;
      }
      
    } // end if
    computedGeometry &= ~(
      THEinterpolationCoordinates | THEinterpoleeGrid     |
      THEinterpoleeLocation       | THEinterpolationPoint |
      THEinterpolationCondition   | THEinverseMap         );
  } // end if

  const int numGrids=numberOfComponentGrids;
  // n1 : newNumberOfGrids - oldNumberOfGrids
  // n2 : newNumberOfMultigridLevels - oldNumber
  const Integer n1 = numGrids - (interpolationIsImplicit.getBound(0)-interpolationIsImplicit.getBase(0)+1),
    n2 = numberOfMultigridLevels -  (interpolationIsImplicit.getBound(2)-interpolationIsImplicit.getBase(2)+1);
  if (n1 || n2)
  {
    if (numGrids && numberOfMultigridLevels) 
    {
      interpolationIsImplicit          .resize(numGrids,numGrids,numberOfMultigridLevels);
      interpolationWidth               .resize(3, numGrids,numGrids,numberOfMultigridLevels);
      interpolationOverlap             .resize(3, numGrids,numGrids,numberOfMultigridLevels);
      maximumHoleCuttingDistance.resize(2,3,numGrids);
      interpolationPreference          .resize(numGrids,numGrids,numberOfMultigridLevels);
      mayInterpolate                   .resize(numGrids,numGrids,numberOfMultigridLevels);
      mayCutHoles                      .resize(numGrids,numGrids);
      sharedSidesMayCutHoles           .resize(numGrids,numGrids);
      multigridCoarseningRatio         .resize(3, numGrids,numberOfMultigridLevels);
      multigridProlongationWidth       .resize(3, numGrids,numberOfMultigridLevels);
      multigridRestrictionWidth        .resize(3, numGrids,numberOfMultigridLevels);
    } 
    else 
    {
      interpolationIsImplicit          .redim(0);
      interpolationWidth               .redim(0);
      interpolationOverlap             .redim(0);
      maximumHoleCuttingDistance.redim(0);
      interpolationPreference          .redim(0);
      mayInterpolate                   .redim(0);
      mayCutHoles                      .redim(0);
      sharedSidesMayCutHoles           .redim(0);
      multigridCoarseningRatio         .redim(0);
      multigridProlongationWidth       .redim(0);
      multigridRestrictionWidth        .redim(0);
    } // end if

    if (n1 > 0 && numGrids > 0) 
    {
      const Range three = 3, newComponentGrids(numGrids - n1, numGrids - 1), allComponentGrids = numGrids;
      mayCutHoles(newComponentGrids,allComponentGrids)= LogicalFalse;
      mayCutHoles(allComponentGrids,newComponentGrids)= LogicalFalse;
      sharedSidesMayCutHoles(newComponentGrids,allComponentGrids)= LogicalFalse;
      sharedSidesMayCutHoles(allComponentGrids,newComponentGrids)= LogicalFalse;
      maximumHoleCuttingDistance(nullRange,nullRange,newComponentGrids)=sqrt(.1*REAL_MAX);

      if (numberOfMultigridLevels > n2) 
      {
	const Range oldMultigridLevels = numberOfMultigridLevels - n2;
	interpolationIsImplicit(newComponentGrids,allComponentGrids,oldMultigridLevels) = LogicalFalse;
	interpolationIsImplicit(allComponentGrids,newComponentGrids,oldMultigridLevels) = LogicalFalse;

	interpolationWidth(three,newComponentGrids,allComponentGrids,oldMultigridLevels) = 0;
	interpolationWidth(three,allComponentGrids,newComponentGrids,oldMultigridLevels) = 0;

	interpolationOverlap(three,newComponentGrids,allComponentGrids,oldMultigridLevels) = 0.;
	interpolationOverlap(three,allComponentGrids,newComponentGrids,oldMultigridLevels) = 0.;
	interpolationPreference(newComponentGrids,allComponentGrids,oldMultigridLevels) = 0;
	interpolationPreference(allComponentGrids,newComponentGrids,oldMultigridLevels) = 0;
	mayInterpolate(allComponentGrids,newComponentGrids,oldMultigridLevels) = LogicalFalse;
	mayInterpolate(newComponentGrids,allComponentGrids,oldMultigridLevels) = LogicalFalse;
	multigridCoarseningRatio(three,newComponentGrids,oldMultigridLevels)=0;
	multigridProlongationWidth(three,newComponentGrids,oldMultigridLevels)=0;
	multigridRestrictionWidth(three,newComponentGrids,oldMultigridLevels)=0;
      } // end if
    } // end if
    if (n2 > 0 && numGrids > 0) 
    {
      const Range three = 3, newMultigridLevels(numberOfMultigridLevels - n2, numberOfMultigridLevels - 1),
	allComponentGrids = numGrids;
      interpolationIsImplicit(allComponentGrids,allComponentGrids,newMultigridLevels) = LogicalFalse;
      interpolationWidth(three,allComponentGrids,allComponentGrids,newMultigridLevels) = 0;
      interpolationOverlap(three,allComponentGrids,allComponentGrids,newMultigridLevels) = 0.;
      interpolationPreference(allComponentGrids,allComponentGrids,newMultigridLevels) = 0;
      mayInterpolate(allComponentGrids,allComponentGrids,newMultigridLevels) = LogicalFalse;
      multigridCoarseningRatio(three,allComponentGrids,newMultigridLevels)=0;
      multigridProlongationWidth(three,allComponentGrids,newMultigridLevels)=0;
      multigridRestrictionWidth(three,allComponentGrids,newMultigridLevels)=0;
    } // end if
  } // end if

}

void CompositeGridData::
initialize(const Integer& numberOfDimensions_,
	   const Integer& numberOfGrids_) 
{
  numberOfCompleteMultigridLevels=0; // *wdh* 000825
  CompositeGridData::setNumberOfDimensionsAndGrids
    (numberOfDimensions_, numberOfGrids_);
  destroy(~NOTHING & ~GridCollectionData::EVERYTHING);
//
//  Compute a default value for epsilon.
//
  epsilon = Mapping::epsilon();
}

void CompositeGridData::
getInterpolationStencil(const Integer&      k10,
			const Integer&      k20,
			const RealArray&    r,
			const IntegerArray& interpolationStencil,
			const intArray& useBackupRules) 
{
  MappedGrid& g = grid[k20];
  const Real a = -(Real)2. * epsilon, b = (Real)1. - a;
  const Integer base = r.getBase(0), bound = r.getBound(0),
    k1 = componentGridNumber(k10), k2 = componentGridNumber(k20),
    l  = multigridLevelNumber(k10);

#ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
#define g_boundaryCondition(i,j)     g.boundaryCondition ((i),(j))
#define g_gridSpacing(i)             g.gridSpacing       ((i))
#define g_indexRange(i,j)            g.indexRange        ((i),(j))
#define g_extendedIndexRange(i,j)    g.extendedIndexRange((i),(j))
#define g_isCellCentered(i)          g.isCellCentered    ((i))
#define g_isPeriodic(i)              g.isPeriodic        ((i))
#define r_(i,j)                      r                   ((i),(j))
#define interpolationStencil_(i,j,k) interpolationStencil((i),(j),(k))
#define useBackupRules_(i)           useBackupRules      ((i))
#define iw0_(i)                      iw0                 ((i),k1,k2,l)
#else
#define g_boundaryCondition(i,j)     g_boundaryCondition_ [(i) + 2 * (j)]
#define g_gridSpacing(i)             g_gridSpacing_       [(i)]
#define g_indexRange(i,j)            g_indexRange_        [(i) + 2 * (j)]
#define g_extendedIndexRange(i,j)    g_extendedIndexRange_[(i) + 2 * (j)]
#define g_isCellCentered(i)          g_isCellCentered_    [(i)]
#define g_isPeriodic(i)              g_isPeriodic_        [(i)]
#define r_(i,j)                      r__                  [(i) + r_s * (j)]
#define interpolationStencil_(i,j,k) \
                     interpolationStencil__[(i) + iS_s1 * (j) + iS_s2 * (k)]
#define useBackupRules_(i)           useBackupRules__     [(i)]
#define iw0_(i)                      iw0__                [(i)]
    Integer *g_boundaryCondition_   = g.boundaryCondition() .getDataPointer(),
            *g_indexRange_          = g.indexRange()        .getDataPointer(),
            *g_extendedIndexRange_  = g.extendedIndexRange().getDataPointer(),
            *g_isCellCentered_      = g.isCellCentered()    .getDataPointer(),
            *g_isPeriodic_          = g.isPeriodic()        .getDataPointer(),
            *interpolationStencil__ = interpolationStencil  .getDataPointer(),
            *useBackupRules__       = useBackupRules        .getDataPointer(),
            *interpolationWidth__   = &interpolationWidth(0,k1,k2,l);
//        *backupInterpolationWidth__ = &backupInterpolationWidth(0,k1,k2,l);
    Real    *g_gridSpacing_         = g.gridSpacing()       .getDataPointer(),
            *r__                    = r                     .getDataPointer();
    const Integer r_s = &r(base,1) - &r(base,0),
      iS_s1 = &interpolationStencil(base,1,0) -
              &interpolationStencil(base,0,0),
      iS_s2 = &interpolationStencil(base,0,1) -
              &interpolationStencil(base,0,0);
    r__                    = &r_(-base,0);
    interpolationStencil__ = &interpolationStencil_(-base,0,0);
    useBackupRules__       = &useBackupRules_(-base);
#endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES

    // wdh: we can interpolate from extended index range
    real rBound[3][2];
    Integer kd;
    for (kd=0; kd<3; kd++) 
    {
      rBound[kd][0]=a+(g_extendedIndexRange(0,kd)-g_indexRange(0,kd))*g_gridSpacing(kd);   
      rBound[kd][1]=b+(g_extendedIndexRange(1,kd)-g_indexRange(1,kd))*g_gridSpacing(kd);
    }

    for (Integer i=base; i<=bound; i++) {
#ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
        IntegerArray& iw0 = interpolationWidth;
//        IntegerArray& iw0 = useBackupRules_(i) ?
//          backupInterpolationWidth : interpolationWidth;
#else
        Integer* iw0__ = interpolationWidth__;
//        Integer* iw0__ = useBackupRules_(i) ?
//          backupInterpolationWidth__ : interpolationWidth__;
#endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
        for ( kd=0; kd<3; kd++) if (kd < numberOfDimensions) {
// *wdh     if (r_(i,kd) < a || r_(i,kd) > b) {
            if (r_(i,kd) < rBound[kd][0] || r_(i,kd) > rBound[kd][1]) {
                interpolationStencil_(i,0,kd) =
                interpolationStencil_(i,1,kd) = INTEGER_MAX;
            } else {
                Real rr = r_(i,kd) / g_gridSpacing(kd) + g_indexRange(0,kd);
                interpolationStencil_(i,0,kd) =
                  Integer(floor(rr - (Real).5 * iw0_(kd) +
                  (g_isCellCentered(kd) ? (Real).5 : (Real)1.)));
                interpolationStencil_(i,1,kd) =
                  Integer(floor(rr + (Real).5 * iw0_(kd) -
                  (g_isCellCentered(kd) ? (Real).5 : (Real)0.)));
                if (!g_isPeriodic(kd)) {
                    if (interpolationStencil_(i,0,kd) < g_extendedIndexRange(0,kd) &&
                      g_boundaryCondition(0,kd)) {
//                      Point is close to a BC side.  One-sided interpolation used.
                        interpolationStencil_(i,0,kd) = g_extendedIndexRange(0,kd);
                        interpolationStencil_(i,1,kd) = interpolationStencil_(i,0,kd)
                          + (iw0_(kd) - 1);
                    } // end if
                    if (interpolationStencil_(i,1,kd) > g_extendedIndexRange(1,kd) &&
                      g_boundaryCondition(1,kd)) {
//                      Point is close to a BC side.  One-sided interpolation used.
                        interpolationStencil_(i,1,kd) = g_extendedIndexRange(1,kd);
                        interpolationStencil_(i,0,kd) = interpolationStencil_(i,1,kd)
                          - (iw0_(kd) - 1);
                    } // end if
                } // end if
            } // end if
        } else if (kd <= interpolationStencil.getBound(2)) {
            interpolationStencil_(i,0,kd) = g_extendedIndexRange(0,kd);
            interpolationStencil_(i,1,kd) = g_extendedIndexRange(1,kd);
        } // end if, end for
    } // end for
#undef g_boundaryCondition
#undef g_gridSpacing
#undef g_indexRange
#undef g_extendedIndexRange
#undef g_isCellCentered
#undef g_isPeriodic
#undef interpolationStencil_
#undef useBackupRules_
#undef iw0_
#undef r_
}
void CompositeGridData::getInterpolationStencil(
  const MappedGrid&   g, // The unrefined grid corresponding to grid[k20].
  const Integer&      k10,
  const Integer&      k20,
  const RealArray&    r,
  const IntegerArray& interpolationStencil,
  const intArray& useBackupRules) {
    MappedGrid& g2 = grid[k20];
    const Real a = -(Real)2. * epsilon, b = (Real)1. - a;
    const Integer base = r.getBase(0), bound = r.getBound(0),
      k1 = componentGridNumber(k10), k2 = componentGridNumber(k20),
      l  = multigridLevelNumber(k10);

#ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
#define refinementFactor_(i,j)       refinementFactor     ((i),(j))
#define g_boundaryCondition(i,j)     g.boundaryCondition  ((i),(j))
#define g_discretizationWidth(i)     g.discretizationWidth((i))
#define g_indexRange(i,j)            g.indexRange         ((i),(j))
#define g_numberOfGhostPoints(i,j)   g.numberOfGhostPoints((i),(j))
#define g_isCellCentered(i)          g.isCellCentered     ((i))
#define g_isPeriodic(i)              g.isPeriodic         ((i))
#define g2_indexRange(i,j)           g2.indexRange        ((i),(j))
#define g2_extendedIndexRange(i,j)   g2.extendedIndexRange((i),(j))
#define g2_gridSpacing(i)            g2.gridSpacing       ((i))
#define g2_useGhostPoints            g2.useGhostPoints()
#define r_(i,j)                      r                    ((i),(j))
#define interpolationStencil_(i,j,k) interpolationStencil ((i),(j),(k))
#define useBackupRules_(i)           useBackupRules       ((i))
#define iw0_(i)                      iw0                  ((i),k1,k2,l)
#else
#define refinementFactor_(i,j)       refinementFactor__    [(i) + 3 * (j)]
#define g_boundaryCondition(i,j)     g_boundaryCondition_  [(i) + 2 * (j)]
#define g_discretizationWidth(i)     g_discretizationWidth_[(i)]
#define g_indexRange(i,j)            g_indexRange_         [(i) + 2 * (j)]
#define g_numberOfGhostPoints(i,j)   g_numberOfGhostPoints_[(i) + 2 * (j)]
#define g_isCellCentered(i)          g_isCellCentered_     [(i)]
#define g_isPeriodic(i)              g_isPeriodic_         [(i)]
#define g2_indexRange(i,j)           g2_indexRange_        [(i) + 2 * (j)]
#define g2_extendedIndexRange(i,j)   g2_extendedIndexRange_[(i) + 2 * (j)]
#define g2_gridSpacing(i)            g2_gridSpacing_       [(i)]
#define g2_useGhostPoints            g2_useGhostPoints_
#define r_(i,j)                      r__                   [(i) + r_s * (j)]
#define interpolationStencil_(i,j,k) \
                     interpolationStencil__[(i) + iS_s1 * (j) + iS_s2 * (k)]
#define useBackupRules_(i)           useBackupRules__      [(i)]
#define iw0_(i)                      iw0__                 [(i)]
    Integer *refinementFactor__     = refinementFactor       .getDataPointer(),
            *g_boundaryCondition_   = g.boundaryCondition()  .getDataPointer(),
            *g_discretizationWidth_ = g.discretizationWidth().getDataPointer(),
            *g_indexRange_          = g.indexRange()         .getDataPointer(),
            *g_numberOfGhostPoints_ = g.numberOfGhostPoints().getDataPointer(),
            *g_isCellCentered_      = g.isCellCentered()     .getDataPointer(),
            *g_isPeriodic_          = g.isPeriodic()         .getDataPointer(),
            *g2_indexRange_         = g2.indexRange()        .getDataPointer(),
            *g2_extendedIndexRange_ = g2.extendedIndexRange().getDataPointer(),
            *interpolationStencil__ = interpolationStencil   .getDataPointer(),
            *useBackupRules__       = useBackupRules         .getDataPointer(),
            *interpolationWidth__   = &interpolationWidth(0,k1,k2,l);
//        *backupInterpolationWidth__ = &backupInterpolationWidth(0,k1,k2,l);
    Logical g2_useGhostPoints_      = g2.useGhostPoints();
    Real    *g2_gridSpacing_        = g2.gridSpacing()       .getDataPointer(),
            *r__                    = r                      .getDataPointer();
    const Integer r_s = &r(base,1) - &r(base,0),
      iS_s1 = &interpolationStencil(base,1,0) -
              &interpolationStencil(base,0,0),
      iS_s2 = &interpolationStencil(base,0,1) -
              &interpolationStencil(base,0,0);
    r__                    = &r_(-base,0);
    interpolationStencil__ = &interpolationStencil_(-base,0,0);
    useBackupRules__       = &useBackupRules_(-base);
#endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES

#define g_extendedIndexRange(i,j) g_extendedIndexRange_[(j)][(i)]
    Integer g_extendedIndexRange_[3][2], kd, ks;
    for (kd=0; kd<3; kd++) {
        for (ks=0; ks<2; ks++) g_extendedIndexRange(ks,kd) =
          g_indexRange(ks,kd) * refinementFactor_(kd,k20);
        if (g_isCellCentered(kd) || g_isPeriodic(kd))
          g_extendedIndexRange(1,kd) += refinementFactor(kd,k20) - 1;
    } // end for
    for (kd=0; kd<numberOfDimensions; kd++) for (ks=0; ks<2; ks++)
      if (g_boundaryCondition(ks,kd) == 0 && g2_useGhostPoints)
        g_extendedIndexRange(ks,kd) =
          max0(g_extendedIndexRange(0,kd) - g_numberOfGhostPoints(0,kd),
          min0(g_extendedIndexRange(1,kd) + g_numberOfGhostPoints(1,kd),
          g_extendedIndexRange(ks,kd) +
          (2 * ks - 1) * (g_discretizationWidth(kd) - 1) / 2));


    // wdh: we can interpolate from extended index range
    real rBound[3][2];
    for (kd=0; kd<3; kd++) 
    {
      rBound[kd][0]=a+(g2_extendedIndexRange(0,kd)-g2_indexRange(0,kd))*g2_gridSpacing(kd);   
      rBound[kd][1]=b+(g2_extendedIndexRange(1,kd)-g2_indexRange(1,kd))*g2_gridSpacing(kd);
    }

    for (Integer i=base; i<=bound; i++) {
#ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
        IntegerArray& iw0 = interpolationWidth;
//         IntegerArray& iw0 = useBackupRules_(i) ?
//           backupInterpolationWidth : interpolationWidth;
#else
        Integer* iw0__ = interpolationWidth__;
//         Integer* iw0__ = useBackupRules_(i) ?
//           backupInterpolationWidth__ : interpolationWidth__;
#endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
        for (kd=0; kd<3; kd++) if (kd < numberOfDimensions) {
// *wdh     if (r_(i,kd) < a || r_(i,kd) > b) {
            if (r_(i,kd) < rBound[kd][0] || r_(i,kd) > rBound[kd][1]) {
                interpolationStencil_(i,0,kd) =
                interpolationStencil_(i,1,kd) = INTEGER_MAX;
            } else {
                Real rr = r_(i,kd) / g2_gridSpacing(kd) + g2_indexRange(0,kd);
                interpolationStencil_(i,0,kd) =
                  Integer(floor(rr - (Real).5 * iw0_(kd) +
                  (g_isCellCentered(kd) ? (Real).5 : (Real)1.)));
                interpolationStencil_(i,1,kd) =
                  Integer(floor(rr + (Real).5 * iw0_(kd) -
                  (g_isCellCentered(kd) ? (Real).5 : (Real)0.)));
                if (!g_isPeriodic(kd)) {
                    if (interpolationStencil_(i,0,kd) < g_extendedIndexRange(0,kd) &&
                      g_boundaryCondition(0,kd)) {
//                      Point is close to a BC side.  One-sided interpolation used.
                        interpolationStencil_(i,0,kd) = g_extendedIndexRange(0,kd);
                        interpolationStencil_(i,1,kd) = interpolationStencil_(i,0,kd)
                          + (iw0_(kd) - 1);
                    } // end if
                    if (interpolationStencil_(i,1,kd) > g_extendedIndexRange(1,kd) &&
                      g_boundaryCondition(1,kd)) {
//                      Point is close to a BC side.  One-sided interpolation used.
                        interpolationStencil_(i,1,kd) = g_extendedIndexRange(1,kd);
                        interpolationStencil_(i,0,kd) = interpolationStencil_(i,1,kd)
                          - (iw0_(kd) - 1);
                    } // end if
                } // end if
            } // end if
        } else if (kd <= interpolationStencil.getBound(2)) {
            interpolationStencil_(i,0,kd) = g_extendedIndexRange(0,kd);
            interpolationStencil_(i,1,kd) = g_extendedIndexRange(1,kd);
        } // end if, end for
    } // end for
#undef refinementFactor_
#undef g_boundaryCondition
#undef g_discretizationWidth
#undef g_indexRange
#undef g_extendedIndexRange
#undef g_numberOfGhostPoints
#undef g_isCellCentered
#undef g_isPeriodic
#undef g2_indexRange
#undef g2_gridSpacing
#undef g2_useGhostPoints
#undef interpolationStencil_
#undef useBackupRules_
#undef iw0_
#undef r_
}



Logical CompositeGridData::canInterpolate(
  const Integer&      k10,
  const Integer&      k20,
  const realArray&    r,
  const intArray& ok,
  const intArray& useBackupRules,
  const Logical       checkForOneSided) {
//
//  Determine whether points on grid k1 at r in the coordinates of grids k2
//  can be interpolated from grids k2.
//
    MappedGrid& g = grid[k20];
    Integer iv1[3], &i1=iv1[0], &i2=iv1[1], &i3=iv1[2], ks, kd, iab_[2*3];
    Logical isOneSided, oneSided[3][2], returnValue = LogicalTrue, invalid;
    IntegerArray iab2(1,2,3); RealArray rA(1,numberOfDimensions);
// *wdh* 980607    const Real a = -(Real)100. * epsilon, b = (Real)1. - a;
    const Real a = -(Real)2. * epsilon, b = (Real)1. - a;
    const Integer base = r.getBase(0), bound = r.getBound(0),
      k1 = componentGridNumber(k10), k2 = componentGridNumber(k20),
      l  = multigridLevelNumber(k10);

    assert( k10>=0 && k10<numberOfGrids);
    assert( k20>=0 && k20<numberOfGrids);

    assert( k1>=0 && k1<numberOfGrids);
    assert( k2>=0 && k2<numberOfGrids);
    

#define iab(i,j) iab_[(i) + 2 * (j)]
#ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
#define g_boundaryCondition(i,j)  g.boundaryCondition  ((i),(j))
#define g_dimension(i,j)          g.dimension          ((i),(j))
#define g_discretizationWidth(i)  g.discretizationWidth((i))
#define g_gridSpacing(i)          g.gridSpacing        ((i))
#define g_indexRange(i,j)         g.indexRange         ((i),(j))
#define g_extendedIndexRange(i,j) g.extendedIndexRange ((i),(j))
#define g_isCellCentered(i)       g.isCellCentered     ((i))
#define g_isPeriodic(i)           g.isPeriodic         ((i))
#define g_mask(i,j,k)             g.mask()             ((i),(j),(k))
#define r_(i,j)                   r                    ((i),(j))
#define useBackupRules_(i)        useBackupRules       ((i))
#define ok_(i)                    ok                   ((i))
#define iw0_(i)                   iw0                  ((i),k1,k2,l)
#define ov0_(i)                   ov0                  ((i),k1,k2,l)
#else
#define g_boundaryCondition(i,j)  g_boundaryCondition_  [(i) + 2 * (j)]
#define g_dimension(i,j)          g_dimension_          [(i) + 2 * (j)]
#define g_discretizationWidth(i)  g_discretizationWidth_[(i)]
#define g_gridSpacing(i)          g_gridSpacing_        [(i)]
#define g_indexRange(i,j)         g_indexRange_         [(i) + 2 * (j)]
#define g_extendedIndexRange(i,j) g_extendedIndexRange_ [(i) + 2 * (j)]
#define g_isCellCentered(i)       g_isCellCentered_     [(i)]
#define g_isPeriodic(i)           g_isPeriodic_         [(i)]
#define g_mask(i,j,k)             g_mask_               [(i)+i10*(j)+j10*(k)]
#define r_(i,j)                   r__                   [(i) + r_s * (j)]
#define useBackupRules_(i)        useBackupRules__      [(i)]
#define ok_(i)                    ok__                  [(i)]
#define iw0_(i)                   iw0__                 [(i)]
#define ov0_(i)                   ov0__                 [(i)]
    Integer *g_boundaryCondition_   = g.boundaryCondition()  .getDataPointer(),
            *g_dimension_           = g.dimension()          .getDataPointer(),
            *g_discretizationWidth_ = g.discretizationWidth().getDataPointer(),
            *g_indexRange_          = g.indexRange()         .getDataPointer(),
//            *g_extendedIndexRange_  = g.extendedIndexRange() .getDataPointer(),
            *g_isCellCentered_      = g.isCellCentered()     .getDataPointer(),
            *g_isPeriodic_          = g.isPeriodic()         .getDataPointer(),
            *g_mask_                = g.mask()               .getDataPointer(),
            *useBackupRules__       = useBackupRules         .getDataPointer(),
            *ok__                   = ok                     .getDataPointer(),
            *interpolationWidth__   = &interpolationWidth(0,k1,k2,l);
    
//        *backupInterpolationWidth__ = &backupInterpolationWidth(0,k1,k2,l);
    Real    *g_gridSpacing_         = g.gridSpacing()        .getDataPointer(),
            *r__                    = r                      .getDataPointer(),
            *interpolationOverlap__ = &interpolationOverlap(0,k1,k2,l);
//      *backupInterpolationOverlap__ = &backupInterpolationOverlap(0,k1,k2,l);
    const Integer i10 = g_dimension(1,0) - g_dimension(0,0) + 1,
           j10 = i10 * (g_dimension(1,1) - g_dimension(0,1) + 1),
      r_s = &r(base,1) - &r(base,0);
    g_mask_ = &g_mask(-g_dimension(0,0),-g_dimension(0,1),-g_dimension(0,2));
    r__               = &r_(-base,0);
    ok__              = &ok_(-base);
    useBackupRules__  = &useBackupRules_(-base);
#endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
    const Integer *g_I1 = g.I1(), *g_I2 = g.I2(), *g_I3 = g.I3();
    
    // wdh: we can interpolate from extended index range
    real rBound[3][2];
    for (kd=0; kd<3; kd++) 
    {
      rBound[kd][0]=a+(g.extendedRange(0,kd)-g_indexRange(0,kd))*g_gridSpacing(kd);   
      rBound[kd][1]=b+(g.extendedRange(1,kd)-g_indexRange(1,kd))*g_gridSpacing(kd);
    }

    for (Integer i=base; i<=bound; i++) if (ok_(i)) {
//
//      Determine the stencil of points to check.
//
#ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
        IntegerArray& iw0 = interpolationWidth;
        RealArray& ov0 = interpolationOverlap;
//         IntegerArray& iw0 = useBackupRules_(i) ?
//           backupInterpolationWidth : interpolationWidth;
//         RealArray& ov0 = useBackupRules_(i) ?
//           backupInterpolationOverlap : interpolationOverlap;
#else
        Integer* iw0__ = interpolationWidth__;
        Real* ov0__ = interpolationOverlap__;
//         Integer* iw0__ = useBackupRules_(i) ?
//           backupInterpolationWidth__ : interpolationWidth__;
//         Real* ov0__ = useBackupRules_(i) ?
//           backupInterpolationOverlap__ : interpolationOverlap__;
#endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
        invalid = isOneSided = LogicalFalse;
        for (kd=0; kd<3; kd++) 
        {
            oneSided[kd][0] = oneSided[kd][1] = LogicalFalse;
            if (kd < numberOfDimensions) 
            {
	      // *wdh if (invalid = r_(i,kd) < a || r_(i,kd) > b) break;
	        if (invalid = r_(i,kd) < rBound[kd][0] || r_(i,kd) > rBound[kd][1]) break;
                Real rr = r_(i,kd) / g_gridSpacing(kd) + g_indexRange(0,kd);

                // real overlap=ov0_(kd); // *wdh*
		
                iab(0,kd) = Integer(floor(rr - ov0_(kd) -
                  (g_isCellCentered(kd) ? (Real).5 : (Real)0.)));
                iab(1,kd) = Integer(floor(rr + ov0_(kd) +
                  (g_isCellCentered(kd) ? (Real).5 : (Real)1.)));
                if (!g_isPeriodic(kd)) {
                    if (iab(0,kd) < g.extendedRange(0,kd)) {
//                      Check if point is too close to an interpolated side.
                        if (invalid = !g_boundaryCondition(0,kd)) break;
//                      One-sided interpolation is used close to a boundary.
                        isOneSided = oneSided[kd][0] = LogicalTrue;
                        iab(0,kd) = g.extendedRange(0,kd);
                        iab(1,kd) = iab(0,kd) +
                          Integer(floor((Real).5 * iw0_(kd) + ov0_(kd) + (Real).5));
                    } // end if
                    if (iab(1,kd) > g.extendedRange(1,kd)) {
//                      Check if point is too close to an interpolated side.
                        if (invalid = !g_boundaryCondition(1,kd)) break;
//                      One-sided interpolation is used close to a boundary.
                        isOneSided = oneSided[kd][1] = LogicalTrue;
                        iab(1,kd) = g.extendedRange(1,kd);
                        iab(0,kd) = iab(1,kd) -
                          Integer(floor((Real).5 * iw0_(kd) + ov0_(kd) + (Real).5));
                    } // end if
                } // end if
            } else {
                iab(0,kd) = g.extendedRange(0,kd);
                iab(1,kd) = g.extendedRange(1,kd);
            } // end if
        } // end for
//
//      Check that all points in the stencil are either discretization points
//      or interpolation points.  Backup discretization points and backup
//      interpolation points are also allowed.
//
        if (!invalid) COMPOSITE_GRID_FOR_3(iab_, i1, i2, i3)
          if (invalid = invalid ||
            !(g_mask(g_I1[i1],g_I2[i2],g_I3[i3]) & ISusedPoint)) break;

        if (!invalid && checkForOneSided && isOneSided) {
//
//          Check for one-sided interpolation from BC points
//          that interpolate from the interior of another grid.
//
//          Find the interpolation stencil.
            for (kd=0; kd<numberOfDimensions; kd++) rA(0,kd) = r_(i,kd);
            getInterpolationStencil(k10, k20, rA, iab2, useBackupRules);

#ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
#define     iab2_(i,j,k) iab2((i),(j),(k))
#else
#define     iab2_(i,j,k) iab2__[(j) + 2 * (k)]
#endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
            Integer* iab2__ = iab2.getDataPointer();
            for (kd=0; kd<3; kd++) {
                for (ks=0; ks<2; ks++) if (oneSided[kd][ks]) {
                    Integer iab21=iab2_(0,0,kd), iab22=iab2_(0,1,kd);
                    // *wdh* 021205: check added getInterpolationStencil will return this value for bogus pts
                    if( iab21==INT_MAX ) 
		    {
		      invalid=true;
		      break;
		    }
//                  Restrict the interpolation stencil to points that could be
//                  boundary discretization points of side (kd,ks) of the grid.
                    if (ks == 0) {
                      iab2_(0,0,kd) = g.extendedRange(0,kd);
                      iab2_(0,1,kd) = min0(
                        iab2_(0,1,kd),
                        iab2_(0,0,kd) +
                        (g_discretizationWidth(kd) - 1) / 2 - 1);
                    } else {
                      iab2_(0,1,kd) = g.extendedRange(1,kd);
                      iab2_(0,0,kd) = max0(iab2_(0,0,kd), iab2_(0,1,kd) -
                         (g_discretizationWidth(kd) - 1) / 2 + 1);
                    } // end if
//
//                  Check that all points in the stencil are either
//                  discretization points or interpolation points that are not
//                  interpolated one-sided from another grid.  Backup
//                  discretization points and backup interpolation points that
//                  are not interpolated one-sided from another grid are also
//                  allowed.
//
/* ------ *wdh* 980702
                    COMPOSITE_GRID_FOR_3(iab2__, i1, i2, i3)
                      if (invalid = invalid ||
                        g_mask(g_I1[i1],g_I2[i2],g_I3[i3]) &
                        ISinteriorBoundaryPoint) break;
                    if (invalid) break;
------- */
                    COMPOSITE_GRID_FOR_3(iab2__, i1, i2, i3)
		    {
                      if (invalid = invalid || g_mask(g_I1[i1],g_I2[i2],g_I3[i3]) & ISinteriorBoundaryPoint)
		      {
                        // Make sure that we are not too close to an the interpolation point
			real rDist=0.;
                        real cellCenterederedOffset=g_isCellCentered(kd) ? .5 : 0.;
			for( int dir=0; dir<numberOfDimensions; dir++ )
			  rDist=max(rDist,fabs( r_(i,dir)/g_gridSpacing(dir)
						-(iv1[dir]+cellCenterederedOffset-g_indexRange(Start,dir))));
			if( rDist > ov0_(0) )  // use ov_(0) as the minimum overlap. Normally=.5
			{
			  // printf("CompositeGrid::canInterpolate: near an interior boundary point but rDist=%e"
                          //       ", ov=%6.2e, so this point is ok! \n",rDist,ov0_(0));
                          invalid=FALSE;  // this point is ok after all
			}
                        else
			  break;
		      }
                      if (invalid) break;
		    }
		    
//                  Restore the interpolation stencil;
                    iab2_(0,0,kd) = iab21;
                    iab2_(0,1,kd) = iab22;
                } // end if, end for
                if (invalid) break;
            } // end for
        } // end if

        if (invalid) ok_(i) = returnValue = LogicalFalse;

    } else {
        returnValue = LogicalFalse;
    } // end if, end for
    return returnValue;
#undef iab
#undef g_boundaryCondition
#undef g_dimension
#undef g_discretizationWidth
#undef g_gridSpacing
#undef g_indexRange
#undef g_extendedIndexRange
#undef g_isCellCentered
#undef g_isPeriodic
#undef g_mask
#undef r_
#undef useBackupRules_
#undef ok_
#undef iw0_
#undef ov0_
#undef iab2_
}


#if 0
Logical CompositeGridData::
canInterpolate(
  const Integer&      k10,
  const Integer&      k20,
  const realArray&    r,
  const intArray& ok,
  const intArray& useBackupRules,
  const Logical       checkForOneSided) 
{
  return canInterpolate(k10,k20,r.getLocalArray(),ok.getLocalArray(), useBackupRules.getLocalArray(), checkForOneSided );
}
#endif


// Logical CompositeGridData::canInterpolate(
//   const MappedGrid&   g, // The unrefined grid corresponding to grid[k20].
//   CompositeMask&      g_mask, // Masks on k2 (possibly) and its siblings.
//   const Integer&      k10,
//   const Integer&      k20,
//   const RealArray&    r,
//   const LogicalArray& ok,
//   const LogicalArray& useBackupRules,
//   const Logical       checkForOneSided) {
// //
// //  Determine whether points on grid k1 at r in the coordinates of grids k2
// //  can be interpolated from grids k2.
// //
//     MappedGrid& g2 = grid[k20];
//     Integer iv1[3], &i1=iv1[0], &i2=iv1[1], &i3=iv1[2], ks, kd, iab_[2*3];
//     Logical isOneSided, oneSided[3][2], returnValue = LogicalTrue, invalid;
//     IntegerArray iab2(1,2,3); RealArray rA(1,numberOfDimensions);
//     const Integer base = r.getBase(0), bound = r.getBound(0),
//       k1 = componentGridNumber(k10), k2 = componentGridNumber(k20),
//       l  = multigridLevelNumber(k10);

// #define iab(i,j) iab_[(i) + 2 * (j)]
// #ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
// #define refinementFactor_(i,j)     refinementFactor     ((i),(j))
// #define g_boundaryCondition(i,j)   g.boundaryCondition  ((i),(j))
// #define g_discretizationWidth(i)   g.discretizationWidth((i))
// #define g_indexRange(i,j)          g.indexRange         ((i),(j))
// #define g_numberOfGhostPoints(i,j) g.numberOfGhostPoints((i),(j))
// #define g_isCellCentered(i)        g.isCellCentered     ((i))
// #define g_isPeriodic(i)            g.isPeriodic         ((i))
// #define g2_indexRange(i,j)         g2.indexRange        ((i),(j))
// #define g2_gridSpacing(i)          g2.gridSpacing       ((i))
// #define g2_useGhostPoints          g2.useGhostPoints()
// #define r_(i,j)                    r                    ((i),(j))
// #define useBackupRules_(i)         useBackupRules       ((i))
// #define ok_(i)                     ok                   ((i))
// #define iw0_(i)                    iw0                  ((i),k1,k2,l)
// #define ov0_(i)                    ov0                  ((i),k1,k2,l)
// #else
// #define refinementFactor_(i,j)     refinementFactor__    [(i) + 3 * (j)]
// #define g_boundaryCondition(i,j)   g_boundaryCondition_  [(i) + 2 * (j)]
// #define g_discretizationWidth(i)   g_discretizationWidth_[(i)]
// #define g_indexRange(i,j)          g_indexRange_         [(i) + 2 * (j)]
// #define g_numberOfGhostPoints(i,j) g_numberOfGhostPoints_[(i) + 2 * (j)]
// #define g_isCellCentered(i)        g_isCellCentered_     [(i)]
// #define g_isPeriodic(i)            g_isPeriodic_         [(i)]
// #define g2_indexRange(i,j)         g2_indexRange_        [(i) + 2 * (j)]
// #define g2_gridSpacing(i)          g2_gridSpacing_       [(i)]
// #define g2_useGhostPoints          g2_useGhostPoints_
// #define r_(i,j)                    r__                   [(i) + r_s * (j)]
// #define useBackupRules_(i)         useBackupRules__      [(i)]
// #define ok_(i)                     ok__                  [(i)]
// #define iw0_(i)                    iw0__                 [(i)]
// #define ov0_(i)                    ov0__                 [(i)]
//     Integer *refinementFactor__     = refinementFactor       .getDataPointer(),
//             *g_boundaryCondition_   = g.boundaryCondition()  .getDataPointer(),
//             *g_discretizationWidth_ = g.discretizationWidth().getDataPointer(),
//             *g_indexRange_          = g.indexRange()         .getDataPointer(),
//             *g_numberOfGhostPoints_ = g.numberOfGhostPoints().getDataPointer(),
//             *g_isCellCentered_      = g.isCellCentered()     .getDataPointer(),
//             *g_isPeriodic_          = g.isPeriodic()         .getDataPointer(),
//             *g2_indexRange_         = g2.indexRange()        .getDataPointer(),
//             *useBackupRules__       = useBackupRules         .getDataPointer(),
//             *ok__                   = ok                     .getDataPointer(),
//             *interpolationWidth__   = &interpolationWidth(0,k1,k2,l),
//         *backupInterpolationWidth__ = &backupInterpolationWidth(0,k1,k2,l);
//     Logical g2_useGhostPoints_      = g2.useGhostPoints();
//     Real    *g2_gridSpacing_        = g2.gridSpacing()       .getDataPointer(),
//             *r__                    = r                      .getDataPointer(),
//             *interpolationOverlap__ = &interpolationOverlap(0,k1,k2,l),
//       *backupInterpolationOverlap__ = &backupInterpolationOverlap(0,k1,k2,l);
//     const Integer r_s = &r(base,1) - &r(base,0);
//     r__               = &r_(-base,0);
//     ok__              = &ok_(-base);
//     useBackupRules__  = &useBackupRules_(-base);
// #endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES

// #define g_extendedIndexRange(i,j) g_extendedIndexRange_[(j)][(i)]
//     Integer g_extendedIndexRange_[3][2];
//     for (kd=0; kd<3; kd++) {
//         for (ks=0; ks<2; ks++) g_extendedIndexRange(ks,kd) =
//           g_indexRange(ks,kd) * refinementFactor_(kd,k20);
//         if (g_isCellCentered(kd) || g_isPeriodic(kd))
//           g_extendedIndexRange(1,kd) += refinementFactor(kd,k20) - 1;
//     } // end for
//     for (kd=0; kd<numberOfDimensions; kd++) for (ks=0; ks<2; ks++)
//       if (g_boundaryCondition(ks,kd) == 0 && g2_useGhostPoints)
//         g_extendedIndexRange(ks,kd) =
//           max0(g_extendedIndexRange(0,kd) - g_numberOfGhostPoints(0,kd),
//           min0(g_extendedIndexRange(1,kd) + g_numberOfGhostPoints(1,kd),
//           g_extendedIndexRange(ks,kd) +
//           (2 * ks - 1) * (g_discretizationWidth(kd) - 1) / 2));

//     for (Integer i=base; i<=bound; i++) if (ok_(i)) {
// //
// //      Determine the stencil of points to check.
// //
// #ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
//         IntegerArray& iw0 = useBackupRules_(i) ?
//           backupInterpolationWidth : interpolationWidth;
//         RealArray& ov0 = useBackupRules_(i) ?
//           backupInterpolationOverlap : interpolationOverlap;
// #else
//         Integer* iw0__ = useBackupRules_(i) ?
//           backupInterpolationWidth__ : interpolationWidth__;
//         Real* ov0__ = useBackupRules_(i) ?
//           backupInterpolationOverlap__ : interpolationOverlap__;
// #endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
//         invalid = isOneSided = LogicalFalse;
//         for (kd=0; kd<3; kd++) {
//             oneSided[kd][0] = oneSided[kd][1] = LogicalFalse;
//             if (kd < numberOfDimensions) {
//                 Real rr = r_(i,kd) / g2_gridSpacing(kd) + g2_indexRange(0,kd);
//                 iab(0,kd) = Integer(floor(rr - ov0_(kd) -
//                   (g_isCellCentered(kd) ? (Real).5 : (Real)0.)));
//                 iab(1,kd) = Integer(floor(rr + ov0_(kd) +
//                   (g_isCellCentered(kd) ? (Real).5 : (Real)1.)));
//                 if (!g_isPeriodic(kd)) {
//                     if (iab(0,kd) < g_extendedIndexRange(0,kd)) {
// //                      Check if point is too close to an interpolated side.
//                         if (invalid = !g_boundaryCondition(0,kd)) break;
// //                      One-sided interpolation is used close to a boundary.
//                         isOneSided = oneSided[kd][0] = LogicalTrue;
//                         iab(0,kd) = g_extendedIndexRange(0,kd);
//                         iab(1,kd) = iab(0,kd) +
//                           Integer(floor((Real).5 * iw0_(kd) + ov0_(kd) + (Real).5));
//                     } // end if
//                     if (iab(1,kd) > g_extendedIndexRange(1,kd)) {
// //                      Check if point is too close to an interpolated side.
//                         if (invalid = !g_boundaryCondition(1,kd)) break;
// //                      One-sided interpolation is used close to a boundary.
//                         isOneSided = oneSided[kd][1] = LogicalTrue;
//                         iab(1,kd) = g_extendedIndexRange(1,kd);
//                         iab(0,kd) = iab(1,kd) -
//                           Integer(floor((Real).5 * iw0_(kd) + ov0_(kd) + (Real).5));
//                     } // end if
//                 } // end if
//             } else {
//                 iab(0,kd) = g_extendedIndexRange(0,kd);
//                 iab(1,kd) = g_extendedIndexRange(1,kd);
//             } // end if
//         } // end for
// //
// //      Check that all points in the stencil are either discretization points
// //      or interpolation points.  Backup discretization points and backup
// //      interpolation points are also allowed.
// //
//         if (!invalid) COMPOSITE_GRID_FOR_3(iab_, i1, i2, i3)
//           if (invalid = invalid ||
//             !((const Integer&)g_mask(i1,i2,i3) & ISusedPoint)) break;

//         if (!invalid && checkForOneSided && isOneSided) {
// //
// //          Check for one-sided interpolation from BC points
// //          that interpolate from the interior of another grid.
// //
// //          Find the interpolation stencil.
//             for (kd=0; kd<numberOfDimensions; kd++) rA(0,kd) = r_(i,kd);
//             getInterpolationStencil(g, k10, k20, rA, iab2, useBackupRules);

// #ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
// #define     iab2_(i,j,k) iab2((i),(j),(k))
// #else
// #define     iab2_(i,j,k) iab2__[(j) + 2 * (k)]
// #endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
//             Integer* iab2__ = iab2.getDataPointer();
//             for (kd=0; kd<3; kd++) {
//                 for (ks=0; ks<2; ks++) if (oneSided[kd][ks]) {
//                     Integer iab21=iab2_(0,0,kd), iab22=iab2_(0,1,kd);
// //                  Restrict the interpolation stencil to points that could be
// //                  boundary discretization points of side (kd,ks) of the grid.
//                     if (ks == 0) {
//                       iab2_(0,0,kd) = g_extendedIndexRange(0,kd);
//                       iab2_(0,1,kd) = min0(
//                         iab2_(0,1,kd),
//                         iab2_(0,0,kd) +
//                         (g_discretizationWidth(kd) - 1) / 2 - 1);
//                     } else {
//                       iab2_(0,1,kd) = g_extendedIndexRange(1,kd);
//                       iab2_(0,0,kd) = max0(iab2_(0,0,kd), iab2_(0,1,kd) -
//                          (g_discretizationWidth(kd) - 1) / 2 + 1);
//                     } // end if
// //
// //                  Check that all points in the stencil are either
// //                  discretization points or interpolation points that are not
// //                  interpolated one-sided from another grid.  Backup
// //                  discretization points and backup interpolation points that
// //                  are not interpolated one-sided from another grid are also
// //                  allowed.
// //
//                     COMPOSITE_GRID_FOR_3(iab2__, i1, i2, i3)
//                       if (invalid = invalid || (const Integer&)
//                         g_mask(i1,i2,i3) & ISinteriorBoundaryPoint) break;
//                     if (invalid) break;

// //                  Restore the interpolation stencil;
//                     iab2_(0,0,kd) = iab21;
//                     iab2_(0,1,kd) = iab22;
//                 } // end if, end for
//                 if (invalid) break;
//             } // end for
//         } // end if

//         if (invalid) ok_(i) = returnValue = LogicalFalse;

//     } else {
//         returnValue = LogicalFalse;
//     } // end if, end for
//     return returnValue;
// #undef iab
// #undef refinementFactor_
// #undef g_boundaryCondition
// #undef g_discretizationWidth
// #undef g_indexRange
// #undef g_numberOfGhostPoints
// #undef g_extendedIndexRange
// #undef g_isCellCentered
// #undef g_isPeriodic
// #undef g2_indexRange
// #undef g2_gridSpacing
// #undef g2_useGhostPoints
// #undef r_
// #undef useBackupRules_
// #undef ok_
// #undef iw0_
// #undef ov0_
// #undef iab2_
// }
//
// Check if these boundary discretization points of this grid
// lie at least epsilon inside the parameter space of grid g2.
//
// void CompositeGridData::isInteriorBoundaryPoint(
//   const Integer&      k1,
//   const Integer&      k2,
//   const IntegerArray& i1,
//   const RealArray&    r2,
//   const LogicalArray& ok) {
//     MappedGrid &g1 = grid[k1], &g2 = grid[k2];
//     RealArray r, x; r.redim(r2); x.redim(r2); r = (Real).5;
//     IntegerArray numberOfSides; numberOfSides.redim(ok);
//     Logical areAllBoundaryPoints = LogicalTrue;
//     const Integer base = i1.getBase(0), bound = i1.getBound(0);
//     IntegerArray i1b = i1;
// #ifdef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
// #define g1_boundaryCondition(i,j) g1.boundaryCondition((i),(j))
// #define g1_discretizationWidth(i) g1.discretizationWidth((i))
// #define g1_gridIndexRange(i,j)    g1.gridIndexRange((i),(j))
// #define g1_gridSpacing(i)         g1.gridSpacing((i))
// #define g1_indexRange(i,j)        g1.indexRange((i),(j))
// #define g1_isCellCentered(i)      g1.isCellCentered((i))
// #define g2_boundaryCondition(i,j) g2.boundaryCondition((i),(j))
// #define i1b_(i,j)                 i1b((i),(j))
// #define r_(i,j)                   r((i),(j))
// #define ok_(i)                    ok((i))
// #define numberOfSides_(i) numberOfSides((i))
// #else
// #define g1_boundaryCondition(i,j) g1_boundaryCondition__  [(i) + 2     * (j)]
// #define g1_discretizationWidth(i) g1_discretizationWidth__[(i)              ]
// #define g1_gridIndexRange(i,j)    g1_gridIndexRange__     [(i) + 2     * (j)]
// #define g1_gridSpacing(i)         g1_gridSpacing__        [(i)              ]
// #define g1_indexRange(i,j)        g1_indexRange__         [(i) + 2     * (j)]
// #define g1_isCellCentered(i)      g1_isCellCentered__     [(i)              ]
// #define g2_boundaryCondition(i,j) g2_boundaryCondition_[(i) + 2     * (j)]
// #define i1b_(i,j)                 i1b__                [(i) + i1b_s * (j)]
// #define r_(i,j)                   r__                  [(i) + r_s   * (j)]
// #define ok_(i)                    ok__                 [(i)              ]
// #define numberOfSides_(i)         numberOfSides__      [(i)              ]
//     Integer   *g1_boundaryCondition__ = g1.boundaryCondition().getDataPointer(),
//           *g1_discretizationWidth__ = g1.discretizationWidth().getDataPointer(),
//               *g1_gridIndexRange__    = g1.gridIndexRange()   .getDataPointer(),
//               *g1_indexRange__        = g1.indexRange()       .getDataPointer(),
//               *g1_isCellCentered__    = g1.isCellCentered()   .getDataPointer(),
//               *g2_boundaryCondition_  = g2.boundaryCondition().getDataPointer(),
//               *i1b__                  = i1b                   .getDataPointer();
//     Real      *g1_gridSpacing__       = g1.gridSpacing()      .getDataPointer(),
//               *r__                    = r                     .getDataPointer();
//     LogicalAE *ok__                   = ok                    .getDataPointer(),
//               *numberOfSides__        = numberOfSides         .getDataPointer();
//     const Integer  i1b_s = &i1b(base,1) - &i1b(base,0),
//                    r_s   = &r(base,1)   - &r(base,0);
//     i1b__           = &i1b_(-base,0);
//     r__             = &r_  (-base,0);
//     ok__            = &ok_ (-base);
//     numberOfSides__ = &numberOfSides_(-base);
// #endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES

//     Integer i;
//     for (i=base; i<=bound; i++) if (ok_(i)) {
// //
// //      Find out how many sides of this grid the points lie on.
// //
//         numberOfSides_(i) = 0;
//         for (Integer kd=0; kd<numberOfDimensions; kd++) {
//             r_(i,kd) = (i1b_(i,kd) - g1_gridIndexRange(0,kd)) *
//               g1_gridSpacing(kd);
//             if (g1_isCellCentered(kd))
//               r_(i,kd) += (Real).5 * g1_gridSpacing(kd);
//             if (g1_boundaryCondition(0,kd) > 0 &&
//               i1b_(i,kd) < g1_indexRange(0,kd) +
//               (g1_discretizationWidth(kd) - 1) / 2) {
// //              The point is on the left side.
//                 numberOfSides_(i)++;
//                 if (g1_isCellCentered(kd) ||
//                   i1b_(i,kd) != g1_gridIndexRange(0,kd)) {
// //                  This is not a point on the boundary.
// //                  Compute the corresponding point on the boundary.
//                     i1b_(i,kd) = g1_gridIndexRange(0,kd);
//                     r_(i,kd) = (Real)0.;
//                     areAllBoundaryPoints = LogicalFalse;
//                 } // end if
//             } // end if
//             if (g1_boundaryCondition(1,kd) > 0 &&
//               i1b_(i,kd) > g1_indexRange(1,kd) -
//               (g1_discretizationWidth(kd) - 1) / 2) {
// //              The point is on the right side.
//                 numberOfSides_(i)++;
//                 if (g1_isCellCentered(kd) ||
//                   i1b_(i,kd) != g1_gridIndexRange(1,kd)) {
// //                  This is not a point on the boundary.
// //                  Compute the corresponding point on the boundary.
//                     i1b_(i,kd) = g1_gridIndexRange(1,kd);
//                     r_(i,kd) = (Real)1.;
//                     areAllBoundaryPoints = LogicalFalse;
//                 } // end if
//             } // end if
//         } // end for
//     } // end if, end for

//     if (areAllBoundaryPoints) {
// //
// //      All points lie on the boundary of this grid, and r2
// //      contains their coordinates in the parameter space of g2.
// //
//         r.reference(r2);

//     } else {
// //
// //      Some points do not lie on the boundary of this grid.
// //      Compute the corresponding boundary points and their
// //      coordinates in the parameter space of grid g2.
// //
//         g1.mapping().map(r, x);
//         adjustBoundary(k1, baseGridNumber(k2), i1b, x);
//         g2.mapping().inverseMap(x, r = r2);
//         const Range p(base,bound);
//         LogicalArray ok1(p); ok1 = r(p,0) != (Real)10.;
//         where (ok) ok(p) = ok1;
//     } // end if

// #ifndef DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
//     r__ = r.getDataPointer();
//     r__ = &r_(-base,0);
// #endif // DO_NOT_OPTIMIZE_SCALAR_ARRAY_REFERENCES
//     for (i=base; i<=bound; i++) if (ok_(i)) {
// //      Check if the points interpolate from the interior of grid g2.
// //***** Do not consider interior boundary conditions (e.g. polar axis).
//         for (Integer kd=0; kd<numberOfDimensions; kd++)
//           if ((g2_boundaryCondition(0,kd) > 0 &&
//                          r_(i,kd) < epsilon) ||
//               (g2_boundaryCondition(1,kd) > 0 &&
//               (Real)1. - r_(i,kd) < epsilon)) numberOfSides_(i)--;
// //
// //      The point lies in the interior of grid g2 only if numberOfSides > 0.
// //
//         if (numberOfSides_(i) <= 0) ok_(i) = LogicalFalse;
//     } // end if, end for
// #undef g1_boundaryCondition
// #undef g1_discretizationWidth
// #undef g1_gridIndexRange
// #undef g1_gridSpacing
// #undef g1_indexRange
// #undef g1_isCellCentered
// #undef g2_boundaryCondition
// #undef i1b_
// #undef r_
// #undef ok_
// #undef numberOfSides_
// }

// void CompositeGridData::adjustBoundary(
//   const Integer&      k1,
//   const Integer&      k2,
//   const IntegerArray& i1,
//   const RealArray&    x) {
// //
// //  Adjust the position x of points i1 of grid k1 interpolated from
// //  base grid k2 to take into account mismatch between shared boundaries.
// //
//     if (boundaryAdjustment.getNumberOfElements()) {
//         BoundaryAdjustmentArray& bA12 = boundaryAdjustment(k1,k2);
//         if (bA12.getNumberOfElements()) {
//             const Integer base = i1.getBase(0), bound = i1.getBound(0);
//             Integer jv[3], &j1=jv[0], &j2=jv[1], &j3=jv[2], kd1, kd2, ks1;
//             RealArray x0 = x;
//             for (kd1=numberOfDimensions; kd1<3; kd1++)
//               jv[kd1] = grid[k1].indexRange(0,kd1);
//             for (kd1=0; kd1<numberOfDimensions; kd1++)
//               for (ks1=0; ks1<2; ks1++) {
//                 BoundaryAdjustment& bA = bA12(ks1,kd1);
//                 if ((bA.computedGeometry & (THEinverseMap | THEmask)) ==
//                   (THEinverseMap | THEmask)) {
//                     jv[kd1] = bA.boundaryAdjustment.getBase(kd1);
//                     for (Integer i=base; i<=bound; i++) {
//                         for (kd2=0; kd2<numberOfDimensions; kd2++)
//                           if (kd2 != kd1) jv[kd2] = i1(i,kd2);
//                         Real dot = (Real)0.;
//                         for (kd2=0; kd2<numberOfDimensions; kd2++)
//                           dot += bA.acrossGrid(j1,j2,j3,kd2) *
//                             (bA.oppositeBoundary(j1,j2,j3,kd2) - x0(i,kd2));
//                         for (kd2=0; kd2<numberOfDimensions; kd2++)
//                           x(i,kd2) += dot * bA.boundaryAdjustment(j1,j2,j3,kd2);
//                     } // end for
//                 } // end if
//             } // end for, end for
//         } // end if
//     } // end if
// }

Integer CompositeGridData::updateCollection(
  const Integer&                   what,
  Integer&                         numberOfCollections,
#ifdef USE_STL
  RCVector<CompositeGrid>&         list,
  RCVector<GridCollection>&        gridCollectionList,
  RCVector<GenericGridCollection>& genericGridCollectionList,
#else
  ListOfCompositeGrid&             list,
  ListOfGridCollection&            gridCollectionList,
  ListOfGenericGridCollection&     genericGridCollectionList,
#endif // USE_STL
  IntegerArray&                    number) {
//  Fix up the length of list.
    numberOfCollections = numberOfGrids > 0 ? max(number) + 1 : 0;
#ifdef STL
    if (list.size() > numberOfCollections)
      list.erase(list.begin() + numberOfCollections, list.end());
    if (gridCollectionList.size() > numberOfCollections)
      gridCollectionList.erase(
        gridCollectionList.begin() + numberOfCollections,
        gridCollectionList.end());
    if (genericGridCollectionList.size() > numberOfCollections)
      genericGridCollectionList.erase(
        genericGridCollectionList.begin() + numberOfCollections,
        genericGridCollectionList.end());
#else
    while (list.getLength() > numberOfCollections) list.deleteElement();
    while (gridCollectionList.getLength() > numberOfCollections)
      gridCollectionList.deleteElement();
    while (genericGridCollectionList.getLength() > numberOfCollections)
      genericGridCollectionList.deleteElement();
#endif // STL
    if (numberOfCollections) {
//      Fill lists with appropriately-constructed CompositeGrids.
        IntegerArray nG(numberOfCollections); nG = 0; Integer k, i;
        for (k=0; k<numberOfGrids; k++) nG(number(k))++;
        for (i=0; i<numberOfCollections; i++) {
            if (i < list.getLength())
              list[i].setNumberOfDimensionsAndGrids(numberOfDimensions, nG(i));
#ifdef USE_STL
            else list.push_back(CompositeGrid(numberOfDimensions, nG(i)));
            if (i < gridCollectionList.size())
              gridCollectionList[i].reference(list[i]);
            else gridCollectionList.push_back(list[i]);
            if (i < genericGridCollectionList.size())
              genericGridCollectionList[i].reference(list[i]);
            else genericGridCollectionList.push_back(list[i]);
#else
            else list.addElement(CompositeGrid(numberOfDimensions, nG(i)));
            if (i < gridCollectionList.getLength())
              gridCollectionList[i].reference(list[i]);
            else gridCollectionList.addElement(list[i]);
            if (i < genericGridCollectionList.getLength())
              genericGridCollectionList[i].reference(list[i]);
            else genericGridCollectionList.addElement(list[i]);
#endif // USE_STL
        } // end for
        GridCollectionData::updateCollection(
          what, numberOfCollections, gridCollectionList,
          genericGridCollectionList, number);
        for (i=0; i<numberOfCollections; i++)
          list[i].setNumberOfDimensionsAndGrids(numberOfDimensions, nG(i));
        for (nG=0, k=0; k<numberOfGrids; k++) {
            const Integer j = nG(i = number(k))++;
            const Range three = 3, 
// *wdh	991023        allGrids = numberOfComponentGrids,
// *wdh 991023         allGridsPlusOne = numberOfComponentGrids + 1;
	      allGrids = list[i].numberOfComponentGrids(),
              allGridsPlusOne = list[i].numberOfComponentGrids() + 1;
            list[i].numberOfCompleteMultigridLevels() =
                    i < numberOfCompleteMultigridLevels ? 1 : 0;
            list[i].epsilon() = epsilon;
            list[i].numberOfInterpolationPoints(j) =
                    numberOfInterpolationPoints(k);
//            list[i].numberOfInterpoleePoints(j) =
//                    numberOfInterpoleePoints(k);
            list[i].interpolationIsAllExplicit() = interpolationIsAllExplicit;
            list[i].interpolationIsAllImplicit() = interpolationIsAllImplicit;
            list[i].interpolationIsImplicit(allGrids,allGrids,0) =
                    interpolationIsImplicit(allGrids,allGrids,number(k));
//            list[i].backupInterpolationIsImplicit(allGrids,allGrids,0) =
//                    backupInterpolationIsImplicit(allGrids,allGrids,number(k));
            list[i].interpolationWidth(three,allGrids,allGrids,0) =
                    interpolationWidth(three,allGrids,allGrids,number(k));
//            list[i].backupInterpolationWidth(three,allGrids,allGrids,0) =
//                    backupInterpolationWidth(three,allGrids,allGrids,number(k));
            list[i].interpolationOverlap(three,allGrids,allGrids,0) =
                    interpolationOverlap(three,allGrids,allGrids,number(k));
            list[i].maximumHoleCuttingDistance(nullRange,nullRange,allGrids)=
	      maximumHoleCuttingDistance(nullRange,nullRange,allGrids);
//             list[i].backupInterpolationOverlap(three,allGrids,allGrids,0) =
//                   backupInterpolationOverlap(three,allGrids,allGrids,number(k));
//             list[i].interpolationConditionLimit(allGrids,allGrids,0) =
//                     interpolationConditionLimit(allGrids,allGrids,number(k));
//             list[i].backupInterpolationConditionLimit(allGrids,allGrids,0) =
//                  backupInterpolationConditionLimit(allGrids,allGrids,number(k));
            list[i].interpolationPreference(allGrids,allGrids,0) =
                    interpolationPreference(allGrids,allGrids,number(k));
            list[i].mayInterpolate(allGrids,allGrids,0) =
                    mayInterpolate(allGrids,allGrids,number(k));
//             list[i].mayBackupInterpolate(allGrids,allGrids,0) =
//                     mayBackupInterpolate(allGrids,allGrids,number(k));
            // *wdh* 991023 list[i].mayCutHoles = mayCutHoles;
            list[i].mayCutHoles(allGrids,allGrids) = mayCutHoles(allGrids,allGrids);
            list[i].sharedSidesMayCutHoles(allGrids,allGrids) = sharedSidesMayCutHoles(allGrids,allGrids);
            list[i].multigridCoarseningRatio(three,allGrids,0) =
                    multigridCoarseningRatio(three,allGrids,number(k));
            list[i].multigridProlongationWidth(three,allGrids,0) =
                    multigridProlongationWidth(three,allGrids,number(k));
            list[i].multigridRestrictionWidth(three,allGrids,0) =
                    multigridRestrictionWidth(three,allGrids,number(k));
//            list[i].interpoleeGridRange(allGridsPlusOne,allGrids,0) =
//                    interpoleeGridRange(allGridsPlusOne,allGrids,number(k));
        } // end for
        for (i=0; i<numberOfCollections; i++) {
            const Integer des = ~(computedGeometry | what) &      (
              THEinterpolationCoordinates | THEinterpoleeGrid     |
              THEinterpoleeLocation       | THEinterpolationPoint |
              THEinterpolationCondition   | THEinverseMap         );
            if (des) list[i].destroy(des);
            const Integer upd = (computedGeometry | what) &       (
              THEinterpolationCoordinates | THEinterpoleeGrid     |
              THEinterpoleeLocation       | THEinterpolationPoint |
              THEinterpolationCondition   | THEinverseMap         );
            if (upd) list[i].update(upd, COMPUTEnothing);
        } // end for
        for (nG=0, k=0; k<numberOfGrids; k++) {
            const Integer j = nG(i = number(k))++;
            if ((computedGeometry | what) & THEinterpolationCoordinates)
              list[i]->interpolationCoordinates[j].reference(
                       interpolationCoordinates[k]);
            if ((computedGeometry | what) & THEinterpoleeGrid)
	    {
              list[i]->interpoleeGrid[j].reference(interpoleeGrid[k]);
              list[i]->variableInterpolationWidth[j].reference(variableInterpolationWidth[k]);
	    }
//             if ((computedGeometry | what) & THEinterpoleeGrid)
//               list[i]->interpoleePoint[j].reference(
//                        interpoleePoint[k]);
            if ((computedGeometry | what) & THEinterpoleeLocation)
              list[i]->interpoleeLocation[j].reference(
                       interpoleeLocation[k]);
            if ((computedGeometry | what) & THEinterpolationPoint)
              list[i]->interpolationPoint[j].reference(
                       interpolationPoint[k]);
//             if ((computedGeometry | what) & THEinterpolationCondition)
//               list[i]->interpolationCondition[j].reference(
//                        interpolationCondition[k]);
            if ((computedGeometry | what) & THEinverseMap) {
//                 list[i]->inverseCondition[j].reference(
//                          inverseCondition[k]);
                list[i]->inverseCoordinates[j].reference(
                         inverseCoordinates[k]);
                list[i]->inverseGrid[j].reference(
                         inverseGrid[k]);
//                list[i]->inverseCondition  .updateToMatchGrid(*list[i]);
                list[i]->inverseCoordinates.updateToMatchGrid(*list[i]);
                list[i]->inverseGrid       .updateToMatchGrid(*list[i]);

                for (Integer j2=0; j2<list[i].numberOfBaseGrids(); j2++) {
                    BoundaryAdjustmentArray
                      &bA12j = list[i]->boundaryAdjustment(j,j2),
                      &bA12k = boundaryAdjustment(k,j2);
                    bA12j.redim(bA12k);
                    if (bA12j.getNumberOfElements()) {
                        for (Integer kd=0; kd<numberOfDimensions; kd++)
                          for (Integer ks=0; ks<2; ks++) {
                            BoundaryAdjustment &bAj = bA12j(ks,kd), &bAk = bA12k(ks,kd);
                            bAj.reference(bAk);
			    
//                             bAj.computedGeometry = bAk.computedGeometry;
//                             bAj.boundaryAdjustment.reference(
//                             bAk.boundaryAdjustment);
//                             bAj.acrossGrid        .reference(
//                             bAk.acrossGrid);
//                             bAj.oppositeBoundary  .reference(
//                             bAk.oppositeBoundary);
                        } // end for, end for
                    } // end if
                } // end for


            } // end if
            list[i]->computedGeometry |= (computedGeometry | what) & (
              THEinterpolationCoordinates | THEinterpoleeGrid     |
              THEinterpoleeLocation       | THEinterpolationPoint |
              THEinterpolationCondition   | THEinverseMap         );
            list[i].updateReferences(computedGeometry | what);
        } // end for
    } // end if
    return what & THElists;
}

void 
CompositeGrid::
setHybridConnectivity(const int grid_,
		      intArray * gridIndex2UVertex_,
		      intArray & uVertex2GridIndex_,
		      intArray * gridVertex2UVertex_, // could be build from gridIndex2UVertex ?
		      intArray & boundaryFaceMapping_)
{

  rcData->hybridConnectivity.setCompositeGridHybridConnectivity(grid_,
								gridIndex2UVertex_, uVertex2GridIndex_,
								gridVertex2UVertex_, boundaryFaceMapping_);
  
}

const CompositeGridHybridConnectivity & 
CompositeGrid::
getHybridConnectivity() const
{
  return rcData->hybridConnectivity;
}

#undef COMPOSITE_GRID_FOR_3
