create mappings 
  rectangle 
    set corners 
    0 1 0 .5 
    lines 
*    11 5 
     22 10
*     44 20
    boundary conditions
    -1 -1 2 1
    mappingName 
    upperGrid 
    exit 
  rectangle 
    set corners 
    0 1 -.5 0 
    lines 
     22 10
    boundary conditions
    -1 -1 1 2
    mappingName 
    lowerGrid 
    exit 
  exit this menu 
* 
generate an overlapping grid 
  upperGrid 
  lowerGrid 
  change parameters 
    non-cutting boundary points 
      upperGrid 
        bottom (side=0,axis=1) 
       0 21 0 0 0 0
      lowerGrid 
        top    (side=1,axis=1) 
        0 21  9  9 0 0
      done 
      exit 
    compute overlap 
    exit
  save a grid
  interpMatching2.hdf
  flux boundary manual interpolation
  exit

