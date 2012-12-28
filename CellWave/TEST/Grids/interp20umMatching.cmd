create mappings 
  box 
    set corners 
    0 20 0 20 0. 10 
    lines 
    22 22 10 
    mappingName 
    upperGrid 
    boundary conditions 
    -1 -1 1 1 2 1 
    *   3 3 1 1 2 1 
    exit
  box 
    set corners 
    0 20 0 20 -10 0. 
    lines 
    22 22 10 
    mappingName 
    lowerGrid 
    boundary conditions 
    -1 -1 1 1 1 2 
    *   3 3 1 1 1 2 
    exit
  exit this menu 
generate an overlapping grid 
  upperGrid 
  lowerGrid 
  change parameters 
    non-cutting boundary points 
      upperGrid 
        front  (side=0,axis=2) 
      0 21 0 21 0 0 
      lowerGrid 
        back   (side=1,axis=2) 
      0 21 0 21 9 9 
      done 
      exit 
    compute overlap 
    exit 
  save a grid
  interp20umMatching.hdf
  flux boundary manual interpolation
  exit
