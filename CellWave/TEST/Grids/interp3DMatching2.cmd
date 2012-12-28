create mappings 
  box 
    set corners 
    0 1 0 1 0. 0.5 
    lines 
    22 22 10 
    mappingName 
    upperGrid 
    boundary conditions 
    -1 -1 1 1 2 1 
    *   3 3 1 1 2 1 
    plot 
      colour boundaries by bc number 
      erase 
      exit 
    exit 
  box 
    set corners 
    0 1 0 1 -0.5 0. 
    lines 
    22 22 10 
    mappingName 
    lowerGrid 
    boundary conditions 
    -1 -1 1 1 1 2 
    *   3 3 1 1 1 2 
    plot 
      colour boundaries by bc number 
      exit 
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
  interp3DMatching2.hdf
  flux boundary manual interpolation
  exit
