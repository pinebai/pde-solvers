create mappings
  box
    set corners
    0 1 0 1 0.1 0.5
    lines
    22 22 10
    mappingName
    upperGrid
    boundary conditions
    -1 -1 1 1 1 2
    plot
      colour boundaries by bc number
      erase
      exit
    boundary conditions
    -1 -1 1 1 2 1
    exit
  box
    set corners
    0 1 0 1 -0.5 -0.1
    lines
    22 22 10
    boundary conditions
    -1 -1 1 1 1 2
    plot
      colour boundaries by bc number
      exit
    mappingName
    lowerGrid
    exit
  exit this menu
generate an overlapping grid
  upperGrid
  lowerGrid
  change parameters
    non-cutting boundary points
      upperGrid
        back   (side=1,axis=2)
      0 21 0 21 9 9
      lowerGrid
        front  (side=0,axis=2)
      0 21 0 21 0 0
      done
      exit

    compute overlap
    exit

  exit
