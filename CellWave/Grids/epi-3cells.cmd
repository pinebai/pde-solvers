create mappings
  smoothedPolygon
    vertices
    3
    0 0
    -50 150
    20 20
    curve or area (toggle)
    vertices
    3
    0 0 
    -50 150
    20 300
    vertices
    3
    0 0
    -5 10
    5 20
    lines
    mappingName
    left edge
    exit
  smoothedPolygon
    curve or area (toggle)
    vertices
    2
    20 -5
    25 25
    lines
    21
    mappingName
    rightBoundary
    exit
  tfi
    choose left curve   (r_1=0)
      left edge
    choose right curve  (r_1=1)
      rightBoundary
    exit
  elliptic
    elliptic smoothing
      generate grid
      generate grid
      generate grid
      lines
      exit
    lines
    21 21
    elliptic smoothing
      boundary conditions
        left   (side=0,axis=0)
        slip orthogonal
        exit
      generate grid
      reset elliptic transform
      generate grid
      reset elliptic transform
      boundary conditions
        left   (side=0,axis=0)
        dirichlet
        exit
      generate grid
      exit
    lines
    21
    lines
    21 21
    elliptic smoothing
      generate grid
      test orthogonal
        plot shaded surfaces (3D) toggle
        plot grid lines on boundaries (3D) toggle
        plot non-physical boundaries (toggle)
        plot mapping edges (toggle)
        plot
        plot
        erase
        exit
      generate grid
      plot residual
        exit this menu
      exit
    plot
      clear all:0
      plot
      exit
    exit
  view mappings
    TFIMapping
    exit
  change a mapping
  left edge
    corners
    do not specify positions of the corners
    t-stretch
    0 1
    1.5 1
    0 1
    exit
  view mappings
    TFIMapping
    exit
  change a mapping
  TFIMapping
    choose left curve   (r_1=0)
      left edge
    exit
  change a mapping
  TFIMapping
    exit
  change a mapping
  elliptic-TFIMapping
    transform which mapping?
    TFIMapping
    elliptic smoothing
      generate grid
      exit
    exit
  view mappings
    TFIMapping
    exit
  change a mapping
  TFIMapping
    lines
    21 21
    exit
  elliptic
    transform which mapping?
    TFIMapping
    elliptic smoothing
      generate grid
      exit
    exit
  view mappings
    TFIMapping
    exit
  smoothedPolygon
    vertices
    2
    -20 10
    -25 20
    curve or area (toggle)
    lines
    11
    mappingName
    leftTopEdge
    exit
  smoothedPolygon
    vertices
    3
    -25 -5
    -25 0
    -20 10
    curve or area (toggle)
    lines
    11
    mappingName
    leftBottomEdge
    exit
  reparameterize
    transform which mapping?
    left edge
    restrict parameter space
      set corners
      0.5 1
      exit
    mappingName
    restrictedTop
    exit
  reparameterize
    transform which mapping?
    left edge
    restrict parameter space
      set corners
      0 0.5
      exit
    set corners
    0 0.4
    set corners
    0 0.45
    mappingName
    restrictionBottom
    exit
  change a mapping
  restrictedTop
    restrict parameter space
      set corners
      0.45 1
      exit
    exit
  tfi
    choose left curve   (r_1=0)
      leftTopEdge
    choose right curve  (r_1=1)
      restrictedTop
    lines
    mappingName
    leftTop
    exit
  tfi
    choose left curve   (r_1=0)
      leftBottomEdge
    choose right curve  (r_1=1)
      restrictionBottom
    mappingName
    leftBottom
    exit
  view mappings
    TFIMapping
    leftTop
    leftBottom
    exit
  change a mapping
  TFIMapping
    exit
  change a mapping
  leftTop
    lines
    21 21
    exit
  change a mapping
  leftBottom
    lines
    21 21
    exit
  view mappings
    TFIMapping
    leftTop
    leftBottom
    bigger:0
    bigger:0
    hardcopy file name:0 /home/pfast/epigrid1.ps
    hardcopy save:0
    hardcopy close dialog:0
    exit
  exit this menu
generate an overlapping grid
  TFIMapping
  leftTop
  leftBottom
  done choosing mappings
  compute overlap
  help
  exit
exit
