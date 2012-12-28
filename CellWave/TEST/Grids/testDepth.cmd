create mappings
  smoothedPolygon
    vertices
    3
    0 0
    -5 10
    5 20
    curve or area (toggle)
    mappingName
    left edge
    exit
  smoothedPolygon
    vertices
    2
    20 5
    25 20
    curve or area (toggle)
    mappingName
    right edge
    exit
  tfi
    choose left curve   (r_1=0)
      done
      choose right curve  (r_1=1)
        right edge
      choose left curve   (r_1=0)
        left edge
      set view:0 0 0 0 1 0.999204 -0.00661184 -0.0393333 0.00188039 0.992877 -0.119132 0.0398408 0.118963 0.992099
      lines
      21 21
      exit
    elliptic
      transform which mapping?
      TFIMapping
      elliptic smoothing
        maximum number of iterations
        1
        smooth
        smooth
        smooth
        smooth
        boundary conditions
          left   (side=0,axis=0)
          slip orthogonal
          exit
        smooth
        reset elliptic transform
        boundary conditions
          left   (side=0,axis=0)
          dirichlet
          exit
        smooth
        test orthogonal
          exit
        smooth
        plot control function
          exit this menu
        exit
      exit
    depth mapping
      set view:0 0 0 0 1 0.989099 0.00249121 0.14723 0.114823 0.612919 -0.781759 -0.0921872 0.790142 0.605951
      constant depth
      -20
      reset:0
      set view:0 -0.461012 0.23803 0 2.56491 1 0 0 0 1 0 0 0 1
      reset:0
      set view:0 -0.017301 0.0273598 0 1 0.983113 0.0675918 -0.170061 -0.164379 0.734592 -0.658296 0.0804299 0.675134 0.733298
      exit

    view mappings
      depth-elliptic-TFIMapping
      change plot parameters of last plotted object
        exit
      exit
    change a mapping
    depth-elliptic-TFIMapping
      boundary conditions
      1 2 3 4 5 6
      plot
        colour boundaries by bc number
        exit
      reset:0
      set view:0 0 0 0 1 -0.205737 0.35701 0.911162 0.698983 -0.598032 0.392148 0.684905 0.717566 -0.126507
      reset:0
      set view:0 0 0 0 1 0.995344 0.0211661 0.0940312 0.0446533 0.763312 -0.644485 -0.0854164 0.645683 0.758813
      exit
 
   exit this menu
  exit
