create mappings
  annulus
    start and end angles
    0 0.3
    lines
    21 13
    set view:0 0 0 0 1 0.999527 0.0106735 0.028858 -0.00413014 0.975958 -0.217918 -0.0304902 0.217696 0.97554
    exit
  sweep
    choose reference mapping
    Annulus
    extrude
    0. 1.
    boundary conditions
    2 1 1 1 1 1
    plot
      exit
    mappingName
    Sweep1
    exit
  annulus
    reset:0
    start and end angles
    0.5 0.8
    lines
    30 18
    centre for annulus
    1.5 0.
    boundary conditions
    2 1 1 1
    mappingName
    Annulus2
    exit
  sweep
    choose reference mapping
    Annulus2
    extrude
    0. 1.
    lines
    30 18 20
    mappingName
    Sweep2
    boundary conditions
    2 1 1 1 2 1
    exit
  box
    set corners
    0.5 1.5 -1. -0.5 -1 0
    lines
    30 15 30
    boundary conditions
    1 1 1 1 1 2
    mappingName
    Box
    exit
  view mappings
    Sweep1
    Sweep2
    Box
    exit
  exit this menu
generate an overlapping grid
  Sweep1
  Sweep2
  Box
  done choosing mappings
  change parameters
    prevent hole cutting
      all
      all
      done
    interpolation type
      no change
      done
    prevent interpolation
      all
      all
      done
    exit
  change parameters
    ghost points
      all
      2 2 2 2 2 2
    exit
  compute overlap
  exit
save a grid
sweepSweepBox.hdf
SweepSweepBox
exit
