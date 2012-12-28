*
* simple 2D grid for 1um diameter X. Laevis oocyte 
*
create mappings
*
rectangle
  set corners
    -0.5 0.5 -0.5 0.5
**    -2. 2. -2. 2.
  lines
    32 32 
  boundary conditions
   0 0 0 0 
  mappingName
  square
exit
*
Annulus
  lines
    33 7
*  centre
*    0. 1.
    inner and outer radii
    0.25 0.5
  boundary conditions
    -1 -1 0 1 
exit
*
exit
generate an overlapping grid
    square
    Annulus
  done
  change parameters
    * choose implicit or explicit interpolation
    * interpolation type
    *   implicit for all grids
    ghost points
      all
      2 2 2 2 2 2
  exit
*  display intermediate results
  compute overlap
  exit
*
save an overlapping grid
oocyte_1um.hdf
oocyte
exit

