* Epithelial cell grid 
*
*  -- debugging the corner bug in Neumann bcs & Flux-Jump bcs
*
create mappings 
  ***---- TFICell 7: Top curve=14, bottom curve=13----- 
  *     .. top    curve:   26 
  *     .. bottom curve:   24 
  *   ---top curve 
  spline 
    enter spline points 
    16 
    9.947368         68.923077 
    9.800000         68.564103 
    9.652632         68.205128 
    9.505263         67.846154 
    9.357895         67.487179 
    9.210526         67.128205 
    9.063158         66.769231 
    8.915789         66.410256 
    8.768421         66.051282 
    8.621053         65.692308 
    8.473684         65.333333 
    8.326316         64.974359 
    8.178947         64.615385 
    8.031579         64.256410 
    7.884211         63.897436 
    7.736842         63.538462 
    index 
    lines 
    16 
    shape preserving (toggle) 
    mappingName 
    curve14 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    16 
    0.000000         75.384615 
    0.000000         74.307692 
    0.000000         73.230769 
    0.000000         72.153846 
    0.000000         71.076923 
    0.000000         70.000000 
    0.000000         68.923077 
    0.000000         67.846154 
    0.000000         66.769231 
    0.000000         65.692308 
    0.000000         64.615385 
    0.000000         63.538462 
    0.000000         62.461538 
    0.000000         61.384615 
    0.000000         60.307692 
    0.000000         59.230769 
    index 
    lines 
    16 
    shape preserving (toggle) 
    mappingName 
    curve13 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve13 
    choose top curve    (r_2=1) 
      curve14 
    lines 
    16 11 
    mappingName 
    cell7 
    boundary conditions 
    1 1 1 1 0 0 
    exit 
  * 
 exit this menu 
* 
generate an overlapping grid 
  cell7 
  done choosing mappings 
  compute overlap
  exit
save a grid
neumannCornerCell7.hdf
hirose 2d
exit
