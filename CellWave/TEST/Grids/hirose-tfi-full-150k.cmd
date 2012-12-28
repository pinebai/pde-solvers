* Epithelial cell grid 
*
*   dx, dy = 1um
*   minimum nx, ny = 5, 6
*   .. some cells are still somewhat skewed.
*
**total    points: discretization=  150576
**
**------------------------------------------------------------------------
**                                                     cpu     percentage
**update geometry.....................................4.28e+00    78.07
**interpolate boundaries..............................7.74e-03     0.14
**cut holes...........................................0.00e+00     0.00
**remove exterior points..............................6.16e-02     1.13
**improper interpolation..............................6.04e-01    11.02
**proper interpolation................................1.83e-03     0.03
**all interpolation...................................1.85e-01     3.37
**improve quality of interpolation....................0.00e+00     0.00
**remove redundant points.............................1.34e-01     2.45
**sum of above........................................5.27e+00    96.22
**total...............................................5.48e+00   100.00
**------------------------------------------------------------------------
*****Size of composite grid=22798.56 Kbytes
*
*
create mappings 
  ***---- TFICell 1: Top curve=2, bottom curve=1----- 
  *     .. top    curve:    5 
  *     .. bottom curve:    2    3 
  *   ---top curve 
  spline 
    enter spline points 
    21 
    0.000000         84.000000 
    1.050000         84.000000 
    2.100000         84.000000 
    3.150000         84.000000 
    4.200000         84.000000 
    5.250000         84.000000 
    6.300000         84.000000 
    7.350000         84.000000 
    8.400000         84.000000 
    9.450000         84.000000 
    10.500000         84.000000 
    11.550000         84.000000 
    12.600000         84.000000 
    13.650000         84.000000 
    14.700000         84.000000 
    15.750000         84.000000 
    16.800000         84.000000 
    17.850000         84.000000 
    18.900000         84.000000 
    19.950000         84.000000 
    21.000000         84.000000 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve2 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    21 
    0.000000         75.384615 
    0.828947         74.846154 
    1.657895         74.307692 
    2.486842         73.769231 
    3.315789         73.230769 
    4.144737         72.692308 
    4.973684         72.153846 
    5.802632         71.615385 
    6.631579         71.076923 
    7.460526         70.538462 
    8.289474         70.000000 
    9.118421         69.461538 
    9.947368         68.923077 
    10.914474         68.923077 
    11.881579         68.923077 
    12.848684         68.923077 
    13.815789         68.923077 
    14.782895         68.923077 
    15.750000         68.923077 
    16.717105         68.923077 
    17.684211         68.923077 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve1 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve1 
    choose top curve    (r_2=1) 
      curve2 
    lines 
    21 15 
    mappingName 
    celltop1 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop1 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 1 1 35 
    lines 
    21 15 16 
    mappingName 
    cell1 
    exit 
  ***---- TFICell 2: Top curve=4, bottom curve=3----- 
  *     .. top    curve:   10    9    8 
  *     .. bottom curve:    4    6 
  *   ---top curve 
  spline 
    enter spline points 
    27 
    29.842105         84.000000 
    30.671053         83.192308 
    31.500000         82.384615 
    32.328947         81.576923 
    33.157895         80.769231 
    33.378947         79.835897 
    33.600000         78.902564 
    33.821053         77.969231 
    34.042105         77.035897 
    34.263158         76.102564 
    34.484211         75.169231 
    34.705263         74.235897 
    34.926316         73.302564 
    35.147368         72.369231 
    35.368421         71.435897 
    35.589474         70.502564 
    35.810526         69.569231 
    36.031579         68.635897 
    36.252632         67.702564 
    36.473684         66.769231 
    36.315789         65.846154 
    36.157895         64.923077 
    36.000000         64.000000 
    35.842105         63.076923 
    35.684211         62.153846 
    35.526316         61.230769 
    35.368421         60.307692 
    index 
    lines 
    27 
    shape preserving (toggle) 
    mappingName 
    curve4 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    27 
    21.000000         84.000000 
    20.763158         82.923077 
    20.526316         81.846154 
    20.289474         80.769231 
    20.052632         79.692308 
    19.815789         78.615385 
    19.578947         77.538462 
    19.342105         76.461538 
    19.105263         75.384615 
    18.868421         74.307692 
    18.631579         73.230769 
    18.394737         72.153846 
    18.157895         71.076923 
    17.921053         70.000000 
    17.684211         68.923077 
    18.421053         68.205128 
    19.157895         67.487179 
    19.894737         66.769231 
    20.631579         66.051282 
    21.368421         65.333333 
    22.105263         64.615385 
    22.842105         63.897436 
    23.578947         63.179487 
    24.315789         62.461538 
    25.052632         61.743590 
    25.789474         61.025641 
    26.526316         60.307692 
    index 
    lines 
    27 
    shape preserving (toggle) 
    mappingName 
    curve3 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve3 
    choose top curve    (r_2=1) 
      curve4 
    lines 
    27 8 
    mappingName 
    celltop2 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop2 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 2 2 35 
    lines 
    27 8 16 
    mappingName 
    cell2 
    exit 
  ***---- TFICell 3: Top curve=6, bottom curve=5----- 
  *     .. top    curve:   14 
  *     .. bottom curve:   12 
  *   ---top curve 
  spline 
    enter spline points 
    29 
    29.842105         84.000000 
    30.907895         84.000000 
    31.973684         84.000000 
    33.039474         84.000000 
    34.105263         84.000000 
    35.171053         84.000000 
    36.236842         84.000000 
    37.302632         84.000000 
    38.368421         84.000000 
    39.434211         84.000000 
    40.500000         84.000000 
    41.565789         84.000000 
    42.631579         84.000000 
    43.697368         84.000000 
    44.763158         84.000000 
    45.828947         84.000000 
    46.894737         84.000000 
    47.960526         84.000000 
    49.026316         84.000000 
    50.092105         84.000000 
    51.157895         84.000000 
    52.223684         84.000000 
    53.289474         84.000000 
    54.355263         84.000000 
    55.421053         84.000000 
    56.486842         84.000000 
    57.552632         84.000000 
    58.618421         84.000000 
    59.684211         84.000000 
    index 
    lines 
    29 
    shape preserving (toggle) 
    mappingName 
    curve6 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    29 
    33.157895         80.769231 
    33.986842         80.769231 
    34.815789         80.769231 
    35.644737         80.769231 
    36.473684         80.769231 
    37.302632         80.769231 
    38.131579         80.769231 
    38.960526         80.769231 
    39.789474         80.769231 
    40.618421         80.769231 
    41.447368         80.769231 
    42.276316         80.769231 
    43.105263         80.769231 
    43.934211         80.769231 
    44.763158         80.769231 
    45.592105         80.769231 
    46.421053         80.769231 
    47.250000         80.769231 
    48.078947         80.769231 
    48.907895         80.769231 
    49.736842         80.769231 
    50.565789         80.769231 
    51.394737         80.769231 
    52.223684         80.769231 
    53.052632         80.769231 
    53.881579         80.769231 
    54.710526         80.769231 
    55.539474         80.769231 
    56.368421         80.769231 
    index 
    lines 
    29 
    shape preserving (toggle) 
    mappingName 
    curve5 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve5 
    choose top curve    (r_2=1) 
      curve6 
    lines 
    29 6 
    mappingName 
    celltop3 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop3 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 1 3 35 
    lines 
    29 6 16 
    mappingName 
    cell3 
    exit 
  ***---- TFICell 4: Top curve=8, bottom curve=7----- 
  *     .. top    curve:    9 
  *     .. bottom curve:   46   15 
  *   ---top curve 
  spline 
    enter spline points 
    17 
    36.473684         66.769231 
    36.266447         67.644231 
    36.059211         68.519231 
    35.851974         69.394231 
    35.644737         70.269231 
    35.437500         71.144231 
    35.230263         72.019231 
    35.023026         72.894231 
    34.815789         73.769231 
    34.608553         74.644231 
    34.401316         75.519231 
    34.194079         76.394231 
    33.986842         77.269231 
    33.779605         78.144231 
    33.572368         79.019231 
    33.365132         79.894231 
    33.157895         80.769231 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve8 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    17 
    47.526316         66.769231 
    48.330144         67.552448 
    49.133971         68.335664 
    49.937799         69.118881 
    50.741627         69.902098 
    51.545455         70.685315 
    52.349282         71.468531 
    53.153110         72.251748 
    53.956938         73.034965 
    54.760766         73.818182 
    55.564593         74.601399 
    56.368421         75.384615 
    56.368421         76.461538 
    56.368421         77.538462 
    56.368421         78.615385 
    56.368421         79.692308 
    56.368421         80.769231 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve7 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve7 
    choose top curve    (r_2=1) 
      curve8 
    lines 
    17 23 
    mappingName 
    celltop4 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop4 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 4 35 
    lines 
    17 23 16 
    mappingName 
    cell4 
    exit 
  ***---- TFICell 5: Top curve=10, bottom curve=9----- 
  *     .. top    curve:   19   18 
  *     .. bottom curve:   13   15   16 
  *   ---top curve 
  spline 
    enter spline points 
    27 
    75.157895         84.000000 
    75.274238         83.206478 
    75.390582         82.412955 
    75.506925         81.619433 
    75.623269         80.825911 
    75.739612         80.032389 
    75.855956         79.238866 
    75.972299         78.445344 
    76.088643         77.651822 
    76.204986         76.858300 
    76.321330         76.064777 
    76.437673         75.271255 
    76.554017         74.477733 
    76.670360         73.684211 
    76.786704         72.890688 
    76.903047         72.097166 
    77.019391         71.303644 
    77.135734         70.510121 
    77.252078         69.716599 
    77.368421         68.923077 
    77.052632         68.153846 
    76.736842         67.384615 
    76.421053         66.615385 
    76.105263         65.846154 
    75.789474         65.076923 
    75.473684         64.307692 
    75.157895         63.538462 
    index 
    lines 
    27 
    shape preserving (toggle) 
    mappingName 
    curve10 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    27 
    59.684211         84.000000 
    58.578947         82.923077 
    57.473684         81.846154 
    56.368421         80.769231 
    56.368421         79.692308 
    56.368421         78.615385 
    56.368421         77.538462 
    56.368421         76.461538 
    56.368421         75.384615 
    57.105263         74.726496 
    57.842105         74.068376 
    58.578947         73.410256 
    59.315789         72.752137 
    60.052632         72.094017 
    60.789474         71.435897 
    61.526316         70.777778 
    62.263158         70.119658 
    63.000000         69.461538 
    63.736842         68.803419 
    64.473684         68.145299 
    65.210526         67.487179 
    65.947368         66.829060 
    66.684211         66.170940 
    67.421053         65.512821 
    68.157895         64.854701 
    68.894737         64.196581 
    69.631579         63.538462 
    index 
    lines 
    27 
    shape preserving (toggle) 
    mappingName 
    curve9 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve9 
    choose top curve    (r_2=1) 
      curve10 
    lines 
    27 15 
    mappingName 
    celltop5 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop5 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 2 5 35 
    lines 
    27 15 16 
    mappingName 
    cell5 
    exit 
  ***---- TFICell 6: Top curve=12, bottom curve=11----- 
  *     .. top    curve:   22 
  *     .. bottom curve:   19 
  *   ---top curve 
  spline 
    enter spline points 
    15 
    84.000000         84.000000 
    84.000000         83.230769 
    84.000000         82.461538 
    84.000000         81.692308 
    84.000000         80.923077 
    84.000000         80.153846 
    84.000000         79.384615 
    84.000000         78.615385 
    84.000000         77.846154 
    84.000000         77.076923 
    84.000000         76.307692 
    84.000000         75.538462 
    84.000000         74.769231 
    84.000000         74.000000 
    84.000000         73.230769 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve12 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    15 
    75.157895         84.000000 
    75.315789         82.923077 
    75.473684         81.846154 
    75.631579         80.769231 
    75.789474         79.692308 
    75.947368         78.615385 
    76.105263         77.538462 
    76.263158         76.461538 
    76.421053         75.384615 
    76.578947         74.307692 
    76.736842         73.230769 
    76.894737         72.153846 
    77.052632         71.076923 
    77.210526         70.000000 
    77.368421         68.923077 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve11 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve11 
    choose top curve    (r_2=1) 
      curve12 
    lines 
    15 8 
    mappingName 
    celltop6 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop6 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 1 6 35 
    lines 
    15 8 16 
    mappingName 
    cell6 
    exit 
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
    celltop7 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop7 
    constant depth 
    -9.080000 
    boundary conditions 
    2 1 2 2 7 35 
    lines 
    16 11 16 
    mappingName 
    cell7 
    exit 
  ***---- TFICell 8: Top curve=16, bottom curve=15----- 
  *     .. top    curve:    6   29 
  *     .. bottom curve:   26   27 
  *   ---top curve 
  spline 
    enter spline points 
    21 
    17.684211         68.923077 
    18.488038         68.139860 
    19.291866         67.356643 
    20.095694         66.573427 
    20.899522         65.790210 
    21.703349         65.006993 
    22.507177         64.223776 
    23.311005         63.440559 
    24.114833         62.657343 
    24.918660         61.874126 
    25.722488         61.090909 
    26.526316         60.307692 
    26.280702         59.350427 
    26.035088         58.393162 
    25.789474         57.435897 
    25.543860         56.478632 
    25.298246         55.521368 
    25.052632         54.564103 
    24.807018         53.606838 
    24.561404         52.649573 
    24.315789         51.692308 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve16 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    21 
    9.947368         68.923077 
    9.505263         67.846154 
    9.063158         66.769231 
    8.621053         65.692308 
    8.178947         64.615385 
    7.736842         63.538462 
    8.105263         62.748718 
    8.473684         61.958974 
    8.842105         61.169231 
    9.210526         60.379487 
    9.578947         59.589744 
    9.947368         58.800000 
    10.315789         58.010256 
    10.684211         57.220513 
    11.052632         56.430769 
    11.421053         55.641026 
    11.789474         54.851282 
    12.157895         54.061538 
    12.526316         53.271795 
    12.894737         52.482051 
    13.263158         51.692308 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve15 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve15 
    choose top curve    (r_2=1) 
      curve16 
    lines 
    21 11 
    mappingName 
    celltop8 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop8 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 8 35 
    lines 
    21 11 16 
    mappingName 
    cell8 
    exit 
  ***---- TFICell 9: Top curve=18, bottom curve=17----- 
  *     .. top    curve:   41 
  *     .. bottom curve:   29   39 
  *   ---top curve 
  spline 
    enter spline points 
    20 
    35.368421         60.307692 
    35.717452         59.854251 
    36.066482         59.400810 
    36.415512         58.947368 
    36.764543         58.493927 
    37.113573         58.040486 
    37.462604         57.587045 
    37.811634         57.133603 
    38.160665         56.680162 
    38.509695         56.226721 
    38.858726         55.773279 
    39.207756         55.319838 
    39.556787         54.866397 
    39.905817         54.412955 
    40.254848         53.959514 
    40.603878         53.506073 
    40.952909         53.052632 
    41.301939         52.599190 
    41.650970         52.145749 
    42.000000         51.692308 
    index 
    lines 
    20 
    shape preserving (toggle) 
    mappingName 
    curve18 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    20 
    26.526316         60.307692 
    26.250000         59.230769 
    25.973684         58.153846 
    25.697368         57.076923 
    25.421053         56.000000 
    25.144737         54.923077 
    24.868421         53.846154 
    24.592105         52.769231 
    24.315789         51.692308 
    24.918660         50.811189 
    25.521531         49.930070 
    26.124402         49.048951 
    26.727273         48.167832 
    27.330144         47.286713 
    27.933014         46.405594 
    28.535885         45.524476 
    29.138756         44.643357 
    29.741627         43.762238 
    30.344498         42.881119 
    30.947368         42.000000 
    index 
    lines 
    20 
    shape preserving (toggle) 
    mappingName 
    curve17 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve17 
    choose top curve    (r_2=1) 
      curve18 
    lines 
    20 14 
    mappingName 
    celltop9 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop9 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 9 35 
    lines 
    20 14 16 
    mappingName 
    cell9 
    exit 
  ***---- TFICell 10: Top curve=20, bottom curve=19----- 
  *     .. top    curve:   44   43 
  *     .. bottom curve:    8   41 
  *   ---top curve 
  spline 
    enter spline points 
    18 
    47.526316         66.769231 
    48.002429         65.775148 
    48.478543         64.781065 
    48.954656         63.786982 
    49.430769         62.792899 
    49.906883         61.798817 
    50.382996         60.804734 
    50.859109         59.810651 
    51.335223         58.816568 
    51.811336         57.822485 
    52.287449         56.828402 
    52.763563         55.834320 
    53.239676         54.840237 
    53.715789         53.846154 
    53.550000         52.769231 
    53.384211         51.692308 
    53.218421         50.615385 
    53.052632         49.538462 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve20 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    18 
    36.473684         66.769231 
    36.289474         65.692308 
    36.105263         64.615385 
    35.921053         63.538462 
    35.736842         62.461538 
    35.552632         61.384615 
    35.368421         60.307692 
    35.971292         59.524476 
    36.574163         58.741259 
    37.177033         57.958042 
    37.779904         57.174825 
    38.382775         56.391608 
    38.985646         55.608392 
    39.588517         54.825175 
    40.191388         54.041958 
    40.794258         53.258741 
    41.397129         52.475524 
    42.000000         51.692308 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve19 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve19 
    choose top curve    (r_2=1) 
      curve20 
    lines 
    18 11 
    mappingName 
    celltop10 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop10 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 10 35 
    lines 
    18 11 16 
    mappingName 
    cell10 
    exit 
  ***---- TFICell 11: Top curve=22, bottom curve=21----- 
  *     .. top    curve:   16 
  *     .. bottom curve:   44 
  *   ---top curve 
  spline 
    enter spline points 
    17 
    56.368421         75.384615 
    57.197368         74.644231 
    58.026316         73.903846 
    58.855263         73.163462 
    59.684211         72.423077 
    60.513158         71.682692 
    61.342105         70.942308 
    62.171053         70.201923 
    63.000000         69.461538 
    63.828947         68.721154 
    64.657895         67.980769 
    65.486842         67.240385 
    66.315789         66.500000 
    67.144737         65.759615 
    67.973684         65.019231 
    68.802632         64.278846 
    69.631579         63.538462 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve22 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    17 
    47.526316         66.769231 
    47.913158         65.961538 
    48.300000         65.153846 
    48.686842         64.346154 
    49.073684         63.538462 
    49.460526         62.730769 
    49.847368         61.923077 
    50.234211         61.115385 
    50.621053         60.307692 
    51.007895         59.500000 
    51.394737         58.692308 
    51.781579         57.884615 
    52.168421         57.076923 
    52.555263         56.269231 
    52.942105         55.461538 
    53.328947         54.653846 
    53.715789         53.846154 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve21 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve21 
    choose top curve    (r_2=1) 
      curve22 
    lines 
    17 18 
    mappingName 
    celltop11 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop11 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 11 35 
    lines 
    17 18 16 
    mappingName 
    cell11 
    exit 
  ***---- TFICell 12: Top curve=24, bottom curve=23----- 
  *     .. top    curve:   50 
  *     .. bottom curve:   18   48 
  *   ---top curve 
  spline 
    enter spline points 
    15 
    84.000000         73.230769 
    84.000000         72.153846 
    84.000000         71.076923 
    84.000000         70.000000 
    84.000000         68.923077 
    84.000000         67.846154 
    84.000000         66.769231 
    84.000000         65.692308 
    84.000000         64.615385 
    84.000000         63.538462 
    84.000000         62.461538 
    84.000000         61.384615 
    84.000000         60.307692 
    84.000000         59.230769 
    84.000000         58.153846 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve24 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    15 
    77.368421         68.923077 
    77.052632         68.153846 
    76.736842         67.384615 
    76.421053         66.615385 
    76.105263         65.846154 
    75.789474         65.076923 
    75.473684         64.307692 
    75.157895         63.538462 
    75.473684         62.923077 
    75.789474         62.307692 
    76.105263         61.692308 
    76.421053         61.076923 
    76.736842         60.461538 
    77.052632         59.846154 
    77.368421         59.230769 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve23 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve23 
    choose top curve    (r_2=1) 
      curve24 
    lines 
    15 7 
    mappingName 
    celltop12 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop12 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 1 12 35 
    lines 
    15 7 16 
    mappingName 
    cell12 
    exit 
  ***---- TFICell 13: Top curve=26, bottom curve=25----- 
  *     .. top    curve:   32   31 
  *     .. bottom curve:   25 
  *   ---top curve 
  spline 
    enter spline points 
    16 
    13.263158         51.692308 
    12.710526         50.615385 
    12.157895         49.538462 
    11.605263         48.461538 
    11.052632         47.384615 
    10.047847         47.188811 
    9.043062         46.993007 
    8.038278         46.797203 
    7.033493         46.601399 
    6.028708         46.405594 
    5.023923         46.209790 
    4.019139         46.013986 
    3.014354         45.818182 
    2.009569         45.622378 
    1.004785         45.426573 
    0.000000         45.230769 
    index 
    lines 
    16 
    shape preserving (toggle) 
    mappingName 
    curve26 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    16 
    7.736842         63.538462 
    7.221053         63.251282 
    6.705263         62.964103 
    6.189474         62.676923 
    5.673684         62.389744 
    5.157895         62.102564 
    4.642105         61.815385 
    4.126316         61.528205 
    3.610526         61.241026 
    3.094737         60.953846 
    2.578947         60.666667 
    2.063158         60.379487 
    1.547368         60.092308 
    1.031579         59.805128 
    0.515789         59.517949 
    -0.000000         59.230769 
    index 
    lines 
    16 
    shape preserving (toggle) 
    mappingName 
    curve25 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve25 
    choose top curve    (r_2=1) 
      curve26 
    lines 
    16 14 
    mappingName 
    celltop13 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop13 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 1 2 13 35 
    lines 
    16 14 16 
    mappingName 
    cell13 
    exit 
  ***---- TFICell 14: Top curve=28, bottom curve=27----- 
  *     .. top    curve:   32   28 
  *     .. bottom curve:   37   38 
  *   ---top curve 
  spline 
    enter spline points 
    15 
    11.052632         47.384615 
    11.605263         48.461538 
    12.157895         49.538462 
    12.710526         50.615385 
    13.263158         51.692308 
    14.368421         51.692308 
    15.473684         51.692308 
    16.578947         51.692308 
    17.684211         51.692308 
    18.789474         51.692308 
    19.894737         51.692308 
    21.000000         51.692308 
    22.105263         51.692308 
    23.210526         51.692308 
    24.315789         51.692308 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve28 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    15 
    17.684211         35.538462 
    19.065789         35.538462 
    20.447368         35.538462 
    21.828947         35.538462 
    23.210526         35.538462 
    23.984211         36.184615 
    24.757895         36.830769 
    25.531579         37.476923 
    26.305263         38.123077 
    27.078947         38.769231 
    27.852632         39.415385 
    28.626316         40.061538 
    29.400000         40.707692 
    30.173684         41.353846 
    30.947368         42.000000 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve27 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve27 
    choose top curve    (r_2=1) 
      curve28 
    lines 
    15 13 
    mappingName 
    celltop14 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop14 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 14 35 
    lines 
    15 13 16 
    mappingName 
    cell14 
    exit 
  ***---- TFICell 15: Top curve=30, bottom curve=29----- 
  *     .. top    curve:   57 
  *     .. bottom curve:   40   55 
  *   ---top curve 
  spline 
    enter spline points 
    24 
    53.052632         49.538462 
    53.004577         48.882943 
    52.956522         48.227425 
    52.908467         47.571906 
    52.860412         46.916388 
    52.812357         46.260870 
    52.764302         45.605351 
    52.716247         44.949833 
    52.668192         44.294314 
    52.620137         43.638796 
    52.572082         42.983278 
    52.524027         42.327759 
    52.475973         41.672241 
    52.427918         41.016722 
    52.379863         40.361204 
    52.331808         39.705686 
    52.283753         39.050167 
    52.235698         38.394649 
    52.187643         37.739130 
    52.139588         37.083612 
    52.091533         36.428094 
    52.043478         35.772575 
    51.995423         35.117057 
    51.947368         34.461538 
    index 
    lines 
    24 
    shape preserving (toggle) 
    mappingName 
    curve30 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    24 
    42.000000         51.692308 
    41.149798         50.946746 
    40.299595         50.201183 
    39.449393         49.455621 
    38.599190         48.710059 
    37.748988         47.964497 
    36.898785         47.218935 
    36.048583         46.473373 
    35.198381         45.727811 
    34.348178         44.982249 
    33.497976         44.236686 
    32.647773         43.491124 
    31.797571         42.745562 
    30.947368         42.000000 
    31.610526         41.246154 
    32.273684         40.492308 
    32.936842         39.738462 
    33.600000         38.984615 
    34.263158         38.230769 
    34.926316         37.476923 
    35.589474         36.723077 
    36.252632         35.969231 
    36.915789         35.215385 
    37.578947         34.461538 
    index 
    lines 
    24 
    shape preserving (toggle) 
    mappingName 
    curve29 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve29 
    choose top curve    (r_2=1) 
      curve30 
    lines 
    24 14 
    mappingName 
    celltop15 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop15 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 15 35 
    lines 
    24 14 16 
    mappingName 
    cell15 
    exit 
  ***---- TFICell 16: Top curve=32, bottom curve=31----- 
  *     .. top    curve:   52   51 
  *     .. bottom curve:   17   47 
  *   ---top curve 
  spline 
    enter spline points 
    32 
    77.368421         59.230769 
    76.858300         58.236686 
    76.348178         57.242604 
    75.838057         56.248521 
    75.327935         55.254438 
    74.817814         54.260355 
    74.307692         53.266272 
    73.797571         52.272189 
    73.287449         51.278107 
    72.777328         50.284024 
    72.267206         49.289941 
    71.757085         48.295858 
    71.246964         47.301775 
    70.736842         46.307692 
    69.754386         46.487179 
    68.771930         46.666667 
    67.789474         46.846154 
    66.807018         47.025641 
    65.824561         47.205128 
    64.842105         47.384615 
    63.859649         47.564103 
    62.877193         47.743590 
    61.894737         47.923077 
    60.912281         48.102564 
    59.929825         48.282051 
    58.947368         48.461538 
    57.964912         48.641026 
    56.982456         48.820513 
    56.000000         49.000000 
    55.017544         49.179487 
    54.035088         49.358974 
    53.052632         49.538462 
    index 
    lines 
    32 
    shape preserving (toggle) 
    mappingName 
    curve32 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    32 
    75.157895         63.538462 
    74.236842         63.538462 
    73.315789         63.538462 
    72.394737         63.538462 
    71.473684         63.538462 
    70.552632         63.538462 
    69.631579         63.538462 
    68.994947         63.150769 
    68.358316         62.763077 
    67.721684         62.375385 
    67.085053         61.987692 
    66.448421         61.600000 
    65.811789         61.212308 
    65.175158         60.824615 
    64.538526         60.436923 
    63.901895         60.049231 
    63.265263         59.661538 
    62.628632         59.273846 
    61.992000         58.886154 
    61.355368         58.498462 
    60.718737         58.110769 
    60.082105         57.723077 
    59.445474         57.335385 
    58.808842         56.947692 
    58.172211         56.560000 
    57.535579         56.172308 
    56.898947         55.784615 
    56.262316         55.396923 
    55.625684         55.009231 
    54.989053         54.621538 
    54.352421         54.233846 
    53.715789         53.846154 
    index 
    lines 
    32 
    shape preserving (toggle) 
    mappingName 
    curve31 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve31 
    choose top curve    (r_2=1) 
      curve32 
    lines 
    32 6 
    mappingName 
    celltop16 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop16 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 16 35 
    lines 
    32 6 16 
    mappingName 
    cell16 
    exit 
  ***---- TFICell 17: Top curve=34, bottom curve=33----- 
  *     .. top    curve:   54 
  *     .. bottom curve:   52 
  *   ---top curve 
  spline 
    enter spline points 
    17 
    84.000000         58.153846 
    84.000000         57.076923 
    84.000000         56.000000 
    84.000000         54.923077 
    84.000000         53.846154 
    84.000000         52.769231 
    84.000000         51.692308 
    84.000000         50.615385 
    84.000000         49.538462 
    84.000000         48.461538 
    84.000000         47.384615 
    84.000000         46.307692 
    84.000000         45.230769 
    84.000000         44.153846 
    84.000000         43.076923 
    84.000000         42.000000 
    84.000000         40.923077 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve34 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    17 
    77.368421         59.230769 
    76.953947         58.423077 
    76.539474         57.615385 
    76.125000         56.807692 
    75.710526         56.000000 
    75.296053         55.192308 
    74.881579         54.384615 
    74.467105         53.576923 
    74.052632         52.769231 
    73.638158         51.961538 
    73.223684         51.153846 
    72.809211         50.346154 
    72.394737         49.538462 
    71.980263         48.730769 
    71.565789         47.923077 
    71.151316         47.115385 
    70.736842         46.307692 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve33 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve33 
    choose top curve    (r_2=1) 
      curve34 
    lines 
    17 14 
    mappingName 
    celltop17 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop17 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 1 17 35 
    lines 
    17 14 16 
    mappingName 
    cell17 
    exit 
  ***---- TFICell 18: Top curve=36, bottom curve=35----- 
  *     .. top    curve:   36   35 
  *     .. bottom curve:   33 
  *   ---top curve 
  spline 
    enter spline points 
    19 
    11.052632         47.384615 
    11.605263         46.397436 
    12.157895         45.410256 
    12.710526         44.423077 
    13.263158         43.435897 
    13.815789         42.448718 
    14.368421         41.461538 
    14.921053         40.474359 
    15.473684         39.487179 
    16.026316         38.500000 
    16.578947         37.512821 
    17.131579         36.525641 
    17.684211         35.538462 
    17.315789         34.641026 
    16.947368         33.743590 
    16.578947         32.846154 
    16.210526         31.948718 
    15.842105         31.051282 
    15.473684         30.153846 
    index 
    lines 
    19 
    shape preserving (toggle) 
    mappingName 
    curve36 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    19 
    0.000000         45.230769 
    0.122807         44.512821 
    0.245614         43.794872 
    0.368421         43.076923 
    0.491228         42.358974 
    0.614035         41.641026 
    0.736842         40.923077 
    0.859649         40.205128 
    0.982456         39.487179 
    1.105263         38.769231 
    1.228070         38.051282 
    1.350877         37.333333 
    1.473684         36.615385 
    1.596491         35.897436 
    1.719298         35.179487 
    1.842105         34.461538 
    1.964912         33.743590 
    2.087719         33.025641 
    2.210526         32.307692 
    index 
    lines 
    19 
    shape preserving (toggle) 
    mappingName 
    curve35 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve35 
    choose top curve    (r_2=1) 
      curve36 
    lines 
    19 13 
    mappingName 
    celltop18 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop18 
    constant depth 
    -9.080000 
    boundary conditions 
    2 1 2 2 18 35 
    lines 
    19 13 16 
    mappingName 
    cell18 
    exit 
  ***---- TFICell 19: Top curve=38, bottom curve=37----- 
  *     .. top    curve:   37   38 
  *     .. bottom curve:   58   59   60   61 
  *   ---top curve 
  spline 
    enter spline points 
    30 
    17.684211         35.538462 
    18.236842         35.538462 
    18.789474         35.538462 
    19.342105         35.538462 
    19.894737         35.538462 
    20.447368         35.538462 
    21.000000         35.538462 
    21.552632         35.538462 
    22.105263         35.538462 
    22.657895         35.538462 
    23.210526         35.538462 
    23.617729         35.878543 
    24.024931         36.218623 
    24.432133         36.558704 
    24.839335         36.898785 
    25.246537         37.238866 
    25.653740         37.578947 
    26.060942         37.919028 
    26.468144         38.259109 
    26.875346         38.599190 
    27.282548         38.939271 
    27.689751         39.279352 
    28.096953         39.619433 
    28.504155         39.959514 
    28.911357         40.299595 
    29.318560         40.639676 
    29.725762         40.979757 
    30.132964         41.319838 
    30.540166         41.659919 
    30.947368         42.000000 
    index 
    lines 
    30 
    shape preserving (toggle) 
    mappingName 
    curve38 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    30 
    15.473684         30.153846 
    16.210526         29.256410 
    16.947368         28.358974 
    17.684211         27.461538 
    18.421053         26.564103 
    19.157895         25.666667 
    19.894737         24.769231 
    21.000000         24.230769 
    22.105263         23.692308 
    23.210526         23.153846 
    24.315789         22.615385 
    25.298246         22.974359 
    26.280702         23.333333 
    27.263158         23.692308 
    28.245614         24.051282 
    29.228070         24.410256 
    30.210526         24.769231 
    31.192982         25.128205 
    32.175439         25.487179 
    33.157895         25.846154 
    33.600000         26.707692 
    34.042105         27.569231 
    34.484211         28.430769 
    34.926316         29.292308 
    35.368421         30.153846 
    35.810526         31.015385 
    36.252632         31.876923 
    36.694737         32.738462 
    37.136842         33.600000 
    37.578947         34.461538 
    index 
    lines 
    30 
    shape preserving (toggle) 
    mappingName 
    curve37 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve37 
    choose top curve    (r_2=1) 
      curve38 
    lines 
    30 10 
    mappingName 
    celltop19 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop19 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 19 35 
    lines 
    30 10 16 
    mappingName 
    cell19 
    exit 
  ***---- TFICell 20: Top curve=40, bottom curve=39----- 
  *     .. top    curve:   56 
  *     .. bottom curve:   62   63 
  *   ---top curve 
  spline 
    enter spline points 
    28 
    37.578947         34.461538 
    38.111111         34.461538 
    38.643275         34.461538 
    39.175439         34.461538 
    39.707602         34.461538 
    40.239766         34.461538 
    40.771930         34.461538 
    41.304094         34.461538 
    41.836257         34.461538 
    42.368421         34.461538 
    42.900585         34.461538 
    43.432749         34.461538 
    43.964912         34.461538 
    44.497076         34.461538 
    45.029240         34.461538 
    45.561404         34.461538 
    46.093567         34.461538 
    46.625731         34.461538 
    47.157895         34.461538 
    47.690058         34.461538 
    48.222222         34.461538 
    48.754386         34.461538 
    49.286550         34.461538 
    49.818713         34.461538 
    50.350877         34.461538 
    50.883041         34.461538 
    51.415205         34.461538 
    51.947368         34.461538 
    index 
    lines 
    28 
    shape preserving (toggle) 
    mappingName 
    curve40 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    28 
    33.157895         25.846154 
    34.162679         25.356643 
    35.167464         24.867133 
    36.172249         24.377622 
    37.177033         23.888112 
    38.181818         23.398601 
    39.186603         22.909091 
    40.191388         22.419580 
    41.196172         21.930070 
    42.200957         21.440559 
    43.205742         20.951049 
    44.210526         20.461538 
    45.108553         20.932692 
    46.006579         21.403846 
    46.904605         21.875000 
    47.802632         22.346154 
    48.700658         22.817308 
    49.598684         23.288462 
    50.496711         23.759615 
    51.394737         24.230769 
    52.292763         24.701923 
    53.190789         25.173077 
    54.088816         25.644231 
    54.986842         26.115385 
    55.884868         26.586538 
    56.782895         27.057692 
    57.680921         27.528846 
    58.578947         28.000000 
    index 
    lines 
    28 
    shape preserving (toggle) 
    mappingName 
    curve39 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve39 
    choose top curve    (r_2=1) 
      curve40 
    lines 
    28 9 
    mappingName 
    celltop20 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop20 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 20 35 
    lines 
    28 9 16 
    mappingName 
    cell20 
    exit 
  ***---- TFICell 21: Top curve=42, bottom curve=41----- 
  *     .. top    curve:   57 
  *     .. bottom curve:   65 
  *   ---top curve 
  spline 
    enter spline points 
    21 
    51.947368         34.461538 
    52.002632         35.215385 
    52.057895         35.969231 
    52.113158         36.723077 
    52.168421         37.476923 
    52.223684         38.230769 
    52.278947         38.984615 
    52.334211         39.738462 
    52.389474         40.492308 
    52.444737         41.246154 
    52.500000         42.000000 
    52.555263         42.753846 
    52.610526         43.507692 
    52.665789         44.261538 
    52.721053         45.015385 
    52.776316         45.769231 
    52.831579         46.523077 
    52.886842         47.276923 
    52.942105         48.030769 
    52.997368         48.784615 
    53.052632         49.538462 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve42 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    21 
    58.578947         28.000000 
    59.186842         28.915385 
    59.794737         29.830769 
    60.402632         30.746154 
    61.010526         31.661538 
    61.618421         32.576923 
    62.226316         33.492308 
    62.834211         34.407692 
    63.442105         35.323077 
    64.050000         36.238462 
    64.657895         37.153846 
    65.265789         38.069231 
    65.873684         38.984615 
    66.481579         39.900000 
    67.089474         40.815385 
    67.697368         41.730769 
    68.305263         42.646154 
    68.913158         43.561538 
    69.521053         44.476923 
    70.128947         45.392308 
    70.736842         46.307692 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve41 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve41 
    choose top curve    (r_2=1) 
      curve42 
    lines 
    21 17 
    mappingName 
    celltop21 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop21 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 21 35 
    lines 
    21 17 16 
    mappingName 
    cell21 
    exit 
  ***---- TFICell 22: Top curve=44, bottom curve=43----- 
  *     .. top    curve:   53 
  *     .. bottom curve:   66   67   68 
  *   ---top curve 
  spline 
    enter spline points 
    29 
    70.736842         46.307692 
    71.210526         46.115385 
    71.684211         45.923077 
    72.157895         45.730769 
    72.631579         45.538462 
    73.105263         45.346154 
    73.578947         45.153846 
    74.052632         44.961538 
    74.526316         44.769231 
    75.000000         44.576923 
    75.473684         44.384615 
    75.947368         44.192308 
    76.421053         44.000000 
    76.894737         43.807692 
    77.368421         43.615385 
    77.842105         43.423077 
    78.315789         43.230769 
    78.789474         43.038462 
    79.263158         42.846154 
    79.736842         42.653846 
    80.210526         42.461538 
    80.684211         42.269231 
    81.157895         42.076923 
    81.631579         41.884615 
    82.105263         41.692308 
    82.578947         41.500000 
    83.052632         41.307692 
    83.526316         41.115385 
    84.000000         40.923077 
    index 
    lines 
    29 
    shape preserving (toggle) 
    mappingName 
    curve44 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    29 
    58.578947         28.000000 
    59.500000         27.282051 
    60.421053         26.564103 
    61.342105         25.846154 
    62.263158         25.128205 
    63.184211         24.410256 
    64.105263         23.692308 
    65.210526         23.756923 
    66.315789         23.821538 
    67.421053         23.886154 
    68.526316         23.950769 
    69.631579         24.015385 
    70.736842         24.080000 
    71.842105         24.144615 
    72.947368         24.209231 
    74.052632         24.273846 
    75.157895         24.338462 
    75.894737         24.912821 
    76.631579         25.487179 
    77.368421         26.061538 
    78.105263         26.635897 
    78.842105         27.210256 
    79.578947         27.784615 
    80.315789         28.358974 
    81.052632         28.933333 
    81.789474         29.507692 
    82.526316         30.082051 
    83.263158         30.656410 
    84.000000         31.230769 
    index 
    lines 
    29 
    shape preserving (toggle) 
    mappingName 
    curve43 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve43 
    choose top curve    (r_2=1) 
      curve44 
    lines 
    29 21 
    mappingName 
    celltop22 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop22 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 1 2 22 35 
    lines 
    29 21 16 
    mappingName 
    cell22 
    exit 
  ***---- TFICell 23: Top curve=46, bottom curve=45----- 
  *     .. top    curve:   70   34 
  *     .. bottom curve:   72   73 
  *   ---top curve 
  spline 
    enter spline points 
    29 
    0.000000         29.076923 
    0.368421         29.615385 
    0.736842         30.153846 
    1.105263         30.692308 
    1.473684         31.230769 
    1.842105         31.769231 
    2.210526         32.307692 
    2.813397         32.209790 
    3.416268         32.111888 
    4.019139         32.013986 
    4.622010         31.916084 
    5.224880         31.818182 
    5.827751         31.720280 
    6.430622         31.622378 
    7.033493         31.524476 
    7.636364         31.426573 
    8.239234         31.328671 
    8.842105         31.230769 
    9.444976         31.132867 
    10.047847         31.034965 
    10.650718         30.937063 
    11.253589         30.839161 
    11.856459         30.741259 
    12.459330         30.643357 
    13.062201         30.545455 
    13.665072         30.447552 
    14.267943         30.349650 
    14.870813         30.251748 
    15.473684         30.153846 
    index 
    lines 
    29 
    shape preserving (toggle) 
    mappingName 
    curve46 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    29 
    0.000000          4.307692 
    0.442105          5.384615 
    0.884211          6.461538 
    1.326316          7.538462 
    1.768421          8.615385 
    2.210526          9.692308 
    2.979405         10.347826 
    3.748284         11.003344 
    4.517162         11.658863 
    5.286041         12.314381 
    6.054920         12.969900 
    6.823799         13.625418 
    7.592677         14.280936 
    8.361556         14.936455 
    9.130435         15.591973 
    9.899314         16.247492 
    10.668192         16.903010 
    11.437071         17.558528 
    12.205950         18.214047 
    12.974828         18.869565 
    13.743707         19.525084 
    14.512586         20.180602 
    15.281465         20.836120 
    16.050343         21.491639 
    16.819222         22.147157 
    17.588101         22.802676 
    18.356979         23.458194 
    19.125858         24.113712 
    19.894737         24.769231 
    index 
    lines 
    29 
    shape preserving (toggle) 
    mappingName 
    curve45 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve45 
    choose top curve    (r_2=1) 
      curve46 
    lines 
    29 24 
    mappingName 
    celltop23 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop23 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 1 23 35 
    lines 
    29 24 16 
    mappingName 
    cell23 
    exit 
  ***---- TFICell 24: Top curve=48, bottom curve=47----- 
  *     .. top    curve:   60   62 
  *     .. bottom curve:   75 
  *   ---top curve 
  spline 
    enter spline points 
    21 
    24.315789         22.615385 
    25.421053         23.019231 
    26.526316         23.423077 
    27.631579         23.826923 
    28.736842         24.230769 
    29.842105         24.634615 
    30.947368         25.038462 
    32.052632         25.442308 
    33.157895         25.846154 
    34.078947         25.397436 
    35.000000         24.948718 
    35.921053         24.500000 
    36.842105         24.051282 
    37.763158         23.602564 
    38.684211         23.153846 
    39.605263         22.705128 
    40.526316         22.256410 
    41.447368         21.807692 
    42.368421         21.358974 
    43.289474         20.910256 
    44.210526         20.461538 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve48 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    21 
    24.315789         15.076923 
    24.923684         14.861538 
    25.531579         14.646154 
    26.139474         14.430769 
    26.747368         14.215385 
    27.355263         14.000000 
    27.963158         13.784615 
    28.571053         13.569231 
    29.178947         13.353846 
    29.786842         13.138462 
    30.394737         12.923077 
    31.002632         12.707692 
    31.610526         12.492308 
    32.218421         12.276923 
    32.826316         12.061538 
    33.434211         11.846154 
    34.042105         11.630769 
    34.650000         11.415385 
    35.257895         11.200000 
    35.865789         10.984615 
    36.473684         10.769231 
    index 
    lines 
    21 
    shape preserving (toggle) 
    mappingName 
    curve47 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve47 
    choose top curve    (r_2=1) 
      curve48 
    lines 
    21 12 
    mappingName 
    celltop24 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop24 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 24 35 
    lines 
    21 12 16 
    mappingName 
    cell24 
    exit 
  ***---- TFICell 25: Top curve=50, bottom curve=49----- 
  *     .. top    curve:   77   76   63 
  *     .. bottom curve:   79   80 
  *   ---top curve 
  spline 
    enter spline points 
    37 
    42.000000          4.307692 
    41.210526          5.230769 
    40.421053          6.153846 
    39.631579          7.076923 
    38.842105          8.000000 
    38.052632          8.923077 
    37.263158          9.846154 
    36.473684         10.769231 
    37.118421         11.576923 
    37.763158         12.384615 
    38.407895         13.192308 
    39.052632         14.000000 
    39.697368         14.807692 
    40.342105         15.615385 
    40.986842         16.423077 
    41.631579         17.230769 
    42.276316         18.038462 
    42.921053         18.846154 
    43.565789         19.653846 
    44.210526         20.461538 
    45.055728         20.904977 
    45.900929         21.348416 
    46.746130         21.791855 
    47.591331         22.235294 
    48.436533         22.678733 
    49.281734         23.122172 
    50.126935         23.565611 
    50.972136         24.009050 
    51.817337         24.452489 
    52.662539         24.895928 
    53.507740         25.339367 
    54.352941         25.782805 
    55.198142         26.226244 
    56.043344         26.669683 
    56.888545         27.113122 
    57.733746         27.556561 
    58.578947         28.000000 
    index 
    lines 
    37 
    shape preserving (toggle) 
    mappingName 
    curve50 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    37 
    59.684211          8.615385 
    59.622807          9.094017 
    59.561404          9.572650 
    59.500000         10.051282 
    59.438596         10.529915 
    59.377193         11.008547 
    59.315789         11.487179 
    59.254386         11.965812 
    59.192982         12.444444 
    59.131579         12.923077 
    59.070175         13.401709 
    59.008772         13.880342 
    58.947368         14.358974 
    58.885965         14.837607 
    58.824561         15.316239 
    58.763158         15.794872 
    58.701754         16.273504 
    58.640351         16.752137 
    58.578947         17.230769 
    58.885965         17.589744 
    59.192982         17.948718 
    59.500000         18.307692 
    59.807018         18.666667 
    60.114035         19.025641 
    60.421053         19.384615 
    60.728070         19.743590 
    61.035088         20.102564 
    61.342105         20.461538 
    61.649123         20.820513 
    61.956140         21.179487 
    62.263158         21.538462 
    62.570175         21.897436 
    62.877193         22.256410 
    63.184211         22.615385 
    63.491228         22.974359 
    63.798246         23.333333 
    64.105263         23.692308 
    index 
    lines 
    37 
    shape preserving (toggle) 
    mappingName 
    curve49 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve49 
    choose top curve    (r_2=1) 
      curve50 
    lines 
    37 18 
    mappingName 
    celltop25 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop25 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 25 35 
    lines 
    37 18 16 
    mappingName 
    cell25 
    exit 
  ***---- TFICell 26: Top curve=52, bottom curve=51----- 
  *     .. top    curve:   79   80 
  *     .. bottom curve:   82   83   84 
  *   ---top curve 
  spline 
    enter spline points 
    18 
    59.684211          8.615385 
    59.546053          9.692308 
    59.407895         10.769231 
    59.269737         11.846154 
    59.131579         12.923077 
    58.993421         14.000000 
    58.855263         15.076923 
    58.717105         16.153846 
    58.578947         17.230769 
    59.192982         17.948718 
    59.807018         18.666667 
    60.421053         19.384615 
    61.035088         20.102564 
    61.649123         20.820513 
    62.263158         21.538462 
    62.877193         22.256410 
    63.491228         22.974359 
    64.105263         23.692308 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve52 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    18 
    76.263158          9.692308 
    77.368421         10.769231 
    78.473684         11.846154 
    79.578947         12.923077 
    79.736842         14.000000 
    79.894737         15.076923 
    80.052632         16.153846 
    80.210526         17.230769 
    80.368421         18.307692 
    80.526316         19.384615 
    80.684211         20.461538 
    79.894737         21.015385 
    79.105263         21.569231 
    78.315789         22.123077 
    77.526316         22.676923 
    76.736842         23.230769 
    75.947368         23.784615 
    75.157895         24.338462 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve51 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve51 
    choose top curve    (r_2=1) 
      curve52 
    lines 
    18 16 
    mappingName 
    celltop26 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop26 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 2 26 35 
    lines 
    18 16 16 
    mappingName 
    cell26 
    exit 
  ***---- TFICell 27: Top curve=54, bottom curve=53----- 
  *     .. top    curve:   85 
  *     .. bottom curve:   84 
  *   ---top curve 
  spline 
    enter spline points 
    11 
    84.000000         31.230769 
    84.000000         30.046154 
    84.000000         28.861538 
    84.000000         27.676923 
    84.000000         26.492308 
    84.000000         25.307692 
    84.000000         24.123077 
    84.000000         22.938462 
    84.000000         21.753846 
    84.000000         20.569231 
    84.000000         19.384615 
    index 
    lines 
    11 
    shape preserving (toggle) 
    mappingName 
    curve54 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    11 
    75.157895         24.338462 
    75.710526         23.950769 
    76.263158         23.563077 
    76.815789         23.175385 
    77.368421         22.787692 
    77.921053         22.400000 
    78.473684         22.012308 
    79.026316         21.624615 
    79.578947         21.236923 
    80.131579         20.849231 
    80.684211         20.461538 
    index 
    lines 
    11 
    shape preserving (toggle) 
    mappingName 
    curve53 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve53 
    choose top curve    (r_2=1) 
      curve54 
    lines 
    11 11 
    mappingName 
    celltop27 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop27 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 1 27 35 
    lines 
    11 11 16 
    mappingName 
    cell27 
    exit 
  ***---- TFICell 28: Top curve=56, bottom curve=55----- 
  *     .. top    curve:   87   72   73 
  *     .. bottom curve:   89   90   74 
  *   ---top curve 
  spline 
    enter spline points 
    33 
    0.000000          0.000000 
    0.000000          1.435897 
    0.000000          2.871795 
    0.000000          4.307692 
    0.442105          5.384615 
    0.884211          6.461538 
    1.326316          7.538462 
    1.768421          8.615385 
    2.210526          9.692308 
    2.947368         10.320513 
    3.684211         10.948718 
    4.421053         11.576923 
    5.157895         12.205128 
    5.894737         12.833333 
    6.631579         13.461538 
    7.368421         14.089744 
    8.105263         14.717949 
    8.842105         15.346154 
    9.578947         15.974359 
    10.315789         16.602564 
    11.052632         17.230769 
    11.789474         17.858974 
    12.526316         18.487179 
    13.263158         19.115385 
    14.000000         19.743590 
    14.736842         20.371795 
    15.473684         21.000000 
    16.210526         21.628205 
    16.947368         22.256410 
    17.684211         22.884615 
    18.421053         23.512821 
    19.157895         24.141026 
    19.894737         24.769231 
    index 
    lines 
    33 
    shape preserving (toggle) 
    mappingName 
    curve56 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    33 
    15.473684          0.000000 
    15.473684          0.861538 
    15.473684          1.723077 
    15.473684          2.584615 
    15.473684          3.446154 
    15.473684          4.307692 
    15.993808          4.941176 
    16.513932          5.574661 
    17.034056          6.208145 
    17.554180          6.841629 
    18.074303          7.475113 
    18.594427          8.108597 
    19.114551          8.742081 
    19.634675          9.375566 
    20.154799         10.009050 
    20.674923         10.642534 
    21.195046         11.276018 
    21.715170         11.909502 
    22.235294         12.542986 
    22.755418         13.176471 
    23.275542         13.809955 
    23.795666         14.443439 
    24.315789         15.076923 
    24.315789         15.830769 
    24.315789         16.584615 
    24.315789         17.338462 
    24.315789         18.092308 
    24.315789         18.846154 
    24.315789         19.600000 
    24.315789         20.353846 
    24.315789         21.107692 
    24.315789         21.861538 
    24.315789         22.615385 
    index 
    lines 
    33 
    shape preserving (toggle) 
    mappingName 
    curve55 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve55 
    choose top curve    (r_2=1) 
      curve56 
    lines 
    33 15 
    mappingName 
    celltop28 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop28 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 1 28 35 
    lines 
    33 15 16 
    mappingName 
    cell28 
    exit 
  ***---- TFICell 29: Top curve=58, bottom curve=57----- 
  *     .. top    curve:   89   90 
  *     .. bottom curve:   92 
  *   ---top curve 
  spline 
    enter spline points 
    18 
    15.473684          0.000000 
    15.473684          1.435897 
    15.473684          2.871795 
    15.473684          4.307692 
    16.105263          5.076923 
    16.736842          5.846154 
    17.368421          6.615385 
    18.000000          7.384615 
    18.631579          8.153846 
    19.263158          8.923077 
    19.894737          9.692308 
    20.526316         10.461538 
    21.157895         11.230769 
    21.789474         12.000000 
    22.421053         12.769231 
    23.052632         13.538462 
    23.684211         14.307692 
    24.315789         15.076923 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve58 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    18 
    26.526316          0.000000 
    27.111455          0.633484 
    27.696594          1.266968 
    28.281734          1.900452 
    28.866873          2.533937 
    29.452012          3.167421 
    30.037152          3.800905 
    30.622291          4.434389 
    31.207430          5.067873 
    31.792570          5.701357 
    32.377709          6.334842 
    32.962848          6.968326 
    33.547988          7.601810 
    34.133127          8.235294 
    34.718266          8.868778 
    35.303406          9.502262 
    35.888545         10.135747 
    36.473684         10.769231 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve57 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve57 
    choose top curve    (r_2=1) 
      curve58 
    lines 
    18 12 
    mappingName 
    celltop29 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop29 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 2 29 35 
    lines 
    18 12 16 
    mappingName 
    cell29 
    exit 
  ***---- TFICell 30: Top curve=60, bottom curve=59----- 
  *     .. top    curve:   93 
  *     .. bottom curve:   77 
  *   ---top curve 
  spline 
    enter spline points 
    15 
    42.000000          0.000000 
    40.894737          0.000000 
    39.789474          0.000000 
    38.684211          0.000000 
    37.578947          0.000000 
    36.473684          0.000000 
    35.368421          0.000000 
    34.263158          0.000000 
    33.157895          0.000000 
    32.052632          0.000000 
    30.947368          0.000000 
    29.842105          0.000000 
    28.736842          0.000000 
    27.631579          0.000000 
    26.526316          0.000000 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve60 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    15 
    42.000000          4.307692 
    41.605263          4.769231 
    41.210526          5.230769 
    40.815789          5.692308 
    40.421053          6.153846 
    40.026316          6.615385 
    39.631579          7.076923 
    39.236842          7.538462 
    38.842105          8.000000 
    38.447368          8.461538 
    38.052632          8.923077 
    37.657895          9.384615 
    37.263158          9.846154 
    36.868421         10.307692 
    36.473684         10.769231 
    index 
    lines 
    15 
    shape preserving (toggle) 
    mappingName 
    curve59 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve59 
    choose top curve    (r_2=1) 
      curve60 
    lines 
    15 14 
    mappingName 
    celltop30 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop30 
    constant depth 
    -9.080000 
    boundary conditions 
    2 2 2 1 30 35 
    lines 
    15 14 16 
    mappingName 
    cell30 
    exit 
  ***---- TFICell 31: Top curve=62, bottom curve=61----- 
  *     .. top    curve:   78 
  *     .. bottom curve:   95 
  *   ---top curve 
  spline 
    enter spline points 
    18 
    42.000000          4.307692 
    43.040248          4.561086 
    44.080495          4.814480 
    45.120743          5.067873 
    46.160991          5.321267 
    47.201238          5.574661 
    48.241486          5.828054 
    49.281734          6.081448 
    50.321981          6.334842 
    51.362229          6.588235 
    52.402477          6.841629 
    53.442724          7.095023 
    54.482972          7.348416 
    55.523220          7.601810 
    56.563467          7.855204 
    57.603715          8.108597 
    58.643963          8.361991 
    59.684211          8.615385 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve62 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    18 
    42.000000          0.000000 
    42.715170          0.000000 
    43.430341          0.000000 
    44.145511          0.000000 
    44.860681          0.000000 
    45.575851          0.000000 
    46.291022          0.000000 
    47.006192          0.000000 
    47.721362          0.000000 
    48.436533          0.000000 
    49.151703          0.000000 
    49.866873          0.000000 
    50.582043          0.000000 
    51.297214          0.000000 
    52.012384          0.000000 
    52.727554          0.000000 
    53.442724          0.000000 
    54.157895          0.000000 
    index 
    lines 
    18 
    shape preserving (toggle) 
    mappingName 
    curve61 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve61 
    choose top curve    (r_2=1) 
      curve62 
    lines 
    18 10 
    mappingName 
    celltop31 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop31 
    constant depth 
    -9.080000 
    boundary conditions 
    2 1 2 2 31 35 
    lines 
    18 10 16 
    mappingName 
    cell31 
    exit 
  ***---- TFICell 32: Top curve=64, bottom curve=63----- 
  *     .. top    curve:   96 
  *     .. bottom curve:   98   99  100 
  *   ---top curve 
  spline 
    enter spline points 
    12 
    54.157895          0.000000 
    54.660287          0.783217 
    55.162679          1.566434 
    55.665072          2.349650 
    56.167464          3.132867 
    56.669856          3.916084 
    57.172249          4.699301 
    57.674641          5.482517 
    58.177033          6.265734 
    58.679426          7.048951 
    59.181818          7.832168 
    59.684211          8.615385 
    index 
    lines 
    12 
    shape preserving (toggle) 
    mappingName 
    curve64 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    12 
    77.368421          0.000000 
    78.105263          1.076923 
    78.842105          2.153846 
    79.578947          3.230769 
    80.131579          4.846154 
    80.684211          6.461538 
    79.947368          7.000000 
    79.210526          7.538462 
    78.473684          8.076923 
    77.736842          8.615385 
    77.000000          9.153846 
    76.263158          9.692308 
    index 
    lines 
    12 
    shape preserving (toggle) 
    mappingName 
    curve63 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve63 
    choose top curve    (r_2=1) 
      curve64 
    lines 
    12 23 
    mappingName 
    celltop32 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop32 
    constant depth 
    -9.080000 
    boundary conditions 
    1 2 2 2 32 35 
    lines 
    12 23 16 
    mappingName 
    cell32 
    exit 
  ***---- TFICell 33: Top curve=66, bottom curve=65----- 
  *     .. top    curve:  100   82   83 
  *     .. bottom curve:  104 
  *   ---top curve 
  spline 
    enter spline points 
    17 
    80.684211          6.461538 
    79.578947          7.269231 
    78.473684          8.076923 
    77.368421          8.884615 
    76.263158          9.692308 
    77.092105         10.500000 
    77.921053         11.307692 
    78.750000         12.115385 
    79.578947         12.923077 
    79.717105         13.865385 
    79.855263         14.807692 
    79.993421         15.750000 
    80.131579         16.692308 
    80.269737         17.634615 
    80.407895         18.576923 
    80.546053         19.519231 
    80.684211         20.461538 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve66 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    17 
    84.000000          7.538462 
    84.000000          8.278846 
    84.000000          9.019231 
    84.000000          9.759615 
    84.000000         10.500000 
    84.000000         11.240385 
    84.000000         11.980769 
    84.000000         12.721154 
    84.000000         13.461538 
    84.000000         14.201923 
    84.000000         14.942308 
    84.000000         15.682692 
    84.000000         16.423077 
    84.000000         17.163462 
    84.000000         17.903846 
    84.000000         18.644231 
    84.000000         19.384615 
    index 
    lines 
    17 
    shape preserving (toggle) 
    mappingName 
    curve65 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve65 
    choose top curve    (r_2=1) 
      curve66 
    lines 
    17 6 
    mappingName 
    celltop33 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop33 
    constant depth 
    -9.080000 
    boundary conditions 
    2 1 2 2 33 35 
    lines 
    17 6 16 
    mappingName 
    cell33 
    exit 
  ***---- TFICell 34: Top curve=68, bottom curve=67----- 
  *     .. top    curve:   98   99 
  *     .. bottom curve:  102 
  *   ---top curve 
  spline 
    enter spline points 
    7 
    77.368421          0.000000 
    78.105263          1.076923 
    78.842105          2.153846 
    79.578947          3.230769 
    79.947368          4.307692 
    80.315789          5.384615 
    80.684211          6.461538 
    index 
    lines 
    7 
    shape preserving (toggle) 
    mappingName 
    curve68 
    exit 
  * 
  *   ---bottom curve 
  spline 
    enter spline points 
    7 
    84.000000          0.000000 
    84.000000          1.256410 
    84.000000          2.512821 
    84.000000          3.769231 
    84.000000          5.025641 
    84.000000          6.282051 
    84.000000          7.538462 
    index 
    lines 
    7 
    shape preserving (toggle) 
    mappingName 
    curve67 
    exit 
  * 
  tfi 
    choose bottom curve (r_2=0) 
      curve67 
    choose top curve    (r_2=1) 
      curve68 
    lines 
    7 6 
    mappingName 
    celltop34 
    exit 
  * 
  depth mapping 
    extend depth from which mapping? 
    celltop34 
    constant depth 
    -9.080000 
    boundary conditions 
    1 1 2 2 34 35 
    lines 
    7 6 16 
    mappingName 
    cell34 
    exit 
  exit this menu 
* 
generate an overlapping grid 
  cell1 
  cell2 
  cell3 
  cell4 
  cell5 
  cell6 
  cell7 
  cell8 
  cell9 
  cell10 
  cell11 
  cell12 
  cell13 
  cell14 
  cell15 
  cell16 
  cell17 
  cell18 
  cell19 
  cell20 
  cell21 
  cell22 
  cell23 
  cell24 
  cell25 
  cell26 
  cell27 
  cell28 
  cell29 
  cell30 
  cell31 
  cell32 
  cell33 
  cell34 
  done choosing mappings 
  change parameters 
    exit 
  set view:0 -0.415179 -0.21875 0 7.22581 1 0 0 0 1 0 0 0 1
  reset:0
  set view:0 -0.28125 -0.65625 0 8 1 0 0 0 1 0 0 0 1
  reset:0
  set view:0 0.457589 0.232143 0 14 1 0 0 0 1 0 0 0 1
  reset:0
  change parameters
    prevent hole cutting
      all
      all
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
hirose-tfi-full-150k.hdf
hirose-tfi
exit
