* Epithelial cell grid 
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
    cell1 
    boundary conditions 
    1 2 2 1 0 0 
    exit 
  * 
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
    cell2 
    boundary conditions 
    1 2 2 2 0 0 
    exit 
  * 
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
    cell3 
    boundary conditions 
    2 2 2 1 0 0 
    exit 
  * 
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
    cell4 
    boundary conditions 
    2 2 1 2 0 0 
    exit 
  * 
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
    cell5 
    boundary conditions 
    1 1 1 1 0 0 
    exit 
  * 
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
    cell6 
    boundary conditions 
    1 2 1 1 0 0 
    exit 
  * 
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
    2 1 1 2 0 0 
    exit 
  * 
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
    cell8 
    boundary conditions 
    2 1 1 2 0 0 
    exit 
  * 
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
    cell9 
    boundary conditions 
    2 1 1 2 0 0 
    exit 
  * 
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
    cell10 
    boundary conditions 
    2 1 2 1 0 0 
    exit 
  * 
  exit this menu 
* 
generate an overlapping grid 
  cell4 
  done choosing mappings 
  compute overlap
  exit
save a grid
test2d-cell4.hdf
hirose 2d
exit
