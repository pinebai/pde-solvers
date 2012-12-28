* Epithelial cell grid
create mappings
***---- TFICell 1: Top curve=2, bottom curve=1-----
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
   curve2
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
   curve1
   exit
*
 tfi
   choose bottom curve (r_2=0)
     curve1
   choose top curve    (r_2=1)
     curve2
   lines
   28 9
   mappingName
   cell1
     boundary conditions
     1 1 1 1 0 0
   exit
   *
exit this menu
*
generate an overlapping grid
cell1
done choosing mappings
change parameters
*** non cutting/non interpolating boundaries
 prevent hole cutting
    all
    all
    done
 prevent interpolation
    all
    all
    done
exit
*
compute overlap
  exit
*
  save a grid
    Grids/debug-simple-hirose-just-cell20.hdf
    epithelial cell with non interpolated boundaries
  exit
