set mxtics 5
set mytics 5
set grid xtics ytics mxtics mytics
#set grid
plot 'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:2 title "1 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:5 title "2 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:8 title "3 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:11 title "4 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:14 title "5 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:17 title "6 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:20 title "7 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:23 title "8 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:26 title "9 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:29 title "10 calcium" with linespoints, \
     'OUT/ref_lirinzelwagner_2D_newprobes.sequence' u 1:32 title "11 calcium" with linespoints
