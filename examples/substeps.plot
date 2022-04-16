#!/usr/bin/gnuplot

set key bottom
set xlabel "time, s"
set ylabel "temperature, K"

fname = "substeps.dat"

plot\
  fname u 1:2 w lp pt 6 title "B1",\
  fname u 1:3 w lp pt 6 title "B2",\

pause -1