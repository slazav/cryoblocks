#!/usr/bin/gnuplot

set key bottom
set xlabel "time, s"
set ylabel "temperature, K"

fname = "zero-c.dat"

plot\
  fname u 1:2 w lp pt 6 title "B1",\
  fname u 1:3 w lp pt 6 title "B2",\
  fname u 1:4 w lp pt 6 title "B3",\
  fname u 1:5 w lp pt 6 title "B4",\
 300

pause -1