#!/usr/bin/gnuplot

set key top
set xlabel "field, T"
set ylabel "temperature, K"

fname = "demag1.dat"

plot\
  fname u 2:3 w l title "N",\
  fname u 2:4 w l title "E",\

pause -1