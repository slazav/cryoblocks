#!/usr/bin/gnuplot

set key top
set xlabel "time, h"
set ylabel "temperature, mK"

fname = "doc_example.dat"


plot []\
  fname u ($1/3600):($3*1000) w l title "T(MC)",\
  fname u ($1/3600):($4*1000) w l title "T(Cu)",\
  fname u ($1/3600):($5*1000) w l title "T(He)",\

pause -1

set terminal png
set output "doc_example.png"
replot