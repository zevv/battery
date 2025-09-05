#!/usr/bin/gnuplot

reset
set grid
set multiplot layout 2, 1
set log x
plot "/tmp/d" u 1:4 w lp, "/tmp/d" u 1:5 w lp axis x1y2

unset log x
plot "/tmp/d" u 2:(-$3) w lp
