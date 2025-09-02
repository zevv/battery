#!/usr/bin/gnuplot -p

reset
set grid
set key off


set multiplot layout 5, 1
set lmargin at screen 0.08
set noxtics

# 0 margin between multiplots
set tmargin 1
set bmargin 1

# title:
set ylabel "I (A)"
plot "/tmp/cell_00.log" u 1 w l, "/tmp/cell_01.log" u 1 w l, "/tmp/cell_02.log" u 1 w l, "/tmp/cell_03.log" u 1 w l

set ylabel "U (V)"
set yrange [2.4:4.3]
plot "/tmp/cell_00.log" u 2 w l, "/tmp/cell_01.log" u 2 w l, "/tmp/cell_02.log" u 2 w l, "/tmp/cell_03.log" u 2 w l

set ylabel "SOC (%)"
set yrange [0.0:1.0]
plot "/tmp/cell_00.log" u 3 w l, "/tmp/cell_01.log" u 3 w l, "/tmp/cell_02.log" u 3 w l, "/tmp/cell_03.log" u 3 w l

set ylabel "T (°C)"
unset yrange
plot "/tmp/cell_00.log" u 4 w l, "/tmp/cell_01.log" u 4 w l, "/tmp/cell_02.log" u 4 w l, "/tmp/cell_03.log" u 4 w l

set ylabel "SOH (%)"
set yrange [0.0:1.0]
plot "/tmp/cell_00.log" u 5 w l, "/tmp/cell_01.log" u 5 w l, "/tmp/cell_02.log" u 5 w l, "/tmp/cell_03.log" u 5 w l

unset multiplot
