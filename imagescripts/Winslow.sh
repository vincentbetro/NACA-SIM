#!/bin/sh
gnuplot << EOF
set terminal png
set output "$1.png"
set xrange[-0.5:1.5]
set yrange[-1.0:1.0]
plot "$1.dat" w l notitle
EOF
