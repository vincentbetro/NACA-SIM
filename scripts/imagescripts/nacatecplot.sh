#!/usr/bin/sh
for file in TecplotEuler_*.dat
do
tecplot -mesa naca0012.lay $file -b -p export.mcr -y `basename $file .dat`.png
done
