#!/usr/bin/sh
for file in TecplotEuler_*.dat
do
  tecplot -mesa /ibrix-scr/vbetro/naca0012.lay $file -b -p /ibrix-scr/vbetro/export.mcr -y `basename $file .dat`.png
done
