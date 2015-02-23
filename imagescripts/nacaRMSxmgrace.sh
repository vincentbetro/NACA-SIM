#!/usr/bin/sh
for file in RMSEulerIter_*.dat
do
  gracebat $file -p /ibrix-scr/vbetro/RMS.params
done

