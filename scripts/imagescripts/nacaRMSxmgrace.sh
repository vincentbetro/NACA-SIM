#!/usr/bin/sh
for file in RMSEulerIter_*.dat
do
  gracebat $file -p RMS.params
done

