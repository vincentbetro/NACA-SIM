#!/bin/bash

alpha=(-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10)

#this allows us to move the negative around
alpha2=(10 9 8 7 6 5 4 3 2 1)

mach=(0.55 0.65 0.75 0.85 0.95 1.05 1.15 1.25)

order=(1 2)

ptiter=(20 40 60 80 100 120 140 160 180 200)

cfl=(50 100 150 200 250 300 350 400 450 500)

type=(-1 0 1)

#will all be generated in /tmp/vbetro and moved accordingly

#generate initial meshes for all alpha
for (( k = 0; k < ${#alpha[@]} ; k++ ))
    do
    #run L-E smoother to rotate mesh
    if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
     ./SMOOTH naca0012.mesh naca0012.mesh << MYSMOOTH
     naca0012_-0${alpha2[$k]}
     0
     1
     0.25
     45
     1
     1
     ${alpha[$k]}
     2 3
     1.0
     1.0e-14
     MYSMOOTH
    elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
     ./SMOOTH naca0012.mesh naca0012.mesh << MYSMOOTH
     naca0012_+0${alpha[$k]}
     0
     1
     0.25
     45
     1
     1
     ${alpha[$k]}
     2 3
     1.0
     1.0e-14
     MYSMOOTH
    elif ((${alpha[$k]} == -10)); then
     ./SMOOTH naca0012.mesh naca0012.mesh << MYSMOOTH
     naca0012_${alpha[$k]}
     0
     1
     0.25
     45
     1
     1
     ${alpha[$k]}
     2 3
     1.0
     1.0e-14
     MYSMOOTH
    elif ((${alpha[$k]} == 10)); then
     ./SMOOTH naca0012.mesh naca0012.mesh << MYSMOOTH
     naca0012_${alpha[$k]}
     0
     1
     0.25
     45
     1
     1
     ${alpha[$k]}
     2 3
     1.0
     1.0e-14
     MYSMOOTH
    fi
 #now, artifically store for all mach and order
  for (( i = 0; i < ${#order[@]} ; i++ ))
    do
     for (( j = 0; j < ${#mach[@]} ; j++ ))
       do
       #now, store as proper file for future
       if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
          cp naca0012_-0${alpha2[$k]}_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_00.dat
       elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
          cp naca0012_+0${alpha[$k]}_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_00.dat
       elif ((${alpha[$k]} == -10)); then
          cp naca0012_${alpha[$k]}_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_00.dat 
       elif ((${alpha[$k]} == 10)); then
          cp naca0012_${alpha[$k]}_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_00.dat
       fi
       #store here for starting cases
       if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
          cp naca0012_-0${alpha2[$k]}_01.mesh naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_00.mesh
       elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
          cp naca0012_+0${alpha[$k]}_01.mesh naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_00.mesh
       elif ((${alpha[$k]} == -10)); then
          cp naca0012_${alpha[$k]}_01.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_00.mesh 
       elif ((${alpha[$k]} == 10)); then
          cp naca0012_${alpha[$k]}_01.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_00.mesh
       fi
       done
     done
  done
#now, remove all .dat and deprecated mesh files
rm -f *.dat *_01.mesh
