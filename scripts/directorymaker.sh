#!/bin/bash

alpha=(-10 -9 -8 -7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 5 6 7 8 9 10)

#this allows us to move the negative around
alpha2=(10 9 8 7 6 5 4 3 2 1)

mach=(0.55 0.65 0.75 0.85 0.95 1.05 1.15 1.25)

order=(1 2)

for (( i = 0; i < ${#order[@]} ; i++ ))
do

  for (( j = 0; j < ${#mach[@]} ; j++ ))
  do

    for (( k = 0; k < ${#alpha[@]} ; k++ ))
    do
      
      if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
          mkdir -p order${order[$i]}/mach${mach[$j]}/alpha-0${alpha2[$k]}
      elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
          mkdir -p order${order[$i]}/mach${mach[$j]}/alpha+0${alpha[$k]}
      elif ((${alpha[$k]} == -10)); then
          mkdir -p order${order[$i]}/mach${mach[$j]}/alpha${alpha[$k]}
      elif ((${alpha[$k]} == 10)); then
          mkdir -p order${order[$i]}/mach${mach[$j]}/alpha+${alpha[$k]}
      fi

    done

  done

done
