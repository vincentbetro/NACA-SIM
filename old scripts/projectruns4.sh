#!/bin/bash

alpha=(-10 -9 -8 -7 -6)

#this allows us to move the negative around
alpha2=(10 9 8 7 6 5 4 3 2 1)

mach=(0.55 0.65 0.75 0.85 0.95 1.05 1.15 1.25)

order=(1 2)

ptiter=(20 40 60 80 100 120 140 160 180 200)

cfl=(50 100 150 200 250 300 350 400 450 500)

type=(-1 0 1)

#will all be generated in /tmp/vbetro and moved accordingly
#do three refinements
for (( f = 1; f < 3 ; f++ ))
  do

  for (( k = 0; k < ${#alpha[@]} ; k++ ))
    do

     for (( i = 0; i < ${#order[@]} ; i++ ))
       do

       for (( j = 0; j < ${#mach[@]} ; j++ ))
         do

         for (( m = 0; m < ${#cfl[@]} ; m++ ))
           do

           for (( n = 0; n < ${#ptiter[@]} ; n++ ))
             do
    
           # run EULER
           if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
             ./EULER << MYSOLVER
             ${alpha[$k]}
             ${mach[$j]}
             1.4
             naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0${type[$f]}.mesh
             2 3
             1
             1.0
             ${cfl[$m]}
             100
             ${order[$i]}
             ${ptiter[$n]}
             10000
             1.0
             10.e-15
             MYSOLVER
           elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then  
             ./EULER << MYSOLVER
             ${alpha[$k]}
             ${mach[$j]}
             1.4
             naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0${type[$f]}.mesh
             2 3
             1
             1.0
             ${cfl[$m]}
             100
             ${order[$i]}
             ${ptiter[$n]}
             10000
             1.0
             10.e-15
             MYSOLVER
           elif ((${alpha[$k]} == -10)); then
             ./EULER << MYSOLVER
             ${alpha[$k]}
             ${mach[$j]}
             1.4
             naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0${type[$f]}.mesh
             2 3
             1
             1.0
             ${cfl[$m]}
             100
             ${order[$i]}
             ${ptiter[$n]}
             10000
             1.0
             10.e-15
             MYSOLVER
           elif ((${alpha[$k]} == 10)); then
             ./EULER << MYSOLVER
             ${alpha[$k]}
             ${mach[$j]}
             1.4
             naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0${type[$f]}.mesh
             2 3
             1
             1.0
             ${cfl[$m]}
             100
             ${order[$i]}
             ${ptiter[$n]}
             10000
             1.0
             10.e-15
             MYSOLVER
           fi

             if (( $? == 0 && ${order[$i]} == 2 )); then

               if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
               ./EULER << MYSOLVER1
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0${type[$f]}.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               1
               ${ptiter[$n]}
               100
               1.0
               1.0e-15
               MYSOLVER1

               ./EULER << MYSOLVER2
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0${type[$f]}_out.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               ${order[$i]}
               ${ptiter[$n]}
               10000
               1.0
               1.0e-15
               MYSOLVER2

               elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then

               ./EULER << MYSOLVER1
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_+0${alpha2[$k]}_0${type[$f]}.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               1
               ${ptiter[$n]}
               100
               1.0
               1.0e-15
               MYSOLVER1

               ./EULER << MYSOLVER2
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_+0${alpha2[$k]}_0${type[$f]}_out.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               ${order[$i]}
               ${ptiter[$n]}
               10000
               1.0
               1.0e-15
               MYSOLVER2

               elif ((${alpha[$k]} == -10)); then

               ./EULER << MYSOLVER1
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_${alpha2[$k]}_0${type[$f]}.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               1
               ${ptiter[$n]}
               100
               1.0
               1.0e-15
               MYSOLVER1

               ./EULER << MYSOLVER2
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_${alpha2[$k]}_0${type[$f]}_out.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               ${order[$i]}
               ${ptiter[$n]}
               10000
               1.0
               1.0e-15
               MYSOLVER2

               elif ((${alpha[$k]} == 10)); then

               ./EULER << MYSOLVER1
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_${alpha2[$k]}_0${type[$f]}.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               1
               ${ptiter[$n]}
               100
               1.0
               1.0e-15
               MYSOLVER1

               ./EULER << MYSOLVER2
               ${alpha[$k]}
               ${mach[$j]}
               1.4
               naca0012_${order[$i]}_${mach[$j]}_${alpha2[$k]}_0${type[$f]}_out.mesh
               2 3
               1
               1.0
               ${cfl[$m]}
               100
               ${order[$i]}
               ${ptiter[$n]}
               10000
               1.0
               1.0e-15
               MYSOLVER2

               fi
               #need to annotate as going from 1st to 2nd order
               mv *.dat `basename *.dat .dat`_frstscnd.dat
             fi

           if (( $? == 0 )); then
             mv *.dat `basename *.dat .dat`_nan.dat
           elif (( $? == -1 )); then
             mv *.dat `basename *.dat .dat`_uncvg.dat
           fi

           #now, move .dat files to right place
           if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
             mv *.dat /tmp/vbetro/order${order[$i]}/mach${mach[$j]}/alpha-0${alpha2[$k]}
           elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
             mv *.dat /tmp/vbetro/order${order[$i]}/mach${mach[$j]}/alpha+0${alpha[$k]}
           elif ((${alpha[$k]} == -10)); then
             mv *.dat /tmp/vbetro/order${order[$i]}/mach${mach[$j]}/alpha${alpha[$k]}
           elif ((${alpha[$k]} == 10)); then
             mv *.dat /tmp/vbetro/order${order[$i]}/mach${mach[$j]}/alpha${alpha[$k]}
           fi

           #only do if cvg and on highest ptiter and lowest CFL
           if (( $? == 1 && ${ptiter[$n]} == 200 && ${cfl[$m]} == 50 )); then

           # Run subdivision
           if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0${type[$f]}_out.mesh naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0${type[$f]}_out.mesh << MYSUBDIV
           naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f
           1
           3
           5
           0
           1.0
           1.0
           1.0
           5.0
           5.0
           5.0
           2.0
           2.0
           2.0
           1
           MYSUBDIV
           elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0${type[$f]}_out.mesh naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0${type[$f]}_out.mesh << MYSUBDIV
           naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f
           1
           3
           5
           0
           1.0
           1.0
           1.0
           5.0
           5.0
           5.0
           2.0
           2.0
           2.0
           1
           MYSUBDIV
           elif ((${alpha[$k]} == -10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0${type[$f]}_out.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0${type[$f]}_out.mesh << MYSUBDIV
           naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f
           1
           3
           5
           0
           1.0
           1.0
           1.0
           5.0
           5.0
           5.0
           2.0
           2.0
           2.0
           1
           MYSUBDIV
           elif ((${alpha[$k]} == 10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0${type[$f]}_out.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0${type[$f]}_out.mesh << MYSUBDIV
           naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f
           1
           3
           5
           0
           1.0
           1.0
           1.0
           5.0
           5.0
           5.0
           2.0
           2.0
           2.0
           1
           MYSUBDIV
           fi 

           #rm .dat files
           rm -f *.dat
           #meshfile will be overwritten after smoother by mv command, except _out.mesh that is deleted

           #run optimization-based smoother
           if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f.mesh naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f.mesh << MYSMOOTH1
           naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f
           0
           2
           500
           1
           1
           0
           2 3
           1.0e-6
           MYSMOOTH1
           elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f.mesh naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f.mesh << MYSMOOTH1
           naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f
           0
           2
           500
           1
           1
           0
           2 3
           1.0e-6
           MYSMOOTH1
           elif ((${alpha[$k]} == -10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.mesh << MYSMOOTH1
           naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f
           0
           2
           500
           1
           1
           0
           2 3
           1.0e-6
           MYSMOOTH1
           elif ((${alpha[$k]} == 10)); then
           ./SMOOTH naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.mesh << MYSMOOTH1
           naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f
           0
           2
           500
           1
           1
           0
           2 3
           1.0e-6
           MYSMOOTH1
           fi

           #move mesh.dat to proper folder, fix notation on mesh file for next run
           if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f.dat
           elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f.dat
           elif ((${alpha[$k]} == -10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.dat
           elif ((${alpha[$k]} == 10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f_01.dat /tmp/vbetro/meshes/naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.dat
           fi

           if ((${alpha[$k]} < 0)) && ((${alpha[$k]} > -10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f_01.mesh naca0012_${order[$i]}_${mach[$j]}_-0${alpha2[$k]}_0$f.mesh
           elif ((${alpha[$k]} >= 0)) && ((${alpha[$k]} < 10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f_01.mesh naca0012_${order[$i]}_${mach[$j]}_+0${alpha[$k]}_0$f.mesh
           elif ((${alpha[$k]} == -10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f_01.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.mesh
           elif ((${alpha[$k]} == 10)); then
             mv naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f_01.mesh naca0012_${order[$i]}_${mach[$j]}_${alpha[$k]}_0$f.mesh
           fi

           #remove out files, including orig .dat from subdivision
           rm -f *_out.mesh *.dat

           fi

           done

         done

      done

    done

  done

done
