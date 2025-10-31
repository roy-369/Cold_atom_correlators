#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=16
Dtau=0.05
Run_no=20
Nmax=256


for U in 10.0
do
      for Mu in $(seq 0.00 0.20 5.00)
      do
             for L in 50 60
             do      

                 python3 ./Data_unequal_time_green_function_average_legacy_code_imaginary_time_longer_sweep.py $N $U $Mu $L $Dtau $Run_no &> un_eq_average_output_imaginary_time_longer_sweep.out
                 python3 ./Data_unequal_time_green_function_average_legacy_code_spline_matsubara_frequency_longer_sweep.py $N $U $Mu $L $Dtau $Run_no $Nmax &> un_eq_average_output_spline_interpolation_matsubara_frequency_longer_sweep.out



      echo "U_$U, Mu_$Mu, L_$L"
      done
      done

     # echo "U_$U, Mu_$Mu, L_$L"
done



