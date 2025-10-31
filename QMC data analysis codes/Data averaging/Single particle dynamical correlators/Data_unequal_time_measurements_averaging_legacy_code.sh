#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=16
Dtau=0.05
Run_no=50
Nmax=256


for U in 3.00 3.50 4.00 4.50 5.00 5.50 6.00 #6.50 7.00 7.50 8.00 8.50 9.00 9.50 10.00 #3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50 7.00
do
      for Mu in 0.00
      do
             for L in 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 
             do      

                 python3 ./Data_unequal_time_green_function_average_legacy_code_imaginary_time.py $N $U $Mu $L $Dtau $Run_no &> un_eq_average_output_imaginary_time_v2.out
                 #python3 ./Data_unequal_time_green_function_average_legacy_code_spline_matsubara_frequency.py $N $U $Mu $L $Dtau $Run_no $Nmax &> un_eq_average_output_spline_interpolation_matsubara_frequency.out

		# python3 ./Data_unequal_time_measurements_average_legacy_code_spline_matsubara_frequency_longer_sweep.py $N $U $Mu $L $Dtau $Run_no $Nmax &> un_eq_average_output_spline_interpolation_matsubara_frequency_longer_sweep.out


      echo "U_$U, Mu_$Mu, L_$L"
      done
      done

     # echo "U_$U, Mu_$Mu, L_$L"
done



