#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=10
dtau=0.05
run_no=25

for U in 3.00
do
      for Mu in 0.00 
      do
             for L in 18 20 22 24 26 28 30 32 34 36 38 40 42 44 46 48 50 52 54 56 58 60 
             do
           
                   python3 ./Data_read_legacy_code_equal_time_half_filling.py $N $U $Mu $L $dtau $run_no &> eq_time_data_read_half_filling_output.out
                   python3 ./Data_read_legacy_code_unequal_time_half_filling.py $N $U $Mu $L $dtau $run_no &> un_eq_time_data_read_half_filling_output.out                    
             done
      echo "U_$U, Mu_$Mu, L_$L"
      done
      #echo "U_$U, Mu_$Mu, L_$L"
done



