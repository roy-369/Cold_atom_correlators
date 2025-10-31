#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=10
Dtau=0.05
Run_no=20

for U in 5.10 5.20 5.30 5.40 #2.00 $(seq 3.20 0.10 5.00) 5.10 5.25 5.50 5.75 6.00 8.00
do
      for Mu in 0.00 
      do
             for L in $(seq 20 2 40)
             do       
                 python3 ./Data_unequal_time_measurements_average_legacy_code.py $N $U $Mu $L $Dtau $Run_no #&> un_eq_average_output.out                  
             done
      echo "U_$U, Mu_$Mu, L_$L"
      done

     # echo "U_$U, Mu_$Mu, L_$L"
done



