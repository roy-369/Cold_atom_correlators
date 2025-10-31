#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

dtau=0.025
run_no=3

for U in 0.00
do
   for Mu in $(seq 0.00 0.10 8.00) 
   do
       for L in $(seq 10 2 24)
       do           
       
        python3 ./Data_equal_time_measurements_connected_correlators_neighbor_averaged_normal_average_legacy_code.py 10 $U $Mu $L $dtau $run_no &> connected_neighbor_averaged_normal_average_output.out

        echo "U_$U, Mu_$Mu, L_$L"
       
        done
   done
done



