#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=16
dtau=0.05
run_no=4

for U in 6.00 #3.00 3.50 4.00 4.50 5.00 5.50 6.00 6.50 7.00 #5.0 7.0 9.0 #1.0 4.0 5.0 6.0 8.0 10.0
do
	for Mu in $(seq -0.20 0.20 2.20) #$(seq -0.20 0.20 2.00) $(seq 4.50 0.50 10.00) #$(seq 2.50 0.50 10.00) #-1.80 -1.60 -1.40 -1.20 -0.80 -0.60 -0.40 -0.20 0.20 0.40 0.60 0.80 1.20 1.40 1.60 1.80
   do
       for L in 10 12 14 16 18 20 22 24 26 28 30 32 34 36 38 40 #42 44 46 48 52 54
       do           
       
        python3 ./Data_equal_time_measurements_connected_correlators_neighbor_averaged_normal_average_legacy_code.py $N $U $Mu $L $dtau $run_no &> connected_neighbor_averaged_normal_average_output.out
        python3 ./Data_equal_time_measurements_neighbor_averaged_normal_average_legacy_code.py $N $U $Mu $L $dtau $run_no &> neighbor_averaged_normal_average_output.out

        echo "U_$U, Mu_$Mu, L_$L"
       
        done
   done
done



