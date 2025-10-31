#!/bin/bash

# Nx 2 - 4 ,Nxstep = 1
# Ny 2 - 4 ,Nystep = 1
# t1 0 - 1, t1Step = 0.5
# U 10 - 20, uStep  = 10

N=16
dtau=0.05
run_no=50

for U in 3.50 4.00 4.50 5.00 #$(seq 3.10 0.10 5.00) #3.00 3.50 4.00 4.50 5.00 5.50 6.00 #6.50 7.00 7.50 8.00 8.50 9.00 9.50 10.00 #6.00 #7.00 7.50 8.00 8.50 9.00 9.50 10.00  #5.25 5.50 5.75 6.00 #3.10 3.20 3.30 3.40 3.50 3.60 3.70 3.80 3.90 4.00 4.10 4.20 4.30 4.40 4.50 4.60 4.70 4.80 4.90 5.00
do
      for Mu in -0.20 0.00 #$(seq -0.00 0.20 2.20) # -0.20 0.00 0.20 #$(seq -10.00 0.50 0.00) #$(seq 4.50 0.50 10.00) 
      do
	     for L in 42 #20 22 24 26 28 30 32 34 36 38 40 42 #44 46 48 50 52 54 56 58 60  #42 44 46 48 50 52 54 56 58 60
             do
           
                   python3 ./Data_read_legacy_code_equal_time.py $N $U $Mu $L $dtau $run_no &> eq_time_data_read_output.out
                   python3 ./Data_read_legacy_code_unequal_time.py $N $U $Mu $L $dtau $run_no &> un_eq_time_data_read_output.out                    
             done
      echo "U_$U, Mu_$Mu, L_$L"
      done
      #echo "U_$U, Mu_$Mu, L_$L"
done



