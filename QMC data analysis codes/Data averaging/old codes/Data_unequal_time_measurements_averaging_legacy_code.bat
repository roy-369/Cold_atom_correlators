#!/bin/sh
#SBATCH --job-name=l18-Mu12.00
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 8:00:00
#SBATCH --mem=1G
#SBATCH --mail-type=FAIL                        
#SBATCH --mail-user=roy.369@osu.edu                             
hostname
#SBATCH --no-requeue
module load gnu/9.1.0
cd $SLURM_SUBMIT_DIR
time

N=10
Dtau=0.05
Run_no=20

for U in 6.00 8.00
do
      for Mu in $(seq -8.00 0.00 8.00) 
      do
             for L in $(seq 8 2 42)
             do       
                 python3 ./Data_unequal_time_measurements_average_legacy_code.py $N $U $Mu $L $Dtau $Run_no &> un_eq_average_output.out                  
             done
      echo "U_$U, Mu_$Mu, L_$L"
      done

     # echo "U_$U, Mu_$Mu, L_$L"
done

time

