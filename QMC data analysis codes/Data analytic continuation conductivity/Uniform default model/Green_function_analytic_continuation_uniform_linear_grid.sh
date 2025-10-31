#!/bin/bash

N=16
dt=0.05
Freq=3000
omega_max=30
omega_min=-30

cd /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Spectral_function_machine_learning/Text_files/Text_files_N_${N}_analytic_continuation

for U in 8.0
do
   cd Text_files_N_${N}_U_${U}_dtau_${dt} || exit 

   for Mu in 6.50
   do
        cd Mu_${Mu} || exit

        for trot_slices in 40
        do
               Beta_val=$(echo $dt*$trot_slices | bc)
               Timestep=$(echo $trot_slices+1 | bc)

               cd dtau_${dt}_L_${trot_slices} || exit

               mkdir -p Spectral_functions_uniform_nfreq_${Freq}_linear_grid
               cd Spectral_functions_uniform_nfreq_${Freq}_linear_grid || exit

               
               for klab in $(seq 0 1 44)
               do
               
                  mkdir -p k_point_$klab
                  cd k_point_$klab
  
                  cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Spectral_function_machine_learning/Text_files/Text_files_N_${N}_analytic_continuation/Text_files_N_${N}_U_${U}_dtau_${dt}/Mu_${Mu}/dtau_${dt}_L_${trot_slices}/Retarded_green_functions_momentum_space_averaged/Retarded_Green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}.dat ./Retarded_Green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat
                  cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Spectral_function_machine_learning/Data_processing_codes/Data_analytic_continuation/Uniform_default_model/time_uniform_linear_grid_temp.param ./time_uniform_linear_grid.param

data_file=Retarded_Green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat
                         
out_file=Spectral_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}_linear_grid.dat
out_file_back_cont=Retarded_green_function_back_continued_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}_linear_grid.dat
out_file_data=Analytic_continuation_data_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}_linear_grid.h5
out_file_prob=Posterior_probabilities_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}_linear_grid.dat
out_file_class_max=Spectral_function_momentum_space_avg_classical_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}_linear_grid.dat
out_file_hist_max=Spectral_function_momentum_space_avg_historic_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}_linear_grid.dat

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=l$trot_slices-Mu$Mu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH -t 8:00:00                    
#SBATCH --mem=4G
#SBATCH --partition trivedi
#SBATCH --mail-type=FAIL                        
#SBATCH --mail-user=roy.369@osu.edu                             
hostname
#SBATCH --no-requeue
ml load CQMP-maxent
cd \$SLURM_SUBMIT_DIR
time
                     
echo "$(pwd)"
               
echo "$data_file"
                         
echo "$out_file"
                            
sed -i 's/beta/'$Beta_val'/g' ./time_uniform_linear_grid.param
sed -i 's/Tstep/'$Timestep'/g' ./time_uniform_linear_grid.param
sed -i 's/Omega/'$Freq'/g' ./time_uniform_linear_grid.param
sed -i 's/grid_type/'$Grid_type'/g' ./time_unifrom_linear_grid.param 
sed -i 's/fname/'$data_file'/g' ./time_uniform_linear_grid.param
sed -i 's/om_max/'$omega_max'/g' ./time_uniform_linear_grid.param
sed -i 's/om_min/'$omega_min'/g' ./time_uniform_linear_grid.param

maxent time_uniform_linear_grid.param &>k_label_${klab}_output_uniform.out
                         
sleep 5
cp time_uniform_linear_grid.out.avspec.dat $out_file
rm time_uniform_linear_grid.out.avspec.dat
cp time_uniform_linear_grid.out.avspec_back.dat $out_file_back_cont
rm time_uniform_linear_grid.out.avspec_back.dat
rm time_uniform_linear_grid.out.chi2.dat
cp time_uniform_linear_grid.out.chispec.dat $out_file_hist_max
rm time_uniform_linear_grid.out.chispec.dat
rm time_uniform_linear_grid.out.chispec_back.dat
rm time_uniform_linear_grid.out.fits.dat
cp time_uniform_linear_grid.out.maxspec.dat $out_file_class_max
rm time_uniform_linear_grid.out.maxspec.dat
rm time_uniform_linear_grid.out.maxspec_back.dat
cp time_uniform_linear_grid.out.out.h5 $out_file_data
rm time_uniform_linear_grid.out.out.h5
cp time_uniform_linear_grid.out.prob.dat $out_file_prob
rm time_uniform_linear_grid.out.prob.dat
rm time_uniform_linear_grid.out.spex.dat
      
time
EOF
)"

			 echo "$rawjob" &> job.bat
			 sbatch job.bat
                    cd ..
                  done 
#loop for k ends here                      
                  cd ..
                  cd ..
               done
#loop for beta ends now
             cd ..
             done
#loop for Mu ends now
        cd ..               
        echo "U_$U, Mu_$Mu, L_$trot_slices"
        done
#loop for U ends now    
                    
