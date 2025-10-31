#!/bin/bash

N=10
dt=0.05
Freq=2500
omega_max=25
omega_min=-25


cd /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation

for U in 4.00 4.50 5.00 5.50 6.00 6.50 7.00
do
   cd Text_files_N_${N}_U_${U}_dtau_${dt} || exit 

   for Mu in 0.00
   do
        cd Mu_${Mu} || exit

        for trot_slices in 20 30 40
        do
               Beta_val=$(echo $dt*$trot_slices | bc)
               Timestep=$(echo $trot_slices+1 | bc)

               cd dtau_${dt}_L_${trot_slices} || exit

               
               mkdir -p Spectral_functions_sum_rule_optimized_gaussian_nfreq_${Freq}
               cd Spectral_functions_sum_rule_optimized_gaussian_nfreq_${Freq} || exit

               cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation/Text_files_N_${N}_U_${U}_dtau_${dt}/Mu_${Mu}/dtau_${dt}_L_${trot_slices}/Retarded_green_functions_momentum_space_averaged/Default_model_mean_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}.dat ./Default_model_mean_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}.dat

               cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation/Text_files_N_${N}_U_${U}_dtau_${dt}/Mu_${Mu}/dtau_${dt}_L_${trot_slices}/Retarded_green_functions_momentum_space_averaged/Default_model_variance_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}.dat ./Default_model_variance_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}.dat
                     
               IFS=$'\n' read -d '' -r -a Mean_arr < "./Default_model_mean_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}.dat"
               IFS=$'\n' read -d '' -r -a Variance_arr < "./Default_model_variance_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}.dat"

               for klab in $(seq 0 1 20)
               do
                      
                      mkdir -p k_point_$klab
                      cd k_point_$klab
                      
                      Mean=${Mean_arr[$klab]}
                      Sig=${Variance_arr[$klab]}
                      #echo "Mean_$Mean, Sigma_$Sig"

                      cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation/Text_files_N_${N}_U_${U}_dtau_${dt}/Mu_${Mu}/dtau_${dt}_L_${trot_slices}/Retarded_green_functions_momentum_space_averaged/Retarded_Green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}.dat ./Retarded_green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat
                     
                      cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Data_processing_codes/Data_analytic_continuation_green_functions/Shifted_gaussian_default_model/time_gaussian_temp.param ./time_gaussian.param
 
                      data_file=Retarded_green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat

                      out_file=Spectral_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_gaussian_sum_rule_optimized_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat

                      out_file_back_cont=Retarded_green_function_back_continued_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_gaussian_sum_rule_optimized_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat

                      out_file_data=Analytic_continuation_data_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_gaussian_sum_rule_optimized_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.h5

                      out_file_prob=Posterior_probabilities_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_gaussian_sum_rule_optimized_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat

                      out_file_class_max=Spectral_function_momentum_space_avg_classical_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_gaussian_sum_rule_optimized_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat

                      out_file_hist_max=Spectral_function_momentum_space_avg_historic_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_gaussian_sum_rule_optimized_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
 
mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=l$trot_slices-Mu$Mu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --partition trivedi
#SBATCH -t 16:00:00                    
#SBATCH --mem=2G
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
                            
sed -i 's/beta/'$Beta_val'/g' ./time_gaussian.param
sed -i 's/Tstep/'$Timestep'/g' ./time_gaussian.param
sed -i 's/Omega/'$Freq'/g' ./time_gaussian.param
sed -i 's/sigma/'$Sig'/g' ./time_gaussian.param
sed -i 's/mean/'$Mean'/g' ./time_gaussian.param
sed -i 's/fname/'$data_file'/g' ./time_gaussian.param
sed -i 's/om_max/'$omega_max'/g' ./time_gaussian.param
sed -i 's/om_min/'$omega_min'/g' ./time_gaussian.param

maxent time_gaussian.param &>k_label_${klab}_output_gaussian.out
                         
sleep 5
cp time_gaussian.out.avspec.dat $out_file
rm time_gaussian.out.avspec.dat
cp time_gaussian.out.avspec_back.dat $out_file_back_cont
rm time_gaussian.out.avspec_back.dat
rm time_gaussian.out.chi2.dat
cp time_gaussian.out.chispec.dat $out_file_hist_max
rm time_gaussian.out.chispec.dat
rm time_gaussian.out.chispec_back.dat
rm time_gaussian.out.fits.dat
cp time_gaussian.out.maxspec.dat $out_file_class_max
rm time_gaussian.out.maxspec.dat
rm time_gaussian.out.maxspec_back.dat
cp time_gaussian.out.out.h5 $out_file_data
rm time_gaussian.out.out.h5
cp time_gaussian.out.prob.dat $out_file_prob
rm time_gaussian.out.prob.dat
rm time_gaussian.out.spex.dat
      
time
EOF
)"

                         echo "$rawjob" &> job.bat
                         sbatch job.bat

                      cd ..
                      done
#loop for k ends here                         
        
                  #echo "Sigma_$Sig"
                  cd ..     
                  #done
#Exit out of gaussian spectral function folder
                 
               cd ..
               done
#loop for beta ends here
               
             cd ..
             done
#loop for Mu ends here
             
        cd ..               
        echo "U_$U, Mu_$Mu, L_$trot_slices"
        done
#loop for U ends here

                 

