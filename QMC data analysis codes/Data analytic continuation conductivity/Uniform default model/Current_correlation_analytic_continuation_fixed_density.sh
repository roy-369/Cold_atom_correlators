#!/bin/bash

N=10
dt=0.05
Freq=3000
omega_max=30.0
omega_min=-30.0
Ph_case=false

cd /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation_current_correlations

for U in 3.20 3.30 3.40 3.50 3.60 3.70 3.80 3.90 4.00 4.10 4.20 4.30 4.40 4.50 4.60 4.70 4.80 4.90 5.00 5.10 5.30 5.40 5.50 5.60 5.70 6.00
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

                  mkdir -p Current_correlation_functions_analytic_continued_uniform_nfreq_${Freq}
                  cd Current_correlation_functions_analytic_continued_uniform_nfreq_${Freq} || exit
                 
                     
                  for klab in 0 #$(seq 0 1 48)
                  do
                      
                      mkdir -p k_point_$klab
                      cd k_point_$klab

                      cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation_current_correlations/Text_files_N_${N}_U_${U}_dtau_${dt}/Mu_${Mu}/dtau_${dt}_L_${trot_slices}/Current_current_correlation_functions_momentum_space_averaged/Current_current_correlation_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}.dat ./Current_current_correlation_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat

                      cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Data_processing_codes/Data_analytic_continuation_current_correlations/Uniform_default_model/time_uniform_temp.param ./time_uniform.param
 
                      data_file=Current_current_correlation_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat

                      out_file=Conductivity_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
                      
                      out_bose_file=Conductivity_momentum_space_avg_bosonic_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat

                      out_file_back_cont=Current_current_correlation_function_back_continued_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
                      out_file_data=Analytic_continuation_data_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.h5
                      out_file_prob=Posterior_probabilities_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
                      out_file_class_max=Conductivity_momentum_space_avg_classical_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat

                      out_bose_file_class_max=Conductivity_momentum_space_avg_bosonic_classical_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat

                      out_file_hist_max=Conductivity_momentum_space_avg_historic_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_uniform_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
 

mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=l$trot_slices-Mu$Mu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --partition trivedi
#SBATCH -t 24:00:00                    
#SBATCH --mem=4G
#SBATCH --mail-type=FAIL                        
#SBATCH --mail-user=roy.369@osu.edu                             
hostname
#SBATCH --no-requeue
ml load CQMP-maxent/2.0
cd \$SLURM_SUBMIT_DIR
time
                     
echo "$(pwd)"
               
echo "$data_file"
                         
echo "$out_file"
                            
sed -i 's/beta/'$Beta_val'/g' ./time_uniform.param
sed -i 's/Tstep/'$Timestep'/g' ./time_uniform.param
sed -i 's/Omega/'$Freq'/g' ./time_uniform.param
sed -i 's/fname/'$data_file'/g' ./time_uniform.param
sed -i 's/om_max/'$omega_max'/g' ./time_uniform.param
sed -i 's/om_min/'$omega_min'/g' ./time_uniform.param
sed -i 's/ph_case/'$Ph_case'/g' ./time_uniform.param

maxent time_uniform.param &>k_label_${klab}_output_uniform.out
                         
sleep 5
cp time_uniform.out.avspec.dat $out_file
rm time_uniform.out.avspec.dat
cp time_uniform.out.avspec_bose.dat $out_bose_file
rm time_uniform.out.avspec_bose.dat
cp time_uniform.out.avspec_back.dat $out_file_back_cont
rm time_uniform.out.avspec_back.dat
rm time_uniform.out.chi2.dat
cp time_uniform.out.chispec.dat $out_file_hist_max
rm time_uniform.out.chispec.dat
rm time_uniform.out.chispec_back.dat
rm time_uniform.out.fits.dat
cp time_uniform.out.maxspec.dat $out_file_class_max
rm time_uniform.out.maxspec.dat
cp time_uniform.out.maxspec_bose.dat $out_bose_file_class_max
rm time_uniform.out.maxspec_bose.dat
rm time_uniform.out.maxspec_back.dat
cp time_uniform.out.out.h5 $out_file_data
rm time_uniform.out.out.h5
cp time_uniform.out.prob.dat $out_file_prob
rm time_uniform.out.prob.dat
rm time_uniform.out.spex.dat
      
time
EOF
)"

                         echo "$rawjob" &> job.bat
                         sbatch job.bat
                      cd ..
                      done
#loop for k ends here                         
        
                  cd ..     
#loop for sigma ends here
                 
               cd ..
               done
#loop for beta ends here
               
             cd ..
             echo "U_$U, Mu_$Mu, L_$trot_slices"
             #sleep 300
             done
#loop for Mu ends here
             
        cd ..               
        echo "U_$U, Mu_$Mu, L_$trot_slices"
        done
#loop for U ends here

                 

