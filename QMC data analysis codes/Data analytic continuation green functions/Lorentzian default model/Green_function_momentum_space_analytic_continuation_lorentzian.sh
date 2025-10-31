#!/bin/bash

N=10
dt=0.05
Freq=2500
omega_max=25.0
omega_min=-25.0


cd /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation

for U in $(seq 7.00 1.00 10.00)
do
   cd Text_files_N_${N}_U_${U}_dtau_${dt} || exit 

   for Mu in 0.00 #$(seq 0.00 0.20 5.00)
   do

        cd Mu_${Mu} || exit

        for trot_slices in 20 30 40 50 60
        do
               Beta_val=$(echo $dt*$trot_slices | bc)
               Timestep=$(echo $trot_slices+1 | bc)

               cd dtau_${dt}_L_${trot_slices} || exit

               for Gam in 0.5 1.0 2.0
               do

                  mkdir -p Spectral_functions_lorentzian_gamma_${Gam}_nfreq_${Freq}
                  cd Spectral_functions_lorentzian_gamma_${Gam}_nfreq_${Freq} || exit
                 
                     
                  for klab in $(seq 0 1 20)
                  do
                      
                      mkdir -p k_point_$klab
                      cd k_point_$klab

                      cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_${N}_analytic_continuation/Text_files_N_${N}_U_${U}_dtau_${dt}/Mu_${Mu}/dtau_${dt}_L_${trot_slices}/Retarded_green_functions_momentum_space_averaged/Retarded_Green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}.dat ./Retarded_green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat
                      cp /fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Data_processing_codes/Data_analytic_continuation_green_functions/Lorentzian_default_model/time_lorentzian_temp.param ./time_lorentzian.param
 
                      data_file=Retarded_green_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_ac.dat
                      out_file=Spectral_function_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_lorentzian_gamma_${Gam}_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
                      out_file_back_cont=Retarded_green_function_back_continued_momentum_space_avg_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_lorentzian_gamma_${Gam}_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
                      out_file_data=Analytic_continuation_data_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_lorentzian_gamma_${Gam}_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.h5
                      out_file_prob=Posterior_probabilities_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_lorentzian_gamma_${Gam}_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
                      out_file_class_max=Spectral_function_momentum_space_avg_classical_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_lorentzian_gamma_${Gam}_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
                      out_file_hist_max=Spectral_function_momentum_space_avg_historic_mxnt_N_${N}_U_${U}_mu_${Mu}_dtau_${dt}_L_${trot_slices}_k_label_${klab}_lorentzian_gamma_${Gam}_nfreq_${Freq}_omega_max_${omega_max}_omega_min_${omega_min}.dat
 
mydir=$(pwd)
rawjob="$(cat <<EOF
#!/bin/sh
#SBATCH --job-name=l$trot_slices-Mu$Mu
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --partition trivedi
#SBATCH -t 12:00:00                    
#SBATCH --mem=1G
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
                            
sed -i 's/beta/'$Beta_val'/g' ./time_lorentzian.param
sed -i 's/Tstep/'$Timestep'/g' ./time_lorentzian.param
sed -i 's/Omega/'$Freq'/g' ./time_lorentzian.param
sed -i 's/gamma/'$Gam'/g' time_lorentzian.param
sed -i 's/fname/'$data_file'/g' ./time_lorentzian.param
sed -i 's/om_max/'$omega_max'/g' ./time_lorentzian.param
sed -i 's/om_min/'$omega_min'/g' ./time_lorentzian.param

maxent time_lorentzian.param &>k_label_${klab}_output_lorentzian.out
                         
sleep 5
cp time_lorentzian.out.avspec.dat $out_file
rm time_lorentzian.out.avspec.dat
cp time_lorentzian.out.avspec_back.dat $out_file_back_cont
rm time_lorentzian.out.avspec_back.dat
rm time_lorentzian.out.chi2.dat
cp time_lorentzian.out.chispec.dat $out_file_hist_max
rm time_lorentzian.out.chispec.dat
rm time_lorentzian.out.chispec_back.dat
rm time_lorentzian.out.fits.dat
cp time_lorentzian.out.maxspec.dat $out_file_class_max
rm time_lorentzian.out.maxspec.dat
rm time_lorentzian.out.maxspec_back.dat
cp time_lorentzian.out.out.h5 $out_file_data
rm time_lorentzian.out.out.h5
cp time_lorentzian.out.prob.dat $out_file_prob
rm time_lorentzian.out.prob.dat
rm time_lorentzian.out.spex.dat
      
time
EOF
)"

                         echo "$rawjob" &> job.bat
                         sbatch job.bat
                      cd ..
                      done
#loop for k ends here                         
        
                  echo "Gamma_$Gam"
                  cd ..     
                  done
#loop for gamma ends here
                 
               cd ..
               done
#loop for beta ends here
               
             cd ..
             echo "U_$U, Mu_$Mu, L_$trot_slices"
             #sleep 100
             done
#loop for Mu ends here
             
        cd ..               
        echo "U_$U, Mu_$Mu, L_$trot_slices"
        done
#loop for U ends here

                 

