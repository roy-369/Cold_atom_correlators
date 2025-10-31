#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:11:44 2022

@author: roy.369
"""


import numpy as np
import pickle

import os
import sys
from scipy.interpolate import griddata
from scipy import integrate



def arr_loc(val,X):
   ind = 0
   for i in range(len(X)):
      if(round(float(X[i]),2) == round(val,2)):
         ind = i
         break
   return ind


def normal_average(data_set,data_set_std):

   data_avg = 0
   data_avg_var = 0
   
   a = 0
   
   for i in range(len(data_set)):
       data_avg = data_avg+data_set[i]
       data_avg_var = data_avg_var+data_set_std[i]*data_set_std[i]
       a= a+1
   

   return data_avg/a, np.sqrt(data_avg_var/(a*a))
   
   
     


def thermodynamic_measurements_average(Text_dir_eqm,N,u,mu,dtau,L,run_no):

   
   #Average_up_sign = []
   #Average_dn_sign = []
   Average_total_sign = []
   Average_density = []
   Average_up_occupancy = []
   Average_dn_occupancy = []
   Average_energy = []
   Average_kinetic_energy = []
   Average_Nup_Ndn = []
   AF_corr_func_xx = []
   AF_corr_func_zz = []
   Ferro_corr_func_xx = []
   Ferro_corr_func_zz = []


   #Average_up_sign_var = []
   #Average_dn_sign_var = []
   Average_total_sign_var = []
   Average_density_var = []
   Average_up_occupancy_var = []
   Average_dn_occupancy_var = []
   Average_energy_var = []
   Average_kinetic_energy_var = []
   Average_Nup_Ndn_var = []
   AF_corr_func_xx_var = []
   AF_corr_func_zz_var = []
   Ferro_corr_func_xx_var = []
   Ferro_corr_func_zz_var = []

   run_counter = 0
   r_counter = 0
   while(run_counter<run_no):
      realization = str(run_counter)

      run_counter=run_counter+1
      thermodynamic_data_file = '%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eqm,N,u,mu,dtau,L,realization)
      if os.path.exists(thermodynamic_data_file):
         with open(thermodynamic_data_file, 'rb') as infile:
                 sys_measure = pickle.load(infile)
         if bool(sys_measure):        
            r_counter=r_counter+1
    

             
            #print(realization,mu,L,u,"realization,mu,L,u")
            #up_sign = sys_measure['Average up sign']
            #down_sign = sys_measure['Average dn sign']
            
            total_sign = sys_measure['Average total sign']
            total_sign_s = total_sign.split(' ')
            Average_total_sign.append(float(total_sign_s[0].strip(' ')))
            Average_total_sign_var.append((float(total_sign_s[-1].strip(' ')))**2)
            
            density = sys_measure['Average density']
            density_s = density.split(' ')  
            Average_density.append(float(density_s[0].strip(' '))) 
            Average_density_var.append((float(density_s[-1].strip(' ')))**2)    


            up_occupancy = sys_measure['Average up occupancy']
            up_occupancy_s = up_occupancy.split(' ')
            Average_up_occupancy.append(float(up_occupancy_s[0].strip(' ')))
            Average_up_occupancy_var.append((float(up_occupancy_s[-1].strip(' ')))**2)


            dn_occupancy = sys_measure['Average dn occupancy']
            dn_occupancy_s = dn_occupancy.split(' ')
            Average_dn_occupancy.append(float(dn_occupancy_s[0].strip(' ')))
            Average_dn_occupancy_var.append((float(dn_occupancy_s[-1].strip(' ')))**2)


            Energy = sys_measure['Average Energy']
            Energy_s = Energy.split(' ')
            Average_energy.append(float(Energy_s[0].strip(' ')))
            Average_energy_var.append((float(Energy_s[-1].strip(' ')))**2)


            Kinetic_energy = sys_measure['Average Kinetic Energy']
            Kinetic_energy_s = Kinetic_energy.split(' ')
            Average_kinetic_energy.append(float(Kinetic_energy_s[0].strip(' ')))
            Average_kinetic_energy_var.append((float(Kinetic_energy_s[-1].strip(' ')))**2)


            Nup_Ndn = sys_measure['Average Nup*Ndn']
            Nup_Ndn_s = Nup_Ndn.split(' ')
            Average_Nup_Ndn.append(float(Nup_Ndn_s[0].strip(' ')))
            Average_Nup_Ndn_var.append((float(Nup_Ndn_s[-1].strip(' ')))**2)


            AF_correlation_function_xx = sys_measure['AF correlation function (xx)']
            AF_correlation_function_xx_s = AF_correlation_function_xx.split(' ')
            #print(AF_correlation_function_xx_s,"AF_xx")
            AF_corr_func_xx.append(float(AF_correlation_function_xx_s[0].strip(' ')))
            AF_corr_func_xx_var.append((float(AF_correlation_function_xx_s[-1].strip(' ')))**2)


            AF_correlation_function_zz = sys_measure['AF correlation function (zz)']
            AF_correlation_function_zz_s = AF_correlation_function_zz.split(' ')
            #print(AF_correlation_function_zz_s,"AF_zz")
            AF_corr_func_zz.append(float(AF_correlation_function_zz_s[0].strip(' '))) 
            AF_corr_func_zz_var.append((float(AF_correlation_function_zz_s[-1].strip(' ')))**2) 


            Ferro_correlation_function_xx = sys_measure['Ferro corr. func. (xx)']
            Ferro_correlation_function_xx_s = Ferro_correlation_function_xx.split(' ')
            Ferro_corr_func_xx.append(float(Ferro_correlation_function_xx_s[0].strip(' ')))
            Ferro_corr_func_xx_var.append((float(Ferro_correlation_function_xx_s[-1].strip(' ')))**2)


            Ferro_correlation_function_zz = sys_measure['Ferro corr. func. (zz)']
            Ferro_correlation_function_zz_s = Ferro_correlation_function_zz.split(' ')
            Ferro_corr_func_zz.append(float(Ferro_correlation_function_zz_s[0].strip(' ')))
            Ferro_corr_func_zz_var.append((float(Ferro_correlation_function_zz_s[-1].strip(' ')))**2)


            #Average_up_sign.append(float(up_sign_s[0].strip(' ')))
            #Average_dn_sign.append(float(down_sign_s[0].strip(' ')))

            
            #Average_up_sign_var.append((float(up_sign_s[8].strip(' ')))**2)
            #Average_dn_sign_var.append((float(down_sign_s[8].strip(' ')))**2)

   if(r_counter>0):
      #Up_sign_avg, Up_sign_std  = weighted_average(Average_up_sign, Average_up_sign_var)
      #Dn_sign_avg, Dn_sign_std  = weighted_average(Average_dn_sign, Average_dn_sign_var)
      Total_sign_avg, Total_sign_std  = normal_average(Average_total_sign, Average_total_sign_var)
      Density_avg, Density_std = normal_average(Average_density, Average_density_var)
      Up_occupancy_avg, Up_occupancy_std = normal_average(Average_up_occupancy, Average_up_occupancy_var)
      Dn_occupancy_avg, Dn_occupancy_std = normal_average(Average_dn_occupancy, Average_dn_occupancy_var)
      Energy_avg, Energy_std = normal_average(Average_energy, Average_energy_var)
      Kinetic_Energy_avg, Kinetic_Energy_std = normal_average(Average_kinetic_energy, Average_kinetic_energy_var)
      Nup_Ndn_avg, Nup_Ndn_std = normal_average(Average_Nup_Ndn, Average_Nup_Ndn_var)
      AF_corr_func_xx_avg, AF_corr_func_xx_std = normal_average(AF_corr_func_xx, AF_corr_func_xx_var)
      AF_corr_func_zz_avg, AF_corr_func_zz_std = normal_average(AF_corr_func_zz, AF_corr_func_zz_var)
      Ferro_corr_func_xx_avg, Ferro_corr_func_xx_std = normal_average(Ferro_corr_func_xx, Ferro_corr_func_xx_var)
      Ferro_corr_func_zz_avg, Ferro_corr_func_zz_std = normal_average(Ferro_corr_func_zz, Ferro_corr_func_zz_var)


      Sys_measure_avg = {}
   
      #Sys_measure_avg['Up sign averaged'] = Up_sign_avg
      #Sys_measure_avg['Down sign averaged'] = Dn_sign_avg
      Sys_measure_avg['Total sign averaged'] = Total_sign_avg
      Sys_measure_avg['Density averaged'] = Density_avg
      Sys_measure_avg['Up spin occupancy averaged'] = Up_occupancy_avg
      Sys_measure_avg['Down spin occupancy averaged'] = Dn_occupancy_avg
      Sys_measure_avg['Total energy averaged'] = Energy_avg   
      Sys_measure_avg['Kinetic energy averaged'] = Kinetic_Energy_avg
      Sys_measure_avg['Doublon number averaged'] = Nup_Ndn_avg
      Sys_measure_avg['XX AF structure factor averaged'] = AF_corr_func_xx_avg
      Sys_measure_avg['ZZ AF structure factor averaged'] = AF_corr_func_zz_avg
      Sys_measure_avg['XX Ferro structure factor averaged'] = Ferro_corr_func_xx_avg
      Sys_measure_avg['ZZ Ferro structure factor averaged'] = Ferro_corr_func_zz_avg


      filename_equal_time_measurements_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,u,mu,dtau,L)
      data_equal_time_measurements_avg = Sys_measure_avg
      with open(filename_equal_time_measurements_avg, 'wb') as outfile:
           pickle.dump(data_equal_time_measurements_avg, outfile, pickle.HIGHEST_PROTOCOL)


      Sys_measure_std = {}

      #Sys_measure_std['Up sign standard deviation'] = Up_sign_std
      #Sys_measure_std['Down sign standard deviation'] = Dn_sign_std
      Sys_measure_std['Total sign standard deviation'] = Total_sign_std
      Sys_measure_std['Density standard deviation'] = Density_std
      Sys_measure_std['Up spin occupancy standard deviation'] = Up_occupancy_std
      Sys_measure_std['Down spin occupancy standard deviation'] = Dn_occupancy_std
      Sys_measure_std['Total energy standard deviation'] = Energy_std   
      Sys_measure_std['Kinetic energy standard deviation'] = Kinetic_Energy_std
      Sys_measure_std['Doublon number standard deviation'] = Nup_Ndn_std
      Sys_measure_std['XX AF structure factor standard deviation'] = AF_corr_func_xx_std
      Sys_measure_std['ZZ AF structure factor standard deviation'] = AF_corr_func_zz_std
      Sys_measure_std['XX Ferro structure factor standard deviation'] = Ferro_corr_func_xx_std
      Sys_measure_std['ZZ Ferro structure factor standard deviation'] = Ferro_corr_func_zz_std
 


      filename_equal_time_measurements_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,u,mu,dtau,L)
      data_equal_time_measurements_std = Sys_measure_std
      with open(filename_equal_time_measurements_std, 'wb') as outfile:
          pickle.dump(data_equal_time_measurements_std, outfile, pickle.HIGHEST_PROTOCOL)
                  
      print(r_counter,mu,L,u,"no of realization,mu,L,u")
   else:
       print("no realization found at", mu)

       
def main(total,cmdargs):
    if(total!=4):
        raise ValueError('missing args')

    N = cmdargs[1]
    Dtau = cmdargs[2]
    Run_no = int(cmdargs[3])

    U = ["1.00","1.50","2.00","2.50","3.00","3.50","4.00","4.50","5.00","5.50","6.00","6.50","7.00","7.50","8.00","8.50","9.00","9.50","10.00"]
    Mu = ["0.00"]
    Trot = ["18","20","22","24","26","28","30","32","34","36","38","40","42","44","46","48","50","52","54","56","58","60"]
    
    
    
    Text_dir_main = "../../Text_files/Text_files_N_%s"%N
       
#==============================================Averaging over multiple realizations ================================================================

    for u_points in range(len(U)):
        for i in range(len(Mu)):
            for j in range(len(Trot)):
                Text_dir_eqm = '%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,U[u_points],Dtau,Mu[i],Dtau,Trot[j])
                thermodynamic_measurements_average(Text_dir_eqm,N,U[u_points],Mu[i],Dtau,Trot[j],Run_no)


    
      
if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)












