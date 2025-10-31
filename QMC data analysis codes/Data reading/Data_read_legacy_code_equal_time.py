#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 26 18:07:21 2022

@author: sayantanroy
"""

import numpy as np
import pickle
import sys
import os
import os.path


def numeric(equation):
    if '+' in equation:
        y = equation.split('+')
        x = int(y[0])+int(y[1])
    return str(x)

def parameter_read_data(data_file):
    file = open(data_file, 'r')
    lines = file.readlines()
    file.close()

    System_parameters = {}
    b = 0
    line_counter = 0
    for i in range(2, len(lines)):
         line_counter = line_counter+1
         if "===========" in lines[i]:
                  b = 1
         if b == 1:
            break

         if lines[i] != '\n' and '#' not in lines[i] and '========' not in lines:
               line = lines[i].strip('\n').strip(' ').split(':')
    #           print(line[0],strip(' '))
    #           print(line[1].strip(' '))
               System_parameters[line[0].strip(' ')] = line[1].strip(' ') 

#    print("=====================================================================================")
    return System_parameters, line_counter




def r_point_upper_half_grid(N):
    r_rel_x = []
    r_rel_y = []
    x_span = int((int(N))/2)+1
    r_pair_1 = 0
    for i in range(x_span):
        for j in range(x_span):
            if(i<=j):
              r_rel_x.append(i)
              r_rel_y.append(j)
              r_pair_1 = r_pair_1+1
    R_relative_1 = np.stack((r_rel_x,r_rel_y),axis = 1)
 
    return R_relative_1,r_pair_1


def r_point_full_grid(N):
    r_rel_x = []
    r_rel_y = []
    x_span = int((int(N))/2)+1
    r_pair_2 = 0
    for i in range(x_span):
        for j in range(x_span):
              r_rel_x.append(i)
              r_rel_y.append(j)
              r_pair_2 = r_pair_2+1
    R_relative_2 = np.stack((r_rel_x,r_rel_y),axis = 1)
    print("no of r points", r_pair_2)

    return R_relative_2,r_pair_2


def k_point_upper_half_grid(N):
    Kx = []
    Ky = []
    x_span = int((int(N))/2)+1
    k_pair = 0
    for i in range(x_span):
        for j in range(x_span):
            if(i<=j):
              Kx.append(i)
              Ky.append(j)
              k_pair = k_pair+1
    BZ = np.stack((Kx,Ky),axis = 1)

    return BZ,k_pair


def thermodynamic_measurement_read_data(Text_dir_eqm,equal_data_file,N,U,Mu,dtau,L,Realization):
    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()

    System_measurement = {}
    a = 0
    b = 0
    for i in range(1, len(lines)):

         if 'Green' in lines[i] and 'function' in lines[i]:
               b = 1  
       
         if b == 1:
            break

         #  Parameter  = "".join[line[0].split()]
         
         if a == 1:
              if lines[i] != '\n' and '#' not in lines[i] and '========' not in lines[i] and '**please note average energy and kinetic energy**' not in lines[i] and '**  include mu terms, both uniform and random  **' not in lines[i]:
                   line = lines[i].strip('\n').strip(' ').split('=')
                   print(line)
                   print(line[0].strip(' '))
                   print(line[1].strip(' '))
                   System_measurement[line[0].strip(' ')] = line[1].strip(' ')
       
         if 'Acceptance' in lines[i] and 'ratio' in lines[i]:
               a = 1
         


    filename_measurements = '%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eqm,N,U,Mu,dtau,L,Realization)
    data_measurements = System_measurement
    with open(filename_measurements, 'wb') as outfile:
        pickle.dump(data_measurements, outfile, pickle.HIGHEST_PROTOCOL)




def equal_time_density_density_correlation_function_real_space_read_data(Text_dir_eq_den_den_r,equal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs_1):


    Equal_Den_up_Den_up_r = np.zeros(r_rel_pairs_1)
    Equal_Den_up_Den_up_r_std = np.zeros(r_rel_pairs_1)

    Equal_Den_up_Den_down_r = np.zeros(r_rel_pairs_1)
    Equal_Den_up_Den_down_r_std = np.zeros(r_rel_pairs_1)

    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()
    
    a = 0
    b = 0

    eq_den_den_corr_index = 0
    r_index_max = r_rel_pairs_1-1
    for i in range(0, len(lines)):

         if  eq_den_den_corr_index>r_index_max:
             b = 1

         if b == 1:
             break

         if "density-density" in lines[i] and "correlation" in lines[i] and 'fn' in lines[i] and "up-up"  in lines[i] and "up-dn" in lines[i] and ':' in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i] and ':' not in lines[i]:

                line = lines[i].strip('\n').split()              
                if(len(line)==8):
                                     
                   
                     Equal_Den_up_Den_up_r[eq_den_den_corr_index]=float(line[2])
                     Equal_Den_up_Den_up_r_std[eq_den_den_corr_index] = float(line[4])
                     Equal_Den_up_Den_down_r[eq_den_den_corr_index]=float(line[5])
                     Equal_Den_up_Den_down_r_std[eq_den_den_corr_index] = float(line[7])
                     
                     eq_den_den_corr_index = eq_den_den_corr_index+1
 
    filename_eq_dudu_r = '%s/Equal_time_Density_up_Density_up_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
    data_eq_dudu_r = Equal_Den_up_Den_up_r
    with open(filename_eq_dudu_r, 'wb') as outfile:
       pickle.dump(data_eq_dudu_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_dudu_r_variation = '%s/Equal_time_Density_up_Density_up_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
    data_eq_dudu_r_variation = Equal_Den_up_Den_up_r_std
    with open(filename_eq_dudu_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_dudu_r_variation, outfile, pickle.HIGHEST_PROTOCOL)


    filename_eq_dudn_r = '%s/Equal_time_Density_up_Density_down_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
    data_eq_dudn_r = Equal_Den_up_Den_down_r
    with open(filename_eq_dudn_r, 'wb') as outfile:
       pickle.dump(data_eq_dudn_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_dudn_r_variation = '%s/Equal_time_Density_up_Density_down_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
    data_eq_dudn_r_variation = Equal_Den_up_Den_down_r_std
    with open(filename_eq_dudn_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_dudn_r_variation, outfile, pickle.HIGHEST_PROTOCOL)

       

def equal_time_doublon_doublon_correlation_function_real_space_read_data(Text_dir_eq_dbn_dbn_r,equal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs_2):


    Equal_Dbn_Dbn_r = np.zeros(r_rel_pairs_2)
    Equal_Dbn_Dbn_r_std = np.zeros(r_rel_pairs_2)

    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()
    
    a = 0
    b = 0

    eq_dbn_dbn_corr_index = 0
    r_index_max = r_rel_pairs_2-1
    for i in range(0, len(lines)):

         if  eq_dbn_dbn_corr_index>r_index_max:
             b = 1

         if b == 1:
             break

         if "nud-nud" in lines[i] and "correlation" in lines[i] and 'function' in lines[i] and ":"  in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i] and ':' not in lines[i]:

                line = lines[i].strip('\n').split()              
                if(len(line)==4):
                                     
                   
                     Equal_Dbn_Dbn_r[eq_dbn_dbn_corr_index]=float(line[2])
                     Equal_Dbn_Dbn_r_std[eq_dbn_dbn_corr_index] = float(line[3])
                   
                     eq_dbn_dbn_corr_index = eq_dbn_dbn_corr_index+1
                     
    filename_eq_dbdb_r = '%s/Equal_time_Doublon_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L,Realization)
    data_eq_dbdb_r = Equal_Dbn_Dbn_r
    with open(filename_eq_dbdb_r, 'wb') as outfile:
       pickle.dump(data_eq_dbdb_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_dbdb_r_variation = '%s/Equal_time_Doublon_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L,Realization)
    data_eq_dbdb_r_variation = Equal_Dbn_Dbn_r_std
    with open(filename_eq_dbdb_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_dbdb_r_variation, outfile, pickle.HIGHEST_PROTOCOL)



def equal_time_moment_moment_correlation_function_real_space_read_data(Text_dir_eq_m2_m2_r,equal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs_2):


    Equal_M2_M2_r = np.zeros(r_rel_pairs_2)
    Equal_M2_M2_r_std = np.zeros(r_rel_pairs_2)

    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()
    
    a = 0
    b = 0

    eq_m2_m2_corr_index = 0
    r_index_max = r_rel_pairs_2-1
    for i in range(0, len(lines)):

         if  eq_m2_m2_corr_index>r_index_max:
             b = 1

         if b == 1:
             break

         if "mi2x-mi2x" in lines[i] and "correlation" in lines[i] and 'function' in lines[i] and ":"  in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i] and ':' not in lines[i]:

                line = lines[i].strip('\n').split()              
                if(len(line)==4):
                                     
                   
                     Equal_M2_M2_r[eq_m2_m2_corr_index]=float(line[2])
                     Equal_M2_M2_r_std[eq_m2_m2_corr_index] = float(line[3])
                   
                     eq_m2_m2_corr_index = eq_m2_m2_corr_index+1
                     
    filename_eq_m2m2_r = '%s/Equal_time_Moment_Moment_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L,Realization)
    data_eq_m2m2_r = Equal_M2_M2_r
    with open(filename_eq_m2m2_r, 'wb') as outfile:
       pickle.dump(data_eq_m2m2_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_m2m2_r_variation = '%s/Equal_time_Moment_Moment_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L,Realization)
    data_eq_m2m2_r_variation = Equal_M2_M2_r_std
    with open(filename_eq_m2m2_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_m2m2_r_variation, outfile, pickle.HIGHEST_PROTOCOL)





def equal_time_density_up_doublon_correlation_function_real_space_read_data(Text_dir_eq_den_dbn_r,equal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs_2):


    Equal_Den_up_Dbn_r = np.zeros(r_rel_pairs_2)
    Equal_Den_up_Dbn_r_std = np.zeros(r_rel_pairs_2)

    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()
    
    a = 0
    b = 0

    eq_den_up_dbn_corr_index = 0
    r_index_max = r_rel_pairs_2-1
    for i in range(0, len(lines)):

         if  eq_den_up_dbn_corr_index>r_index_max:
             b = 1

         if b == 1:
             break

         if "nup-nud" in lines[i] and "correlation" in lines[i] and 'function' in lines[i] and ":"  in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i] and ':' not in lines[i]:

                line = lines[i].strip('\n').split()              
                if(len(line)==4):
                                     
                   
                     Equal_Den_up_Dbn_r[eq_den_up_dbn_corr_index]=float(line[2])
                     Equal_Den_up_Dbn_r_std[eq_den_up_dbn_corr_index] = float(line[3])
                   
                     eq_den_up_dbn_corr_index = eq_den_up_dbn_corr_index+1
                     

    filename_eq_dudb_r = '%s/Equal_time_Density_up_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)
    data_eq_dudb_r = Equal_Den_up_Dbn_r
    with open(filename_eq_dudb_r, 'wb') as outfile:
       pickle.dump(data_eq_dudb_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_dudb_r_variation = '%s/Equal_time_Density_up_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)
    data_eq_dudb_r_variation = Equal_Den_up_Dbn_r_std
    with open(filename_eq_dudb_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_dudb_r_variation, outfile, pickle.HIGHEST_PROTOCOL)
       


       
def equal_time_density_down_doublon_correlation_function_real_space_read_data(Text_dir_eq_den_dbn_r,equal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs_2):


    Equal_Den_down_Dbn_r = np.zeros(r_rel_pairs_2)
    Equal_Den_down_Dbn_r_std = np.zeros(r_rel_pairs_2)

    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()
    
    a = 0
    b = 0

    eq_den_down_dbn_corr_index = 0
    r_index_max = r_rel_pairs_2-1
    for i in range(0, len(lines)):

         if  eq_den_down_dbn_corr_index>r_index_max:
             b = 1

         if b == 1:
             break

         if "nd-nud" in lines[i] and "correlation" in lines[i] and 'function' in lines[i] and ":"  in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i] and ':' not in lines[i]:

                line = lines[i].strip('\n').split()              
                if(len(line)==4):
                                     
                   
                     Equal_Den_down_Dbn_r[eq_den_down_dbn_corr_index]=float(line[2])
                     Equal_Den_down_Dbn_r_std[eq_den_down_dbn_corr_index] = float(line[3])
                   
                     eq_den_down_dbn_corr_index = eq_den_down_dbn_corr_index+1
                     

    filename_eq_dndb_r = '%s/Equal_time_Density_down_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)
    data_eq_dndb_r = Equal_Den_down_Dbn_r
    with open(filename_eq_dndb_r, 'wb') as outfile:
       pickle.dump(data_eq_dndb_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_dndb_r_variation = '%s/Equal_time_Density_down_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)
    data_eq_dndb_r_variation = Equal_Den_down_Dbn_r_std
    with open(filename_eq_dndb_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_dndb_r_variation, outfile, pickle.HIGHEST_PROTOCOL)


def equal_time_spin_spin_xx_correlation_function_real_space_read_data(Text_dir_eq_spin_spin_r,equal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs_2):


    Equal_SxSx_r = np.zeros(r_rel_pairs_2)
    Equal_SxSx_r_std = np.zeros(r_rel_pairs_2)
    
    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()
    
    a = 0
    b = 0

    eq_sxsx_corr_index = 0
    r_index_max = r_rel_pairs_2-1
    for i in range(0, len(lines)):

         if  eq_sxsx_corr_index>r_index_max:
             b = 1

         if b == 1:
             break

         if "xx" in lines[i] and "Spin" in lines[i] and "correlation" in lines[i] and 'function' in lines[i] and ":"  in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i] and ':' not in lines[i]:

                line = lines[i].strip('\n').split()              
                if(len(line)==4):
                                     
                   
                     Equal_SxSx_r[eq_sxsx_corr_index]=float(line[2])
                     Equal_SxSx_r_std[eq_sxsx_corr_index] = float(line[3])
                   
                     eq_sxsx_corr_index = eq_sxsx_corr_index+1
                     

    filename_eq_sxsx_r = '%s/Equal_time_SxSx_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)
    data_eq_sxsx_r = Equal_SxSx_r
    with open(filename_eq_sxsx_r, 'wb') as outfile:
       pickle.dump(data_eq_sxsx_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_sxsx_r_variation = '%s/Equal_time_SxSx_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)
    data_eq_sxsx_r_variation = Equal_SxSx_r_std
    with open(filename_eq_sxsx_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_sxsx_r_variation, outfile, pickle.HIGHEST_PROTOCOL)


def equal_time_spin_spin_zz_correlation_function_real_space_read_data(Text_dir_eq_spin_spin_r,equal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs_2):


    Equal_SzSz_r = np.zeros(r_rel_pairs_2)
    Equal_SzSz_r_std = np.zeros(r_rel_pairs_2)

    file = open(equal_data_file, 'r')
    lines = file.readlines()
    file.close()
    
    a = 0
    b = 0

    eq_szsz_corr_index = 0
    r_index_max = r_rel_pairs_2-1
    for i in range(0, len(lines)):

         if  eq_szsz_corr_index>r_index_max:
             b = 1

         if b == 1:
             break

         if "zz" in lines[i] and "Spin" in lines[i] and "correlation" in lines[i] and 'function' in lines[i] and ":"  in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i] and ':' not in lines[i]:

                line = lines[i].strip('\n').split()              
                if(len(line)==4):
                                     
                   
                     Equal_SzSz_r[eq_szsz_corr_index]=float(line[2])
                     Equal_SzSz_r_std[eq_szsz_corr_index] = float(line[3])
                   
                     eq_szsz_corr_index = eq_szsz_corr_index+1
                     

    filename_eq_szsz_r = '%s/Equal_time_SzSz_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)
    data_eq_szsz_r = Equal_SzSz_r
    with open(filename_eq_szsz_r, 'wb') as outfile:
       pickle.dump(data_eq_szsz_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_eq_szsz_r_variation = '%s/Equal_time_SzSz_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)
    data_eq_szsz_r_variation = Equal_SzSz_r_std
    with open(filename_eq_szsz_r_variation, 'wb') as outfile:
       pickle.dump(data_eq_szsz_r_variation, outfile, pickle.HIGHEST_PROTOCOL)

       


       
       
def main(total,cmdargs):
    if(total!=7):
        raise ValueError('missing args')

    N = cmdargs[1]
    U = cmdargs[2]
    mu = cmdargs[3]
    L = cmdargs[4]
    Dtau = cmdargs[5]
    run_no=int(cmdargs[6])

    print(U,"U")
    print(L,"L")
    Beta = str((float(L))*(float(Dtau)))
    print(Beta,"Beta")


    Text_dir_params = '../../Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Parameters'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_params):
        os.makedirs(Text_dir_params)

    Text_dir_eqm = '../../Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eqm):
        os.makedirs(Text_dir_eqm)
          
    Text_dir_eq_den_den_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_density_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_den_den_r):
        os.makedirs(Text_dir_eq_den_den_r)

    Text_dir_eq_dbn_dbn_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Doublon_doublon_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_dbn_dbn_r):
        os.makedirs(Text_dir_eq_dbn_dbn_r)

    Text_dir_eq_den_dbn_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_doublon_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_den_dbn_r):
        os.makedirs(Text_dir_eq_den_dbn_r)

    Text_dir_eq_m2_m2_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Moment_moment_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_m2_m2_r):
        os.makedirs(Text_dir_eq_m2_m2_r)

    Text_dir_eq_spin_spin_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_spin_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_spin_spin_r):
        os.makedirs(Text_dir_eq_spin_spin_r)


    bz_grid,bz_pair = k_point_upper_half_grid(N)
    r_rel_grid_1,r_rel_pairs_1 = r_point_upper_half_grid(N)
    r_rel_grid_2,r_rel_pairs_2 = r_point_full_grid(N)

        

    run_counter = 0
    while(run_counter<run_no):
         
         realization = str(run_counter)
         run_counter=run_counter+1 
         equal_data_file = '../../N_%s/U_%s/Mu_%s/dtau_%s_L_%s/Realization_%s/rz%sl%su%sdt%smu%sr%s.out' %(N,U,mu,Dtau,L,realization,N,L,U,Dtau,mu,realization)

         if os.path.exists(equal_data_file):

                thermodynamic_measurement_read_data(Text_dir_eqm,equal_data_file,N,U,mu,Dtau,L,realization)
            
                equal_time_density_density_correlation_function_real_space_read_data(Text_dir_eq_den_den_r,equal_data_file,N,U,mu,Dtau,L,realization,r_rel_pairs_1)
                equal_time_doublon_doublon_correlation_function_real_space_read_data(Text_dir_eq_dbn_dbn_r,equal_data_file,N,U,mu,Dtau,L,realization,r_rel_pairs_2)
                equal_time_moment_moment_correlation_function_real_space_read_data(Text_dir_eq_m2_m2_r,equal_data_file,N,U,mu,Dtau,L,realization,r_rel_pairs_2)

                equal_time_density_up_doublon_correlation_function_real_space_read_data(Text_dir_eq_den_dbn_r,equal_data_file,N,U,mu,Dtau,L,realization,r_rel_pairs_2)
                equal_time_density_down_doublon_correlation_function_real_space_read_data(Text_dir_eq_den_dbn_r,equal_data_file,N,U,mu,Dtau,L,realization,r_rel_pairs_2)

                equal_time_spin_spin_xx_correlation_function_real_space_read_data(Text_dir_eq_spin_spin_r,equal_data_file,N,U,mu,Dtau,L,realization,r_rel_pairs_2)
                equal_time_spin_spin_zz_correlation_function_real_space_read_data(Text_dir_eq_spin_spin_r,equal_data_file,N,U,mu,Dtau,L,realization,r_rel_pairs_2)

         else:
                print("File not found",mu,run_counter)
                continue




if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)

