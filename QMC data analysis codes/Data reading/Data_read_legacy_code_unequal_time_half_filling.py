#/usr/bin/env python3
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




def r_point_upper_half_grid(Text_dir,N):

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
    filename_r_grid = '%s/Relative_seperation_co-oridinates_set_1_N_%s.pkl' %(Text_dir,N)
    data_r_grid = R_relative_1
    with open(filename_r_grid, 'wb') as outfile:
       pickle.dump(data_r_grid, outfile, pickle.HIGHEST_PROTOCOL)
    return R_relative_1,r_pair_1


def r_point_full_grid(Text_dir,N):

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
    filename_r_grid = '%s/Relative_seperation_co-oridinates_set_2_N_%s.pkl' %(Text_dir,N)
    data_r_grid = R_relative_2
    with open(filename_r_grid, 'wb') as outfile:
       pickle.dump(data_r_grid, outfile, pickle.HIGHEST_PROTOCOL)
    return R_relative_2,r_pair_2


def k_point_grid_upper_half_bz(Text_dir,N):
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
    print("no of k points", k_pair)
    print(BZ)
    filename_k_grid = '%s/Upper_half_brillouin_zone_co-oridinates_N_%s.pkl' %(Text_dir,N)
    data_k_grid = BZ
    with open(filename_k_grid, 'wb') as outfile:
       pickle.dump(data_k_grid, outfile, pickle.HIGHEST_PROTOCOL)
    return BZ,k_pair


def k_point_grid_full_bz(Text_dir,N):
    Kx = []
    Ky = []
    x_span = int((int(N))/2)+1
    k_pair = 0
    for i in range(x_span):
        for j in range(x_span):
              Kx.append(i)
              Ky.append(j)
              k_pair = k_pair+1
    BZ = np.stack((Kx,Ky),axis = 1)
    print("no of k points", k_pair)
    print(BZ)
    filename_k_grid = '%s/Brillouin_zone_co-oridinates_N_%s.pkl' %(Text_dir,N)
    data_k_grid = BZ
    with open(filename_k_grid, 'wb') as outfile:
       pickle.dump(data_k_grid, outfile, pickle.HIGHEST_PROTOCOL)
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


       

def unequal_time_green_function_real_space_read_data(Text_dir_gf_r,unequal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs):
    
    print("G_r_working")
    G_r = np.zeros((int(L)+1,r_rel_pairs))
    G_r_std = np.zeros((int(L)+1,r_rel_pairs))
    
   
#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    r_index = 0
    r_index_max = r_rel_pairs
    
    for i in range(0, len(lines)):

         if  "G" in lines[i] and "qx" in lines[i] and "qy" in lines[i] and "ti" in lines[i] and ":":
             b = 1

         if b == 1:
             break


         if "G" in lines[i] and "nx" in lines[i] and 'ny' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "nx" in lines[i] and "=" in lines[i] and "ny" in lines[i] and "=" in lines[i]:

                             r_index = r_index+1


                   if "tau" not in lines[i] and "nx" not in lines[i] and "=" not in lines[i] and "ny" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("r_index",r_index)
                             #print(line)   
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     G_r[L_counter][r_index-1]=float(line[1])
                                     G_r_std[L_counter][r_index-1] = float(line[3])

    filename_g_r = '%s/Retarded_Green_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L,Realization)
    data_g_r = G_r
    with open(filename_g_r, 'wb') as outfile:
       pickle.dump(data_g_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_g_r_variation = '%s/Retarded_Green_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L,Realization)
    data_g_r_variation = G_r_std
    with open(filename_g_r_variation, 'wb') as outfile:
       pickle.dump(data_g_r_variation, outfile, pickle.HIGHEST_PROTOCOL)

    filename_g_r_text = '%s/Retarded_Green_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_gf_r,N,U,Mu,dtau,L,Realization)
    filename_g_r_variation_text = '%s/Retarded_Green_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_gf_r,N,U,Mu,dtau,L,Realization)

    np.savetxt(filename_g_r_text,G_r)
    np.savetxt(filename_g_r_variation_text,G_r_std)



def unequal_time_green_function_momentum_space_read_data(Text_dir_gf_k,unequal_data_file,N,U,Mu,dtau,L,Realization,bz_pairs):
    
    print("G_k_working")
    G_k = np.zeros((int(L)+1,bz_pairs))
    G_k_std = np.zeros((int(L)+1,bz_pairs))
    
   
#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    k_index = 0
    k_index_max = bz_pairs
    
    for i in range(0, len(lines)):
         
         if  "G" in lines[i] and "nx" in lines[i] and "ny" in lines[i] and "omega" in lines[i] and "2" in lines[i] and "pi" in lines[i] and "T" in lines[i] :
             b = 1
        
         if b == 1:
             break         


         if "G" in lines[i] and "qx" in lines[i] and 'qy' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:
                   
                  
                   if "qx" in lines[i] and "=" in lines[i] and "qy" in lines[i] and "=" in lines[i]:
                    
                             k_index = k_index+1
                           
                   
                   if "tau" not in lines[i] and "qx" not in lines[i] and "=" not in lines[i] and "qy" not in lines[i] and "=" not in lines[i]: 
                            
                             line = lines[i].strip('\n').split()
                             print("k_index",k_index)
                             #print(line)                
                             if(len(line)==4):
                                     L_counter = int(line[0])                                     
                                     G_k[L_counter][k_index-1]=float(line[1])
                                     G_k_std[L_counter][k_index-1] = float(line[3])
      
    filename_g_k = '%s/Retarded_Green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L,Realization)
    data_g_k = G_k
    with open(filename_g_k, 'wb') as outfile:
       pickle.dump(data_g_k, outfile, pickle.HIGHEST_PROTOCOL)

    filename_g_k_variation = '%s/Retarded_Green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L,Realization)
    data_g_k_variation = G_k_std
    with open(filename_g_k_variation, 'wb') as outfile:
       pickle.dump(data_g_k_variation, outfile, pickle.HIGHEST_PROTOCOL)

    filename_g_k_text = '%s/Retarded_Green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_gf_k,N,U,Mu,dtau,L,Realization)
    filename_g_k_variation_text = '%s/Retarded_Green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_gf_k,N,U,Mu,dtau,L,Realization)

    np.savetxt(filename_g_k_text,G_k)
    np.savetxt(filename_g_k_variation_text,G_k_std)



def unequal_time_density_correlation_function_real_space_read_data(Text_dir_den_corr_r,unequal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs):

    print("Den_r_working")
    Den_r = np.zeros((int(L)+1,r_rel_pairs))
    Den_r_std = np.zeros((int(L)+1,r_rel_pairs))


#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    r_index = 0
    r_index_max = r_rel_pairs

    for i in range(0, len(lines)):

         if "chi_zz" in lines[i] and "nx" in lines[i] and "ny" in lines[i] and "ti" in lines[i] and ":" in lines[i]:
             b = 1

         if b == 1:
             break


         if "den" in lines[i] and "nx" in lines[i] and 'ny' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "nx" in lines[i] and "=" in lines[i] and "ny" in lines[i] and "=" in lines[i]:

                             r_index = r_index+1


                   if "tau" not in lines[i] and "nx" not in lines[i] and "=" not in lines[i] and "ny" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("r_index",r_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Den_r[L_counter][r_index-1]=float(line[1])
                                     Den_r_std[L_counter][r_index-1] = float(line[3])

    filename_den_r = '%s/Density_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L,Realization)
    data_den_r = Den_r
    with open(filename_den_r, 'wb') as outfile:
       pickle.dump(data_den_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_den_r_variation = '%s/Density_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L,Realization)
    data_den_r_variation = Den_r_std
    with open(filename_den_r_variation, 'wb') as outfile:
       pickle.dump(data_den_r_variation, outfile, pickle.HIGHEST_PROTOCOL)


def unequal_time_density_correlation_function_momentum_space_read_data(Text_dir_den_corr_k,unequal_data_file,N,U,Mu,dtau,L,Realization,bz_pairs):

    print("Den_k_working")
    Den_k = np.zeros((int(L)+1,bz_pairs))
    Den_k_std = np.zeros((int(L)+1,bz_pairs))


#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    k_index = 0
    k_index_max = bz_pairs

    for i in range(0, len(lines)):

         if  "chi_xx" in lines[i] and "qx" in lines[i] and "qy" in lines[i] and "ti" in lines[i]:
             b = 1

         if b == 1:
             break


         if "den" in lines[i] and "qx" in lines[i] and 'qy' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "qx" in lines[i] and "=" in lines[i] and "qy" in lines[i] and "=" in lines[i]:

                             k_index = k_index+1


                   if "tau" not in lines[i] and "qx" not in lines[i] and "=" not in lines[i] and "qy" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("k_index",k_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Den_k[L_counter][k_index-1]=float(line[1])
                                     Den_k_std[L_counter][k_index-1] = float(line[3])

    filename_den_k = '%s/Density_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L,Realization)
    data_den_k = Den_k
    with open(filename_den_k, 'wb') as outfile:
       pickle.dump(data_den_k, outfile, pickle.HIGHEST_PROTOCOL)

    filename_den_k_variation = '%s/Density_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L,Realization)
    data_den_k_variation = Den_k_std
    with open(filename_den_k_variation, 'wb') as outfile:
       pickle.dump(data_den_k_variation, outfile, pickle.HIGHEST_PROTOCOL)


def unequal_time_chi_xx_correlation_function_real_space_read_data(Text_dir_chi_xx_corr_r,unequal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs):

    print("Chi_xx_r_working")

    Chi_xx_r = np.zeros((int(L)+1,r_rel_pairs))
    Chi_xx_r_std = np.zeros((int(L)+1,r_rel_pairs))


#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    r_index = 0
    r_index_max = r_rel_pairs

    for i in range(0, len(lines)):

         if "den" in lines[i] and "nx" in lines[i] and "ny" in lines[i] and "ti" in lines[i] and ":" in lines[i]:
             b = 1

         if b == 1:
             break


         if "chi_xx" in lines[i] and "nx" in lines[i] and 'ny' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "nx" in lines[i] and "=" in lines[i] and "ny" in lines[i] and "=" in lines[i]:

                             r_index = r_index+1


                   if "tau" not in lines[i] and "nx" not in lines[i] and "=" not in lines[i] and "ny" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("r_index",r_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Chi_xx_r[L_counter][r_index-1]=float(line[1])
                                     Chi_xx_r_std[L_counter][r_index-1] = float(line[3])

    filename_chi_xx_r = '%s/Spin_xx_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L,Realization)
    data_chi_xx_r = Chi_xx_r
    with open(filename_chi_xx_r, 'wb') as outfile:
       pickle.dump(data_chi_xx_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_chi_xx_r_variation = '%s/Spin_xx_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L,Realization)
    data_chi_xx_r_variation = Chi_xx_r_std
    with open(filename_chi_xx_r_variation, 'wb') as outfile:
       pickle.dump(data_chi_xx_r_variation, outfile, pickle.HIGHEST_PROTOCOL)



def unequal_time_chi_xx_correlation_function_momentum_space_read_data(Text_dir_chi_xx_corr_k,unequal_data_file,N,U,Mu,dtau,L,Realization,bz_pairs):

    print("Chi_xx_k_working")

    Chi_xx_k = np.zeros((int(L)+1,bz_pairs))
    Chi_xx_k_std = np.zeros((int(L)+1,bz_pairs))
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    k_index = 0
    k_index_max = bz_pairs

    for i in range(0, len(lines)):

         if "chi_zz" in lines[i] and "qx" in lines[i] and "qy" in lines[i] and "ti" in lines[i]:
             b = 1

         if b == 1:
             break

         if "chi_xx" in lines[i] and "qx" in lines[i] and 'qy' in lines[i] and "ti" in lines[i] and ":" in lines[i]:
                a = 1


         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "qx" in lines[i] and "=" in lines[i] and "qy" in lines[i] and "=" in lines[i]:

                             k_index = k_index+1


                   if "tau" not in lines[i] and "qx" not in lines[i] and "=" not in lines[i] and "qy" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("k_index",k_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Chi_xx_k[L_counter][k_index-1]=float(line[1])
                                     Chi_xx_k_std[L_counter][k_index-1] = float(line[3])

    filename_chi_xx_k = '%s/Spin_xx_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L,Realization)
    data_chi_xx_k = Chi_xx_k
    with open(filename_chi_xx_k, 'wb') as outfile:
       pickle.dump(data_chi_xx_k, outfile, pickle.HIGHEST_PROTOCOL)

    filename_chi_xx_k_variation = '%s/Spin_xx_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L,Realization)
    data_chi_xx_k_variation = Chi_xx_k_std
    with open(filename_chi_xx_k_variation, 'wb') as outfile:
       pickle.dump(data_chi_xx_k_variation, outfile, pickle.HIGHEST_PROTOCOL)


def unequal_time_chi_zz_correlation_function_real_space_read_data(Text_dir_chi_zz_corr_r,unequal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs):

    print("chi_zz_r_working")

    Chi_zz_r = np.zeros((int(L)+1,r_rel_pairs))
    Chi_zz_r_std = np.zeros((int(L)+1,r_rel_pairs))


#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    r_index = 0
    r_index_max = r_rel_pairs

    for i in range(0, len(lines)):

         if "den" in lines[i] and "qx" in lines[i] and "qy" in lines[i] and "ti" in lines[i] and ":" in lines[i]:
             b = 1

         if b == 1:
             break


         if "chi_zz" in lines[i] and "nx" in lines[i] and 'ny' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "nx" in lines[i] and "=" in lines[i] and "ny" in lines[i] and "=" in lines[i]:

                             r_index = r_index+1


                   if "tau" not in lines[i] and "nx" not in lines[i] and "=" not in lines[i] and "ny" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("r_index",r_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Chi_zz_r[L_counter][r_index-1]=float(line[1])
                                     Chi_zz_r_std[L_counter][r_index-1] = float(line[3])

    filename_chi_zz_r = '%s/Spin_zz_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L,Realization)
    data_chi_zz_r = Chi_zz_r
    with open(filename_chi_zz_r, 'wb') as outfile:
       pickle.dump(data_chi_zz_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_chi_zz_r_variation = '%s/Spin_zz_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L,Realization)
    data_chi_zz_r_variation = Chi_zz_r_std
    with open(filename_chi_zz_r_variation, 'wb') as outfile:
       pickle.dump(data_chi_zz_r_variation, outfile, pickle.HIGHEST_PROTOCOL)

    filename_chi_zz_r_text = '%s/Spin_zz_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L,Realization)
    filename_chi_zz_r_variation_text = '%s/Spin_zz_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L,Realization)

    np.savetxt(filename_chi_zz_r_text,Chi_zz_r)
    np.savetxt(filename_chi_zz_r_variation_text,Chi_zz_r_std)

    

def unequal_time_chi_zz_correlation_function_momentum_space_read_data(Text_dir_chi_zz_corr_k,unequal_data_file,N,U,Mu,dtau,L,Realization,bz_pairs):

    print("Chi_zz_k_working")
    Chi_zz_k = np.zeros((int(L)+1,bz_pairs))
    Chi_zz_k_std = np.zeros((int(L)+1,bz_pairs))
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    k_index = 0
    k_index_max = bz_pairs

    for i in range(0, len(lines)):

         if "chi,z" in lines[i] and ":" in lines[i]:
             b = 1

         if b == 1:
             break


         if "chi_zz" in lines[i] and "qx" in lines[i] and 'qy' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "qx" in lines[i] and "=" in lines[i] and "qy" in lines[i] and "=" in lines[i]:

                             k_index = k_index+1


                   if "tau" not in lines[i] and "qx" not in lines[i] and "=" not in lines[i] and "qy" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("k_index",k_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Chi_zz_k[L_counter][k_index-1]=float(line[1])
                                     Chi_zz_k_std[L_counter][k_index-1] = float(line[3])

    filename_chi_zz_k = '%s/Spin_zz_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L,Realization)
    data_chi_zz_k = Chi_zz_k
    with open(filename_chi_zz_k, 'wb') as outfile:
       pickle.dump(data_chi_zz_k, outfile, pickle.HIGHEST_PROTOCOL)

    filename_chi_zz_k_variation = '%s/Spin_zz_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L,Realization)
    data_chi_zz_k_variation = Chi_zz_k_std
    with open(filename_chi_zz_k_variation, 'wb') as outfile:
       pickle.dump(data_chi_zz_k_variation, outfile, pickle.HIGHEST_PROTOCOL)

    filename_chi_zz_k_text = '%s/Spin_zz_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L,Realization)
    filename_chi_zz_k_variation_text = '%s/Spin_zz_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L,Realization)

    np.savetxt(filename_chi_zz_k_text,Chi_zz_k)
    np.savetxt(filename_chi_zz_k_variation_text,Chi_zz_k_std)


def unequal_time_current_correlator_real_space_read_data(Text_dir_cc_r,unequal_data_file,N,U,Mu,dtau,L,Realization,r_rel_pairs):


    Curr_r = np.zeros((int(L)+1,r_rel_pairs))
    Curr_r_std = np.zeros((int(L)+1,r_rel_pairs))


#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    r_index = 0
    r_index_max = r_rel_pairs

    for i in range(0, len(lines)):

         if "current" in lines[i] and "qx" in lines[i] and "qy" in lines[i] and "ti" in lines[i] and ":" in lines[i] :
             b = 1

         if b == 1:
             break

         if "current" in lines[i] and "nx" in lines[i] and 'ny' in lines[i] and "ti" in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "nx" in lines[i] and "=" in lines[i] and "ny" in lines[i] and "=" in lines[i]:

                             r_index = r_index+1


                   if "tau" not in lines[i] and "nx" not in lines[i] and "=" not in lines[i] and "ny" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("r_index",r_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Curr_r[L_counter][r_index-1]=float(line[1])
                                     Curr_r_std[L_counter][r_index-1] = float(line[3])

    filename_curr_r = '%s/Current_current_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_cc_r,N,U,Mu,dtau,L,Realization)
    data_curr_r = Curr_r
    with open(filename_curr_r, 'wb') as outfile:
       pickle.dump(data_curr_r, outfile, pickle.HIGHEST_PROTOCOL)

    filename_curr_r_variation = '%s/Current_current_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_cc_r,N,U,Mu,dtau,L,Realization)
    data_curr_r_variation = Curr_r_std
    with open(filename_curr_r_variation, 'wb') as outfile:
       pickle.dump(data_curr_r_variation, outfile, pickle.HIGHEST_PROTOCOL)

    filename_curr_r_text = '%s/Current_current_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_cc_r,N,U,Mu,dtau,L,Realization)
    filename_curr_r_variation_text = '%s/Current_current_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_cc_r,N,U,Mu,dtau,L,Realization)

    np.savetxt(filename_curr_r_text,Curr_r)
    np.savetxt(filename_curr_r_variation_text,Curr_r_std)

    

def unequal_time_current_correlator_momentum_space_read_data(Text_dir_cc_k,unequal_data_file,N,U,Mu,dtau,L,Realization,bz_pairs):


    Curr_k = np.zeros((int(L)+1,bz_pairs))
    Curr_k_std = np.zeros((int(L)+1,bz_pairs))


#    print(spin_species,"Spin_species")
    file = open(unequal_data_file, 'r')
    lines = file.readlines()
    file.close()
    a = 0
    b = 0
    k_index = 0
    k_index_max = bz_pairs

    for i in range(0, len(lines)):

         if  "current" in lines[i] and "nx" in lines[i] and "ny" in lines[i] and "omega" in lines[i] and "2" in lines[i] and "pi" in lines[i] and "T" in lines[i] :
             b = 1

         if b == 1:
             break

         if "current" in lines[i] and "qx" in lines[i] and 'qy' in lines[i] and "ti"  in lines[i] and ":" in lines[i]:
                a = 1

         if a == 1:
             if lines[i] != '\n' and '#' not in lines[i]:


                   if "qx" in lines[i] and "=" in lines[i] and "qy" in lines[i] and "=" in lines[i]:

                             k_index = k_index+1


                   if "tau" not in lines[i] and "qx" not in lines[i] and "=" not in lines[i] and "qy" not in lines[i] and "=" not in lines[i]:

                             line = lines[i].strip('\n').split()
                             print("k_index",k_index)
                             #print(line)
                             if(len(line)==4):
                                     L_counter = int(line[0])
                                     Curr_k[L_counter][k_index-1]=float(line[1])
                                     Curr_k_std[L_counter][k_index-1] = float(line[3])

    filename_curr_k = '%s/Current_current_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_cc_k,N,U,Mu,dtau,L,Realization)
    data_curr_k = Curr_k
    with open(filename_curr_k, 'wb') as outfile:
       pickle.dump(data_curr_k, outfile, pickle.HIGHEST_PROTOCOL)

    filename_curr_k_variation = '%s/Current_current_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_cc_k,N,U,Mu,dtau,L,Realization)
    data_curr_k_variation = Curr_k_std
    with open(filename_curr_k_variation, 'wb') as outfile:
       pickle.dump(data_curr_k_variation, outfile, pickle.HIGHEST_PROTOCOL)       
       
    filename_curr_k_text = '%s/Current_current_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_cc_k,N,U,Mu,dtau,L,Realization)
    filename_curr_k_variation_text = '%s/Current_current_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_cc_k,N,U,Mu,dtau,L,Realization)

    np.savetxt(filename_curr_k_text,Curr_k)
    np.savetxt(filename_curr_k_variation_text,Curr_k_std)
    

def main(total,cmdargs):
    if(total!=7):
        raise ValueError('missing args')

    N = cmdargs[1]
    U = cmdargs[2]
    #t_prime = cmdargs[3]
    mu = cmdargs[3]
    L = cmdargs[4]
    Dtau = cmdargs[5]
    run_no = cmdargs[6]
    print(U,"U")
    print(L,"L")
    Beta = str((float(L))*(float(Dtau)))
    print(Beta,"Beta")
    print(mu,"Mu")

    Text_dir_params = '../../Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Parameters'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_params):
          os.makedirs(Text_dir_params)

    Text_dir_eqm = '../../Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eqm):
          os.makedirs(Text_dir_eqm)
          
    Text_dir_gf_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_r):
         os.makedirs(Text_dir_gf_r)

    Text_dir_gf_k = '../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_k):
         os.makedirs(Text_dir_gf_k)


    Text_dir_den_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_den_r):
         os.makedirs(Text_dir_den_r)

    Text_dir_den_k = '../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_den_k):
         os.makedirs(Text_dir_den_k)

    Text_dir_chi_xx_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_xx_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_chi_xx_r):
         os.makedirs(Text_dir_chi_xx_r)

    Text_dir_chi_xx_k = '../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_xx_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_chi_xx_k):
         os.makedirs(Text_dir_chi_xx_k)

    Text_dir_chi_zz_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_zz_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_chi_zz_r):
         os.makedirs(Text_dir_chi_zz_r)

    Text_dir_chi_zz_k = '../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_zz_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_chi_zz_k):
         os.makedirs(Text_dir_chi_zz_k)

    Text_dir_cc_r = '../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_cc_r):
         os.makedirs(Text_dir_cc_r)

    Text_dir_cc_k = '../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_cc_k):
         os.makedirs(Text_dir_cc_k)


    bz_grid_1,bz_pair_number_1 = k_point_grid_upper_half_bz(Text_dir_params,N)
    bz_grid_2,bz_pair_number_2 = k_point_grid_full_bz(Text_dir_params,N)

    r_rel_grid_1,r_rel_pair_number_1 = r_point_upper_half_grid(Text_dir_params,N)
    r_rel_grid_2,r_rel_pair_number_2 = r_point_full_grid(Text_dir_params,N)

        

    run_counter = 0
    while(run_counter<int(run_no)):

        realization = str(run_counter)
          
        equal_data_file = "../../N_%s_half_filling/U_%s/Mu_%s/dtau_%s_L_%s/Realization_%s/rz%sl%su%sdt%smu%sr%s.out"%(N,U,mu,Dtau,L,realization,N,L,U,Dtau,mu,realization)
        unequal_data_file = "../../N_%s_half_filling/U_%s/Mu_%s/dtau_%s_L_%s/Realization_%s/gz%sl%su%sdt%smu%sr%s.out"%(N,U,mu,Dtau,L,realization,N,L,U,Dtau,mu,realization)
        run_counter=run_counter+1
        if os.path.exists(equal_data_file):
           
           thermodynamic_measurement_read_data(Text_dir_eqm,equal_data_file,N,U,mu,Dtau,L,realization)

           unequal_time_green_function_real_space_read_data(Text_dir_gf_r,unequal_data_file,N,U,mu,Dtau,L,realization,r_rel_pair_number_1)
           unequal_time_green_function_momentum_space_read_data(Text_dir_gf_k,unequal_data_file,N,U,mu,Dtau,L,realization,bz_pair_number_1)

           unequal_time_density_correlation_function_real_space_read_data(Text_dir_den_r,unequal_data_file,N,U,mu,Dtau,L,realization,r_rel_pair_number_1)
           unequal_time_density_correlation_function_momentum_space_read_data(Text_dir_den_k,unequal_data_file,N,U,mu,Dtau,L,realization,bz_pair_number_1)

           unequal_time_chi_xx_correlation_function_real_space_read_data(Text_dir_chi_xx_r,unequal_data_file,N,U,mu,Dtau,L,realization,r_rel_pair_number_1)
           unequal_time_chi_xx_correlation_function_momentum_space_read_data(Text_dir_chi_xx_k,unequal_data_file,N,U,mu,Dtau,L,realization,bz_pair_number_1)

           unequal_time_chi_zz_correlation_function_real_space_read_data(Text_dir_chi_zz_r,unequal_data_file,N,U,mu,Dtau,L,realization,r_rel_pair_number_1)
           unequal_time_chi_zz_correlation_function_momentum_space_read_data(Text_dir_chi_zz_k,unequal_data_file,N,U,mu,Dtau,L,realization,bz_pair_number_1)

           #unequal_time_current_correlator_real_space_read_data(Text_dir_cc_r,unequal_data_file,N,U,mu,Dtau,L,realization,r_rel_pair_number_2)
           #unequal_time_current_correlator_momentum_space_read_data(Text_dir_cc_k,unequal_data_file,N,U,mu,Dtau,L,realization,bz_pair_number_2) 


        else:
           print("File not found",mu,run_counter)
           continue
 

if __name__ == '__main__':
    sys.argv 
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)

