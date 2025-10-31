#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:11:44 2022

@author: roy.369
"""


import numpy as np
import pickle5 as pickle
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import griddata
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colorbar
from matplotlib import rc
mpl.rcParams['axes.linewidth'] = 5

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx


def r_point_grid_upper_half(N):

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


def k_point_grid_upper_half_bz(N):
    K_label = []
    Kx = []
    Ky = []
    x_span = int((int(N))/2)+1
    k_pair = 0
    for i in range(x_span):
        for j in range(x_span):
            if(j>=i):
              K_label.append(k_pair)  
              Kx.append(i)
              Ky.append(j)
              k_pair = k_pair+1

              
    BZ = np.stack((K_label,Kx,Ky),axis = 1)
    
    return BZ,k_pair


def k_point_grid_full_bz(N):
    K_label = []
    Kx = []
    Ky = []
    x_span = int((int(N))/2)+1
    k_pair = 0
    for i in range(x_span):
        for j in range(x_span):
              K_label.append(k_pair)
              Kx.append(i)
              Ky.append(j)
              k_pair = k_pair+1
    BZ = np.stack((K_label,Kx,Ky),axis = 1)

    return BZ,k_pair




def unequal_time_current_correlation_function_momentum_space_average(Text_dir_curr_curr_corr_k,Text_dir_curr_curr_corr_k_avg_ac,N,U,Mu,dtau,L,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_curr_k_0 = '%s/Current_current_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L)

   if os.path.exists(filename_curr_k_0):

       with open(filename_curr_k_0, 'rb') as infile:
           current_correlation_k_0 = pickle.load(infile)

       filename_current_k_variation_0 = '%s/Current_current_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L)
       with open(filename_current_k_variation_0, 'rb') as infile:
           current_correlation_k_std_0 = pickle.load(infile)


       Current_correlation_k = current_correlation_k_0.copy()
       Current_correlation_k_var = np.power(current_correlation_k_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            filename_curr_k = '%s/Current_current_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_curr_k):

               r_counter=r_counter+1
               with open(filename_curr_k, 'rb') as infile:
                   current_correlation_k = pickle.load(infile)

               filename_current_k_variation = '%s/Current_current_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L,realization)
               with open(filename_current_k_variation, 'rb') as infile:
                   current_correlation_k_std = pickle.load(infile)


               Current_correlation_X_k = np.add(Current_correlation_k,current_correlation_k)
               Current_correlation_k = Current_correlation_X_k.copy()

               Current_correlation_Y_k_var = np.power(current_correlation_k_std,2)
               Current_correlation_X_k_var = np.add(Current_correlation_k_var,Current_correlation_Y_k_var)
               Current_correlation_k_var = Current_correlation_X_k_var.copy()


       Current_correlation_k_avg = Current_correlation_k/r_counter

       Current_correlation_k_std_avg = np.sqrt(Current_correlation_k_var/(r_counter*r_counter))

       filename_current_k_avg = '%s/Current_current_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L)
       data_current_k_avg = Current_correlation_k_avg
       with open(filename_current_k_avg, 'wb') as outfile:
           pickle.dump(data_current_k_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_current_k_std_avg = '%s/Current_current_correlation_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L)
       data_current_k_std_avg = Current_correlation_k_std_avg
       with open(filename_current_k_std_avg, 'wb') as outfile:
           pickle.dump(data_current_k_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_current_k_text_avg = '%s/Current_current_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L)
       filename_current_k_variation_text_avg = '%s/Current_current_correlation_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_k,N,U,Mu,dtau,L)

       np.savetxt(filename_current_k_text_avg,Current_correlation_k_avg)
       np.savetxt(filename_current_k_variation_text_avg,Current_correlation_k_std_avg)
       
       print("run total",run_counter)
       
              
       #=======================================Saving data for analytic continuation ==========================================================
       
       F_bz,k_points = k_point_grid_full_bz(N)
       print("No of k points", k_points)
       print("Shape of Curr_k_avg",Current_correlation_k_avg.shape)
       
       for kk in range(k_points):
           kx = F_bz[kk,1]
           ky = F_bz[kk,2]
           Curr_k_data = np.copy(Current_correlation_k_avg[:,kk])
           Curr_k_std_data = np.copy(Current_correlation_k_std_avg[:,kk])
           Curr_data = np.stack((Tau,-1*Curr_k_data,Curr_k_std_data),axis = 1)
           filename_curr_k_point = '%s/Current_current_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_curr_curr_corr_k_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_curr_k_point,Curr_data)


       print("run total",run_counter)
   else:
       print("Error, momentum space current correlation function file not found")
       
 



def main(total,cmdargs):
    if(total!=4):
        raise ValueError('missing args')

    N = cmdargs[1]
    Dtau = cmdargs[2]
    run_no=int(cmdargs[3])

    u = ["5.30","5.40","5.50","5.60","5.70","6.00"] #["3.20","3.30","3.40","3.50","3.60","3.70","3.80","3.90","4.00","4.10","4.20","4.30","4.40","4.50","4.60","4.80","4.90","5.00","5.10","5.20","5.30","5.40","5.50","5.60","5.70","6.00"]
#    Trot = ["20","22","24","26","28","30","32","34","36","38","40","42","44","46","48","50","52","54","56","58","60"] 
    Trot = ["20","30","40"] #,"50","60"]

    mu = "0.00"
    #Text_dir_main = "../../Text_files"
    

    for k in range(len(Trot)):
            L = Trot[k]

            for j in range(len(u)):
            

                U = u[j]
                Text_dir_eqm = "../../../Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements"%(N,N,U,Dtau,mu,Dtau,L)

                filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U,mu,Dtau,L)
                with open(filename_eqm_avg, 'rb') as infile:
                     sys_measure_avg = pickle.load(infile)

                Nden = sys_measure_avg['Density averaged']
                print("avg sign, density, mu", Nden, mu)


                Text_dir_curr_curr_corr_k = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)

                Text_dir_curr_curr_corr_k_mf = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space_matsubara_frequency'%(N,N,U,Dtau,mu,Dtau,L)


    #===========================Directories for correlations for analytic continuation averages, imaginary time space =========================================
        
                Text_dir_curr_k_avg_ac = '../../../Text_files_fixed_density/Text_files_N_%s_analytic_continuation_current_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
                if not os.path.exists(Text_dir_curr_k_avg_ac):
                   os.makedirs(Text_dir_curr_k_avg_ac)
 
    
    # ====================================== Unequal time data averaging ==================================================================================
    #unequal_time_current_correlation_function_real_space_average(Text_dir_curr_curr_corr_r,Text_dir_curr_curr_corr_r_avg,N,U,mu,Dtau,L)
                unequal_time_current_correlation_function_momentum_space_average(Text_dir_curr_curr_corr_k,Text_dir_curr_k_avg_ac,N,U,mu,Dtau,L,run_no)

    # ================================================= Matsubara frequency space data averaging   ================================================================


if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)
