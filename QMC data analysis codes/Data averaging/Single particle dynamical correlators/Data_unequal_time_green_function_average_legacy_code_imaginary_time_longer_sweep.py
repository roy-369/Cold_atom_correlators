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



def unequal_time_greens_function_real_space_average(Text_dir_gf_r,Text_dir_gf_r_avg_ac,N,U,Mu,dtau,L,run_no):

       
   filename_gf_r_0 = '%s/Retarded_green_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_gf_r_0):
       
       with open(filename_gf_r_0, 'rb') as infile:
           green_function_r_0 = pickle.load(infile)

       filename_gf_r_variation_0 = '%s/Retarded_green_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       with open(filename_gf_r_variation_0, 'rb') as infile:
           green_function_r_std_0 = pickle.load(infile)


       Green_function_r = green_function_r_0.copy()
       Green_function_r_var = np.power(green_function_r_std_0,2)
  

       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
            run_counter=run_counter+1
            
            filename_gf_r = '%s/Retarded_green_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_gf_r):

               r_counter=r_counter+1
               with open(filename_gf_r, 'rb') as infile:
                   green_function_r = pickle.load(infile)

               filename_gf_r_variation = '%s/Retarded_green_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L,realization)
               with open(filename_gf_r_variation, 'rb') as infile:
                   green_function_r_std = pickle.load(infile)
 
 
               Green_function_X_r = np.add(Green_function_r,green_function_r)
               Green_function_r = Green_function_X_r.copy()

               Green_function_Y_r_var = np.power(green_function_r_std,2)
               Green_function_X_r_var = np.add(Green_function_r_var,Green_function_Y_r_var)
               Green_function_r_var = Green_function_X_r_var.copy()


       Green_function_r_avg = Green_function_r/r_counter
   
       Green_function_r_std_avg = np.sqrt(Green_function_r_var/(r_counter*r_counter))

       filename_eq_gf_r_avg = '%s/Retarded_green_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       data_eq_gf_r_avg = Green_function_r_avg
       with open(filename_eq_gf_r_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_r_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_eq_gf_r_std_avg = '%s/Retarded_green_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       data_eq_gf_r_std_avg = Green_function_r_std_avg
       with open(filename_eq_gf_r_std_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_r_std_avg, outfile, pickle.HIGHEST_PROTOCOL)



       filename_g_r_avg_text = '%s/Retarded_green_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       filename_g_r_variation_avg_text = '%s/Retarded_green_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_r,N,U,Mu,dtau,L)

       np.savetxt(filename_g_r_avg_text,Green_function_r_avg)
       np.savetxt(filename_g_r_variation_avg_text,Green_function_r_std_avg)


       print("run total",run_counter)

       #====================================Saving data for analytic continuation========================================

       timeslices = int(L)+1
       Tau = np.zeros(timeslices)
       for tt in range(timeslices):
           Tau[tt] = float(dtau)*tt

       UH_r,r_points = r_point_grid_upper_half(N)
       print("No of r points", r_points)
       print("Shape of G_r_avg", Green_function_r_avg.shape)


       for rr in range(r_points):
           rx = UH_r[rr,0]
           ry = UH_r[rr,1]

           G_r_data = np.copy(Green_function_r_avg[:,rr])
           G_r_std_data = np.copy(Green_function_r_std_avg[:,rr])

           G_r_data[-1] = -1*(1+G_r_data[0])
           G_r_std_data[-1] = G_r_std_data[0]

           Gf_data = np.stack((Tau,G_r_data,G_r_std_data),axis = 1)
           filename_g_r_point = '%s/Retarded_green_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_label_%s.dat' %(Text_dir_gf_r_avg_ac,N,U,Mu,dtau,L,str(rr))
           np.savetxt(filename_g_r_point,Gf_data)

       print("run total",run_counter)

   else:
       print("Error, real space green function file not found")


def unequal_time_greens_function_momentum_space_average(Text_dir_gf_k,Text_dir_gf_k_avg_ac,N,U,Mu,dtau,L,n_den,run_no):


   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_gf_k_0 = '%s/Retarded_green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)

   if os.path.exists(filename_gf_k_0):

       with open(filename_gf_k_0, 'rb') as infile:
           green_function_k_0 = pickle.load(infile)

       filename_gf_k_variation_0 = '%s/Retarded_green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       with open(filename_gf_k_variation_0, 'rb') as infile:
           green_function_k_std_0 = pickle.load(infile)


       Green_function_k = green_function_k_0.copy()
       Green_function_k_var = np.power(green_function_k_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1

            filename_gf_k = '%s/Retarded_green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_gf_k):

               r_counter=r_counter+1
               with open(filename_gf_k, 'rb') as infile:
                   green_function_k = pickle.load(infile)

               filename_gf_k_variation = '%s/Retarded_green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L,realization)
               with open(filename_gf_k_variation, 'rb') as infile:
                   green_function_k_std = pickle.load(infile)


               Green_function_X_k = np.add(Green_function_k,green_function_k)
               Green_function_k = Green_function_X_k.copy()

               Green_function_Y_k_var = np.power(green_function_k_std,2)
               Green_function_X_k_var = np.add(Green_function_k_var,Green_function_Y_k_var)
               Green_function_k_var = Green_function_X_k_var.copy()


       Green_function_k_avg = Green_function_k/r_counter

       Green_function_k_std_avg = np.sqrt(Green_function_k_var/(r_counter*r_counter))

       filename_eq_gf_k_avg = '%s/Retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       data_eq_gf_k_avg = Green_function_k_avg
       with open(filename_eq_gf_k_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_k_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_eq_gf_k_std_avg = '%s/Retarded_green_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       data_eq_gf_k_std_avg = Green_function_k_std_avg
       with open(filename_eq_gf_k_std_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_k_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_g_k_text_avg = '%s/Retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       filename_g_k_variation_text_avg = '%s/Retarded_green_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k,N,U,Mu,dtau,L)

       np.savetxt(filename_g_k_text_avg,Green_function_k_avg)
       np.savetxt(filename_g_k_variation_text_avg,Green_function_k_std_avg)

       #=======================================Saving fata for analytic continuation ==========================================================
       UH_bz,k_points = k_point_grid_upper_half_bz(N)
       print("No of k points", k_points)
       print("Shape of G_k_avg", Green_function_k_avg.shape)
       
       M_1 = np.zeros(k_points)
       M_2 = np.zeros(k_points)
       Mean_ac = np.zeros(k_points)
       Variance_ac = np.zeros(k_points)
       
       for kk in range(k_points):
           kx = (2*np.pi/int(N))*UH_bz[kk,1]
           ky = (2*np.pi/int(N))*UH_bz[kk,2]
           
           G_k_data = np.copy(Green_function_k_avg[:,kk])
           G_k_std_data = np.copy(Green_function_k_std_avg[:,kk])
           
           G_k_data[-1] = -1*(1+G_k_data[0])
           G_k_std_data[-1] = G_k_std_data[0]

           Gf_data = np.stack((Tau,G_k_data,G_k_std_data),axis = 1)
           filename_g_k_point = '%s/Retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_gf_k_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_g_k_point,Gf_data)
           
           ep_k = -2*(np.cos(kx)+np.cos(ky))
           m1 = ep_k-float(Mu)-0.5*float(U)+0.5*n_den*float(U)
           m2 = (ep_k-float(Mu)-0.5*float(U))**2+float(U)*(ep_k-float(Mu)-0.5*float(U))*n_den+0.5*float(U)*float(U)*n_den
           M_1[kk] = m1
           M_2[kk] = m2
           Mean_ac[kk] = m1
           Variance_ac[kk] = 2*np.sqrt(m2-m1*m1)
       np.savetxt("%s/Default_model_mean_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_avg_ac,N,U,Mu,dtau,L),Mean_ac)
       np.savetxt("%s/Default_model_variance_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_avg_ac,N,U,Mu,dtau,L),Variance_ac)
       print("run total",run_counter)
   else:
       print("Error, momentum space green function file not found")





def main(total,cmdargs):
    if(total!=7):
        raise ValueError('missing args')

    N = cmdargs[1]
    U = cmdargs[2]
    mu = cmdargs[3]
    L = cmdargs[4]
    Dtau = cmdargs[5]
    run_no = int(cmdargs[6])
    print(U,"U")
    print(L,"L")
    Beta = str((float(L))*(float(Dtau)))
    print(Beta,"Beta")


    Text_dir_eqm = "../../../Text_files_longer_sweep/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements"%(N,N,U,Dtau,mu,Dtau,L)

    filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U,mu,Dtau,L)
    with open(filename_eqm_avg, 'rb') as infile:
         sys_measure_avg = pickle.load(infile)

    Nden = sys_measure_avg['Density averaged']
    print("avg sign, density, mu", Nden, mu)


    #==========================Diretories for real space correlation averages =============================================

    Text_dir_gf_r = '../../../Text_files_longer_sweep/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)


    #==================== Directories for momentum space correlation averages ===================================================

    Text_dir_gf_k = '../../../Text_files_longer_sweep/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)

    #===========================Directories for correlations for analytic continuation averages, imaginary time space =========================================

    Text_dir_gf_r_avg_ac = '../../../Text_files_longer_sweep/Text_files_N_%s_analytic_continuation/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_real_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_r_avg_ac):
        os.makedirs(Text_dir_gf_r_avg_ac)

    Text_dir_gf_k_avg_ac = '../../../Text_files_longer_sweep/Text_files_N_%s_analytic_continuation/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_k_avg_ac):
        os.makedirs(Text_dir_gf_k_avg_ac)
  
    
    unequal_time_greens_function_real_space_average(Text_dir_gf_r,Text_dir_gf_r_avg_ac,N,U,mu,Dtau,L,run_no)
    unequal_time_greens_function_momentum_space_average(Text_dir_gf_k,Text_dir_gf_k_avg_ac,N,U,mu,Dtau,L,Nden,run_no)



if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)
