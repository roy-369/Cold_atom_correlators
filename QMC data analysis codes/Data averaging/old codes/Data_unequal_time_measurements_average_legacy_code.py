#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:11:44 2022

@author: roy.369
"""


import numpy as np
import pickle

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

       
   filename_gf_r_0 = '%s/Retarded_Green_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_gf_r_0):
       
       with open(filename_gf_r_0, 'rb') as infile:
           green_function_r_0 = pickle.load(infile)

       filename_gf_r_variation_0 = '%s/Retarded_Green_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       with open(filename_gf_r_variation_0, 'rb') as infile:
           green_function_r_std_0 = pickle.load(infile)


       Green_function_r = green_function_r_0.copy()
       Green_function_r_var = np.power(green_function_r_std_0,2)
  

       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
            run_counter=run_counter+1
            
            filename_gf_r = '%s/Retarded_Green_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_gf_r):

               r_counter=r_counter+1
               with open(filename_gf_r, 'rb') as infile:
                   green_function_r = pickle.load(infile)

               filename_gf_r_variation = '%s/Retarded_Green_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L,realization)
               with open(filename_gf_r_variation, 'rb') as infile:
                   green_function_r_std = pickle.load(infile)
 
 
               Green_function_X_r = np.add(Green_function_r,green_function_r)
               Green_function_r = Green_function_X_r.copy()

               Green_function_Y_r_var = np.power(green_function_r_std,2)
               Green_function_X_r_var = np.add(Green_function_r_var,Green_function_Y_r_var)
               Green_function_r_var = Green_function_X_r_var.copy()


       Green_function_r_avg = Green_function_r/r_counter
   
       Green_function_r_std_avg = np.sqrt(Green_function_r_var/(r_counter*r_counter))

       filename_eq_gf_r_avg = '%s/Retarded_Green_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       data_eq_gf_r_avg = Green_function_r_avg
       with open(filename_eq_gf_r_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_r_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_eq_gf_r_std_avg = '%s/Retarded_Green_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       data_eq_gf_r_std_avg = Green_function_r_std_avg
       with open(filename_eq_gf_r_std_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_r_std_avg, outfile, pickle.HIGHEST_PROTOCOL)



       filename_g_r_avg_text = '%s/Retarded_Green_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_r,N,U,Mu,dtau,L)
       filename_g_r_variation_avg_text = '%s/Retarded_Green_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_r,N,U,Mu,dtau,L)

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
           for t in range(len(Tau)):
               if G_r_std_data[t] < 0.000001:
                  G_r_std_data[t] = 0.000001 

           Gf_data = np.stack((Tau,G_r_data,G_r_std_data),axis = 1)
           filename_g_r_point = '%s/Retarded_Green_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_label_%s.dat' %(Text_dir_gf_r_avg_ac,N,U,Mu,dtau,L,str(rr))
           np.savetxt(filename_g_r_point,Gf_data)

       print("run total",run_counter)

   else:
       print("Error, real space green function file not found")


def unequal_time_greens_function_momentum_space_average(Text_dir_gf_k,Text_dir_gf_k_avg_ac,N,U,Mu,dtau,L,n_den,run_no):


   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_gf_k_0 = '%s/Retarded_Green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)

   if os.path.exists(filename_gf_k_0):

       with open(filename_gf_k_0, 'rb') as infile:
           green_function_k_0 = pickle.load(infile)

       filename_gf_k_variation_0 = '%s/Retarded_Green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       with open(filename_gf_k_variation_0, 'rb') as infile:
           green_function_k_std_0 = pickle.load(infile)


       Green_function_k = green_function_k_0.copy()
       Green_function_k_var = np.power(green_function_k_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1

            filename_gf_k = '%s/Retarded_Green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_gf_k):

               r_counter=r_counter+1
               with open(filename_gf_k, 'rb') as infile:
                   green_function_k = pickle.load(infile)

               filename_gf_k_variation = '%s/Retarded_Green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L,realization)
               with open(filename_gf_k_variation, 'rb') as infile:
                   green_function_k_std = pickle.load(infile)


               Green_function_X_k = np.add(Green_function_k,green_function_k)
               Green_function_k = Green_function_X_k.copy()

               Green_function_Y_k_var = np.power(green_function_k_std,2)
               Green_function_X_k_var = np.add(Green_function_k_var,Green_function_Y_k_var)
               Green_function_k_var = Green_function_X_k_var.copy()


       Green_function_k_avg = Green_function_k/r_counter

       Green_function_k_std_avg = np.sqrt(Green_function_k_var/(r_counter*r_counter))

       filename_eq_gf_k_avg = '%s/Retarded_Green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       data_eq_gf_k_avg = Green_function_k_avg
       with open(filename_eq_gf_k_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_k_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_eq_gf_k_std_avg = '%s/Retarded_Green_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       data_eq_gf_k_std_avg = Green_function_k_std_avg
       with open(filename_eq_gf_k_std_avg, 'wb') as outfile:
           pickle.dump(data_eq_gf_k_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_g_k_text_avg = '%s/Retarded_Green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       filename_g_k_variation_text_avg = '%s/Retarded_Green_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k,N,U,Mu,dtau,L)

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
           Gf_data = np.stack((Tau,G_k_data,G_k_std_data),axis = 1)
           filename_g_k_point = '%s/Retarded_Green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_gf_k_avg_ac,N,U,Mu,dtau,L,str(kk))
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


def unequal_time_density_correlation_function_real_space_average(Text_dir_den_corr_r,N,U,Mu,dtau,L,run_no):

   filename_den_r_0 = '%s/Density_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_den_r_0):

       with open(filename_den_r_0, 'rb') as infile:
           density_correlation_r_0 = pickle.load(infile)

       filename_den_r_variation_0 = '%s/Density_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L)
       with open(filename_den_r_variation_0, 'rb') as infile:
           density_correlation_r_std_0 = pickle.load(infile)


       Density_correlation_r = density_correlation_r_0.copy()
       Density_correlation_r_var = np.power(density_correlation_r_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
            run_counter=run_counter+1
            
            filename_den_r = '%s/Density_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_den_r):

               r_counter=r_counter+1
               with open(filename_den_r, 'rb') as infile:
                   density_correlation_r = pickle.load(infile)

               filename_den_r_variation = '%s/Density_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L,realization)
               with open(filename_den_r_variation, 'rb') as infile:
                   density_correlation_r_std = pickle.load(infile)


               Density_correlation_X_r = np.add(Density_correlation_r,density_correlation_r)
               Density_correlation_r = Density_correlation_X_r.copy()

               Density_correlation_Y_r_var = np.power(density_correlation_r_std,2)
               Density_correlation_X_r_var = np.add(Density_correlation_r_var,Density_correlation_Y_r_var)
               Density_correlation_r_var = Density_correlation_X_r_var.copy()


       Density_correlation_r_avg = Density_correlation_r/r_counter

       Density_correlation_r_std_avg = np.sqrt(Density_correlation_r_var/(r_counter*r_counter))

       filename_den_r_avg = '%s/Density_correlation_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L)
       data_den_r_avg = Density_correlation_r_avg
       with open(filename_den_r_avg, 'wb') as outfile:
           pickle.dump(data_den_r_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_den_r_std_avg = '%s/Density_correlation_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_den_corr_r,N,U,Mu,dtau,L)
       data_den_r_std_avg = Density_correlation_r_std_avg
       with open(filename_den_r_std_avg, 'wb') as outfile:
           pickle.dump(data_den_r_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       print("run total",run_counter)
   else:
       print("Error, real space density correlation function file not found")



def unequal_time_density_correlation_function_momentum_space_average(Text_dir_den_corr_k,Text_dir_den_corr_k_avg_ac,N,U,Mu,dtau,L,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_den_k_0 = '%s/Density_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L)

   if os.path.exists(filename_den_k_0):

       with open(filename_den_k_0, 'rb') as infile:
           density_correlation_k_0 = pickle.load(infile)

       filename_den_k_variation_0 = '%s/Density_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L)
       with open(filename_den_k_variation_0, 'rb') as infile:
           density_correlation_k_std_0 = pickle.load(infile)


       Density_correlation_k = density_correlation_k_0.copy()
       Density_correlation_k_var = np.power(density_correlation_k_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            filename_den_k = '%s/Density_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_den_k):

               r_counter=r_counter+1
               with open(filename_den_k, 'rb') as infile:
                   density_correlation_k = pickle.load(infile)

               filename_den_k_variation = '%s/Density_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L,realization)
               with open(filename_den_k_variation, 'rb') as infile:
                   density_correlation_k_std = pickle.load(infile)


               Density_correlation_X_k = np.add(Density_correlation_k,density_correlation_k)
               Density_correlation_k = Density_correlation_X_k.copy()

               Density_correlation_Y_k_var = np.power(density_correlation_k_std,2)
               Density_correlation_X_k_var = np.add(Density_correlation_k_var,Density_correlation_Y_k_var)
               Density_correlation_k_var = Density_correlation_X_k_var.copy()


       Density_correlation_k_avg = Density_correlation_k/r_counter

       Density_correlation_k_std_avg = np.sqrt(Density_correlation_k_var/(r_counter*r_counter))

       filename_den_k_avg = '%s/Density_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L)
       data_den_k_avg = Density_correlation_k_avg
       with open(filename_den_k_avg, 'wb') as outfile:
           pickle.dump(data_den_k_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_den_k_std_avg = '%s/Density_correlation_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_den_corr_k,N,U,Mu,dtau,L)
       data_den_k_std_avg = Density_correlation_k_std_avg
       with open(filename_den_k_std_avg, 'wb') as outfile:
           pickle.dump(data_den_k_std_avg, outfile, pickle.HIGHEST_PROTOCOL)
       
       
       #=======================================Saving fata for analytic continuation ==========================================================
       
       UH_bz,k_points = k_point_grid_upper_half_bz(N)
       print("No of k points", k_points)
       print("Shape of Den_k_avg",Density_correlation_k_avg.shape)
       
       for kk in range(k_points):
           kx = UH_bz[kk,1]
           ky = UH_bz[kk,2]
           Den_k_data = np.copy(Density_correlation_k_avg[:,kk])
           Den_k_std_data = np.copy(Density_correlation_k_std_avg[:,kk])
           Den_data = np.stack((Tau,Den_k_data,Den_k_std_data),axis = 1)
           filename_den_k_point = '%s/Density_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_den_corr_k_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_den_k_point,Den_data)


       print("run total",run_counter)
   else:
       print("Error, momentum space density correlation function file not found")



def unequal_time_chi_xx_correlation_function_real_space_average(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L,run_no):

   filename_chi_xx_r_0 = '%s/Spin_xx_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_chi_xx_r_0):

       with open(filename_chi_xx_r_0, 'rb') as infile:
           chi_xx_correlation_r_0 = pickle.load(infile)

       filename_chi_xx_r_variation_0 = '%s/Spin_xx_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L)
       with open(filename_chi_xx_r_variation_0, 'rb') as infile:
           chi_xx_correlation_r_std_0 = pickle.load(infile)


       Chi_xx_correlation_r = chi_xx_correlation_r_0.copy()
       Chi_xx_correlation_r_var = np.power(chi_xx_correlation_r_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1

            filename_chi_xx_r = '%s/Spin_xx_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_chi_xx_r):

               r_counter=r_counter+1
               with open(filename_chi_xx_r, 'rb') as infile:
                   chi_xx_correlation_r = pickle.load(infile)

               filename_chi_xx_r_variation = '%s/Spin_xx_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L,realization)
               with open(filename_chi_xx_r_variation, 'rb') as infile:
                   chi_xx_correlation_r_std = pickle.load(infile)


               Chi_xx_correlation_X_r = np.add(Chi_xx_correlation_r,chi_xx_correlation_r)
               Chi_xx_correlation_r = Chi_xx_correlation_X_r.copy()

               Chi_xx_correlation_Y_r_var = np.power(chi_xx_correlation_r_std,2)
               Chi_xx_correlation_X_r_var = np.add(Chi_xx_correlation_r_var,Chi_xx_correlation_Y_r_var)
               Chi_xx_correlation_r_var = Chi_xx_correlation_X_r_var.copy()


       Chi_xx_correlation_r_avg = Chi_xx_correlation_r/r_counter

       Chi_xx_correlation_r_std_avg = np.sqrt(Chi_xx_correlation_r_var/(r_counter*r_counter))

       filename_chi_xx_r_avg = '%s/Spin_xx_correlation_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L)
       data_chi_xx_r_avg = Chi_xx_correlation_r_avg
       with open(filename_chi_xx_r_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_r_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_chi_xx_r_std_avg = '%s/Spin_xx_correlation_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_r,N,U,Mu,dtau,L)
       data_chi_xx_r_std_avg = Chi_xx_correlation_r_std_avg
       with open(filename_chi_xx_r_std_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_r_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       print("run total",run_counter)
   else:
       print("Error, real space Chi XX function file not found")



def unequal_time_chi_xx_correlation_function_momentum_space_average(Text_dir_chi_xx_corr_k,Text_dir_chi_xx_corr_k_avg_ac,N,U,Mu,dtau,L,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_chi_xx_k_0 = '%s/Spin_xx_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L)

   if os.path.exists(filename_chi_xx_k_0):

       with open(filename_chi_xx_k_0, 'rb') as infile:
           chi_xx_correlation_k_0 = pickle.load(infile)

       filename_chi_xx_k_variation_0 = '%s/Spin_xx_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L)
       with open(filename_chi_xx_k_variation_0, 'rb') as infile:
           chi_xx_correlation_k_std_0 = pickle.load(infile)


       Chi_xx_correlation_k = chi_xx_correlation_k_0.copy()
       Chi_xx_correlation_k_var = np.power(chi_xx_correlation_k_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            filename_chi_xx_k = '%s/Spin_xx_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_chi_xx_k):

               r_counter=r_counter+1
               with open(filename_chi_xx_k, 'rb') as infile:
                   chi_xx_correlation_k = pickle.load(infile)

               filename_chi_xx_k_variation = '%s/Spin_xx_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L,realization)
               with open(filename_chi_xx_k_variation, 'rb') as infile:
                   chi_xx_correlation_k_std = pickle.load(infile)


               Chi_xx_correlation_X_k = np.add(Chi_xx_correlation_k,chi_xx_correlation_k)
               Chi_xx_correlation_k = Chi_xx_correlation_X_k.copy()

               Chi_xx_correlation_Y_k_var = np.power(chi_xx_correlation_k_std,2)
               Chi_xx_correlation_X_k_var = np.add(Chi_xx_correlation_k_var,Chi_xx_correlation_Y_k_var)
               Chi_xx_correlation_k_var = Chi_xx_correlation_X_k_var.copy()


       Chi_xx_correlation_k_avg = Chi_xx_correlation_k/r_counter

       Chi_xx_correlation_k_std_avg = np.sqrt(Chi_xx_correlation_k_var/(r_counter*r_counter))

       filename_chi_xx_k_avg = '%s/Spin_xx_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L)
       data_chi_xx_k_avg = Chi_xx_correlation_k_avg
       with open(filename_chi_xx_k_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_k_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_chi_xx_k_std_avg = '%s/Spin_xx_correlation_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_k,N,U,Mu,dtau,L)
       data_chi_xx_k_std_avg = Chi_xx_correlation_k_std_avg
       with open(filename_chi_xx_k_std_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_k_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       print("run total",run_counter)

       
       #=======================================Saving fata for analytic continuation ==========================================================
       
       UH_bz,k_points = k_point_grid_upper_half_bz(N)
       print("No of k points", k_points)
       print("Shape of Chi_xx_k_avg",Chi_xx_correlation_k_avg.shape)
       
       for kk in range(k_points):
           kx = UH_bz[kk,1]
           ky = UH_bz[kk,2]
           Chi_xx_k_data = np.copy(Chi_xx_correlation_k_avg[:,kk])
           Chi_xx_k_std_data = np.copy(Chi_xx_correlation_k_std_avg[:,kk])
           Chi_xx_data = np.stack((Tau,Chi_xx_k_data,Chi_xx_k_std_data),axis = 1)
           filename_chi_xx_k_point = '%s/Spin_xx_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_chi_xx_corr_k_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_chi_xx_k_point,Chi_xx_data)


       print("run total",run_counter)
   else:
       print("Error, momentum space chi xx correlation function file not found")
       


def unequal_time_chi_zz_correlation_function_real_space_average(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L,run_no):

   filename_chi_zz_r_0 = '%s/Spin_zz_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_chi_zz_r_0):

       with open(filename_chi_zz_r_0, 'rb') as infile:
           chi_zz_correlation_r_0 = pickle.load(infile)

       filename_chi_zz_r_variation_0 = '%s/Spin_zz_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L)
       with open(filename_chi_zz_r_variation_0, 'rb') as infile:
           chi_zz_correlation_r_std_0 = pickle.load(infile)


       Chi_zz_correlation_r = chi_zz_correlation_r_0.copy()
       Chi_zz_correlation_r_var = np.power(chi_zz_correlation_r_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            filename_chi_zz_r = '%s/Spin_zz_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_chi_zz_r):

               r_counter=r_counter+1
               with open(filename_chi_zz_r, 'rb') as infile:
                   chi_zz_correlation_r = pickle.load(infile)

               filename_chi_zz_r_variation = '%s/Spin_zz_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L,realization)
               with open(filename_chi_zz_r_variation, 'rb') as infile:
                   chi_zz_correlation_r_std = pickle.load(infile)


               Chi_zz_correlation_X_r = np.add(Chi_zz_correlation_r,chi_zz_correlation_r)
               Chi_zz_correlation_r = Chi_zz_correlation_X_r.copy()

               Chi_zz_correlation_Y_r_var = np.power(chi_zz_correlation_r_std,2)
               Chi_zz_correlation_X_r_var = np.add(Chi_zz_correlation_r_var,Chi_zz_correlation_Y_r_var)
               Chi_zz_correlation_r_var = Chi_zz_correlation_X_r_var.copy()


       Chi_zz_correlation_r_avg = Chi_zz_correlation_r/r_counter

       Chi_zz_correlation_r_std_avg = np.sqrt(Chi_zz_correlation_r_var/(r_counter*r_counter))

       filename_chi_zz_r_avg = '%s/Spin_zz_correlation_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L)
       data_chi_zz_r_avg = Chi_zz_correlation_r_avg
       with open(filename_chi_zz_r_avg, 'wb') as outfile:
           pickle.dump(data_chi_zz_r_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_chi_zz_r_std_avg = '%s/Spin_zz_correlation_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L)
       data_chi_zz_r_std_avg = Chi_zz_correlation_r_std_avg
       with open(filename_chi_zz_r_std_avg, 'wb') as outfile:
           pickle.dump(data_chi_zz_r_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_chi_zz_r_text_avg = '%s/Spin_zz_correlation_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L)
       filename_chi_zz_r_variation_text_avg = '%s/Spin_zz_correlation_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_zz_corr_r,N,U,Mu,dtau,L)

       np.savetxt(filename_chi_zz_r_text_avg,Chi_zz_correlation_r_avg)
       np.savetxt(filename_chi_zz_r_variation_text_avg,Chi_zz_correlation_r_std_avg)

       
       print("run total",run_counter)
   else:
       print("Error, real space Chi ZZ function file not found")



def unequal_time_chi_zz_correlation_function_momentum_space_average(Text_dir_chi_zz_corr_k,Text_dir_chi_zz_corr_k_avg_ac,N,U,Mu,dtau,L,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_chi_zz_k_0 = '%s/Spin_zz_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L)

   if os.path.exists(filename_chi_zz_k_0):

       with open(filename_chi_zz_k_0, 'rb') as infile:
           chi_zz_correlation_k_0 = pickle.load(infile)

       filename_chi_zz_k_variation_0 = '%s/Spin_zz_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L)
       with open(filename_chi_zz_k_variation_0, 'rb') as infile:
           chi_zz_correlation_k_std_0 = pickle.load(infile)


       Chi_zz_correlation_k = chi_zz_correlation_k_0.copy()
       Chi_zz_correlation_k_var = np.power(chi_zz_correlation_k_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            filename_chi_zz_k = '%s/Spin_zz_correlation_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_chi_zz_k):

               r_counter=r_counter+1
               with open(filename_chi_zz_k, 'rb') as infile:
                   chi_zz_correlation_k = pickle.load(infile)

               filename_chi_zz_k_variation = '%s/Spin_zz_correlation_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L,realization)
               with open(filename_chi_zz_k_variation, 'rb') as infile:
                   chi_zz_correlation_k_std = pickle.load(infile)


               Chi_zz_correlation_X_k = np.add(Chi_zz_correlation_k,chi_zz_correlation_k)
               Chi_zz_correlation_k = Chi_zz_correlation_X_k.copy()

               Chi_zz_correlation_Y_k_var = np.power(chi_zz_correlation_k_std,2)
               Chi_zz_correlation_X_k_var = np.add(Chi_zz_correlation_k_var,Chi_zz_correlation_Y_k_var)
               Chi_zz_correlation_k_var = Chi_zz_correlation_X_k_var.copy()


       Chi_zz_correlation_k_avg = Chi_zz_correlation_k/r_counter

       Chi_zz_correlation_k_std_avg = np.sqrt(Chi_zz_correlation_k_var/(r_counter*r_counter))

       filename_chi_zz_k_avg = '%s/Spin_zz_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L)
       data_chi_zz_k_avg = Chi_zz_correlation_k_avg
       with open(filename_chi_zz_k_avg, 'wb') as outfile:
           pickle.dump(data_chi_zz_k_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_chi_zz_k_std_avg = '%s/Spin_zz_correlation_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L)
       data_chi_zz_k_std_avg = Chi_zz_correlation_k_std_avg
       with open(filename_chi_zz_k_std_avg, 'wb') as outfile:
           pickle.dump(data_chi_zz_k_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_chi_zz_k_text_avg = '%s/Spin_zz_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L)
       filename_chi_zz_k_variation_text_avg = '%s/Spin_zz_correlation_function_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_zz_corr_k,N,U,Mu,dtau,L)

       np.savetxt(filename_chi_zz_k_text_avg,Chi_zz_correlation_k_avg)
       np.savetxt(filename_chi_zz_k_variation_text_avg,Chi_zz_correlation_k_std_avg)
       
       print("run total",run_counter)

       
       #=======================================Saving fata for analytic continuation ==========================================================
       
       UH_bz,k_points = k_point_grid_upper_half_bz(N)
       print("No of k points", k_points)
       print("Shape of Chi_zz_k_avg",Chi_zz_correlation_k_avg.shape)
       
       for kk in range(k_points):
           kx = UH_bz[kk,1]
           ky = UH_bz[kk,2]
           Chi_zz_k_data = np.copy(Chi_zz_correlation_k_avg[:,kk])
           Chi_zz_k_std_data = np.copy(Chi_zz_correlation_k_std_avg[:,kk])
           Chi_zz_data = np.stack((Tau,Chi_zz_k_data,Chi_zz_k_std_data),axis = 1)
           filename_chi_zz_k_point = '%s/Spin_zz_correlation_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_chi_zz_corr_k_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_chi_zz_k_point,Chi_zz_data)


       print("run total",run_counter)
   else:
       print("Error, momentum space chi zz correlation function file not found")




def unequal_time_current_correlation_function_real_space_average(Text_dir_curr_curr_corr_r,N,U,Mu,dtau,L,run_no):

   filename_curr_r_0 = '%s/Current_current_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_curr_r_0):

       with open(filename_curr_r_0, 'rb') as infile:
           current_correlation_r_0 = pickle.load(infile)

       filename_current_r_variation_0 = '%s/Current_current_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_r,N,U,Mu,dtau,L)
       with open(filename_current_r_variation_0, 'rb') as infile:
           current_correlation_r_std_0 = pickle.load(infile)


       Current_correlation_r = current_correlation_r_0.copy()
       Current_correlation_r_var = np.power(current_correlation_r_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            filename_curr_r = '%s/Current_current_correlation_function_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_curr_r):

               r_counter=r_counter+1
               with open(filename_curr_r, 'rb') as infile:
                   current_correlation_r = pickle.load(infile)

               filename_current_r_variation = '%s/Current_current_correlation_function_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_r,N,U,Mu,dtau,L,realization)
               with open(filename_current_r_variation, 'rb') as infile:
                   current_correlation_r_std = pickle.load(infile)


               Current_correlation_X_r = np.add(Current_correlation_r,current_correlation_r)
               Current_correlation_r = Current_correlation_X_r.copy()

               Current_correlation_Y_r_var = np.power(current_correlation_r_std,2)
               Current_correlation_X_r_var = np.add(Current_correlation_r_var,Current_correlation_Y_r_var)
               Current_correlation_r_var = Current_correlation_X_r_var.copy()


       Current_correlation_r_avg = Current_correlation_r/r_counter

       Current_correlation_r_std_avg = np.sqrt(Current_correlation_r_var/(r_counter*r_counter))

       filename_current_r_avg = '%s/Current_current_correlation_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_r_avg,N,U,Mu,dtau,L)
       data_current_r_avg = Current_correlation_r_avg
       with open(filename_current_r_avg, 'wb') as outfile:
           pickle.dump(data_current_r_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_current_r_std_avg = '%s/Current_current_correlation_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_r,N,U,Mu,dtau,L)
       data_current_r_std_avg = Current_correlation_r_std_avg
       with open(filename_current_r_std_avg, 'wb') as outfile:
           pickle.dump(data_current_r_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_current_r_text_avg = '%s/Current_current_correlation_function_real_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_r,N,U,Mu,dtau,L)
       filename_current_r_variation_text_avg = '%s/Current_current_correlation_function_real_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_r_avg,N,U,Mu,dtau,L)

       np.savetxt(filename_current_r_text_avg,Current_correlation_r_avg)
       np.savetxt(filename_current_r_variation_text_avg,Current_correlation_r_std_avg)
       
       print("run total",run_counter)
   else:
       print("Error, real space current current correlation function file not found")


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
       
 
#============================================ Averaging data in Matsubara frequency space ===============================

def matsubara_frequency_green_function_momentum_space_average(Text_dir_gf_k_mf,Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,n_den,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_gf_k_real_0 = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)

   if os.path.exists(filename_gf_k_real_0):

       with open(filename_gf_k_real_0, 'rb') as infile:
           green_function_k_real_0 = pickle.load(infile)

       filename_gf_k_real_variation_0 = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       with open(filename_gf_k_real_variation_0, 'rb') as infile:
           green_function_k_real_std_0 = pickle.load(infile)

       filename_gf_k_imag_0 = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       with open(filename_gf_k_imag_0, 'rb') as infile:
           green_function_k_imag_0 = pickle.load(infile)

       filename_gf_k_imag_variation_0 = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       with open(filename_gf_k_imag_variation_0, 'rb') as infile:
           green_function_k_imag_std_0 = pickle.load(infile)


       Green_function_k_real = green_function_k_real_0.copy()
       Green_function_k_real_var = np.power(green_function_k_real_std_0,2)

       Green_function_k_imag = green_function_k_imag_0.copy()
       Green_function_k_imag_var = np.power(green_function_k_imag_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            
            filename_gf_k_real = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_gf_k_real):

               r_counter=r_counter+1
               with open(filename_gf_k_real, 'rb') as infile:
                   green_function_k_real = pickle.load(infile)

               filename_gf_k_real_variation = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_gf_k_real_variation, 'rb') as infile:
                   green_function_k_real_std = pickle.load(infile)

               filename_gf_k_imag = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_gf_k_imag, 'rb') as infile:
                   green_function_k_imag = pickle.load(infile)

               filename_gf_k_imag_variation = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_gf_k_imag_variation, 'rb') as infile:
                   green_function_k_imag_std = pickle.load(infile)
                   
                   
               Green_function_X_k_real = np.add(Green_function_k_real,green_function_k_real)
               Green_function_k_real = Green_function_X_k_real.copy()

               Green_function_Y_k_real_var = np.power(green_function_k_real_std,2)
               Green_function_X_k_real_var = np.add(Green_function_k_real_var,Green_function_Y_k_real_var)
               Green_function_k_real_var = Green_function_X_k_real_var.copy()

               Green_function_X_k_imag = np.add(Green_function_k_imag,green_function_k_imag)
               Green_function_k_imag = Green_function_X_k_imag.copy()

               Green_function_Y_k_imag_var = np.power(green_function_k_imag_std,2)
               Green_function_X_k_imag_var = np.add(Green_function_k_imag_var,Green_function_Y_k_imag_var)
               Green_function_k_imag_var = Green_function_X_k_imag_var.copy()
               

       Green_function_k_real_avg = Green_function_k_real/r_counter

       Green_function_k_real_std_avg = np.sqrt(Green_function_k_real_var/(r_counter*r_counter))

       Green_function_k_imag_avg = Green_function_k_imag/r_counter

       Green_function_k_imag_std_avg = np.sqrt(Green_function_k_imag_var/(r_counter*r_counter))
       
       #=======================================Saving real part=======================================================================
       
       filename_gf_k_real_avg = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       data_gf_k_real_avg = Green_function_k_real_avg
       with open(filename_gf_k_real_avg, 'wb') as outfile:
           pickle.dump(data_gf_k_real_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_gf_k_real_std_avg = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       data_gf_k_real_std_avg = Green_function_k_real_std_avg
       with open(filename_gf_k_real_std_avg, 'wb') as outfile:
           pickle.dump(data_gf_k_real_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_gf_k_real_text_avg = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       filename_gf_k_real_variation_text_avg = '%s/Retarded_green_function_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_gf_k_real_text_avg,Green_function_k_real_avg)
       np.savetxt(filename_gf_k_real_variation_text_avg,Green_function_k_real_std_avg)
       
       #=====================================Saving imaginary part =====================================================================
       
       filename_gf_k_imag_avg = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       data_gf_k_imag_avg = Green_function_k_imag_avg
       with open(filename_gf_k_imag_avg, 'wb') as outfile:
           pickle.dump(data_gf_k_imag_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_gf_k_imag_std_avg = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       data_gf_k_imag_std_avg = Green_function_k_imag_std_avg
       with open(filename_gf_k_imag_std_avg, 'wb') as outfile:
           pickle.dump(data_gf_k_imag_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_gf_k_imag_text_avg = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)
       filename_gf_k_imag_variation_text_avg = '%s/Retarded_green_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_gf_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_gf_k_imag_text_avg,Green_function_k_imag_avg)
       np.savetxt(filename_gf_k_imag_variation_text_avg,Green_function_k_imag_std_avg)
       
       
       #======================================Saving data for analytic continuation ====================================================
       print("run total",run_counter)
       beta = float(dtau)*float(L)
       
       F_bz,k_points = k_point_grid_upper_half_bz(N)
       print("No of k points", k_points)
       print("Shape of Gf_k_avg",Green_function_k_real_avg.shape)
       Omega_n = np.zeros(int(int(L)/2+1))
       
       for n in range(len(Omega_n)):
           Omega_n[n] = 2*np.pi*(n+1)/beta
           
       Mean_ac = np.zeros(k_points)
       Variance_ac = np.zeros(k_points)
       Norm_self_energy = np.zeros(k_points)
       
       for kk in range(k_points):
           kx = F_bz[kk,1]
           ky = F_bz[kk,2]
           
           Gf_k_real_data = Green_function_k_real_avg[:,kk]
           Gf_k_real_std_data = Green_function_k_real_std_avg[:,kk]
           Gf_k_imag_data = Green_function_k_imag_avg[:,kk]
           Gf_k_imag_std_data = Green_function_k_imag_std_avg[:,kk]
           

           Gf_data_real = np.stack((Omega_n,Gf_k_real_data,Gf_k_real_std_data),axis = 1)
           filename_gf_k_point_real = '%s/Retarded_green_function_real_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_gf_k_point_real,Gf_data_real)
           
           Gf_data_imag = np.stack((Omega_n,Gf_k_imag_data,Gf_k_imag_std_data),axis = 1)
           filename_gf_k_point_imag = '%s/Retarded_green_function_imaginary_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_gf_k_point_imag,Gf_data_imag)
           
           Gf_data_ac = np.stack((Omega_n,Gf_k_real_data,Gf_k_real_std_data,Gf_k_imag_data,Gf_k_imag_std_data),axis = 1)
           filename_gf_k_point = '%s/Retarded_green_function_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_gf_k_point,Gf_data_ac)
           
           ep_k = -2*(np.cos(kx)+np.cos(ky))
           m1 = ep_k-float(Mu)-0.5*float(U)+0.5*n_den*float(U)
           m2 = (ep_k-float(Mu)-0.5*float(U))**2+float(U)*(ep_k-float(Mu)-0.5*float(U))*n_den+0.5*float(U)*float(U)*n_den

           Mean_ac[kk] = m1
           Variance_ac[kk] = 2*np.sqrt(m2-m1*m1)
           Norm_self_energy[kk] = float(U)*float(U)*n_den*(1-n_den)
           
       np.savetxt("%s/Default_model_mean_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L),Mean_ac)
       np.savetxt("%s/Default_model_variance_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L),Variance_ac)
       np.savetxt("%s/Self_energy_norm_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat"%(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L),Norm_self_energy)
       print("run total",run_counter)
   else:
       print("Error, momentum space green function file not found")



def matsubara_frequency_self_energy_momentum_space_average(Text_dir_sigma_k_mf,Text_dir_sigma_k_mf_avg_ac,N,U,Mu,dtau,L,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_sigma_k_real_0 = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)

   if os.path.exists(filename_sigma_k_real_0):

       with open(filename_sigma_k_real_0, 'rb') as infile:
           self_energy_k_real_0 = pickle.load(infile)

       filename_sigma_k_real_variation_0 = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       with open(filename_sigma_k_real_variation_0, 'rb') as infile:
           self_energy_k_real_std_0 = pickle.load(infile)

       filename_sigma_k_imag_0 = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       with open(filename_sigma_k_imag_0, 'rb') as infile:
           self_energy_k_imag_0 = pickle.load(infile)

       filename_sigma_k_imag_variation_0 = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       with open(filename_sigma_k_imag_variation_0, 'rb') as infile:
           self_energy_k_imag_std_0 = pickle.load(infile)


       Self_energy_k_real = self_energy_k_real_0.copy()
       Self_energy_k_real_var = np.power(self_energy_k_real_std_0,2)

       Self_energy_k_imag = self_energy_k_imag_0.copy()
       Self_energy_k_imag_var = np.power(self_energy_k_imag_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            
            filename_sigma_k_real = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_sigma_k_real):

               r_counter=r_counter+1
               with open(filename_sigma_k_real, 'rb') as infile:
                   self_energy_k_real = pickle.load(infile)

               filename_sigma_k_real_variation = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_sigma_k_real_variation, 'rb') as infile:
                   self_energy_k_real_std = pickle.load(infile)

               filename_sigma_k_imag = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_sigma_k_imag, 'rb') as infile:
                  self_energy_k_imag = pickle.load(infile)

               filename_sigma_k_imag_variation = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_sigma_k_imag_variation, 'rb') as infile:
                   self_energy_k_imag_std = pickle.load(infile)
                   
                   
               Self_energy_X_k_real = np.add(Self_energy_k_real,self_energy_k_real)
               Self_energy_k_real = Self_energy_X_k_real.copy()

               Self_energy_Y_k_real_var = np.power(self_energy_k_real_std,2)
               Self_energy_X_k_real_var = np.add(Self_energy_k_real_var,Self_energy_Y_k_real_var)
               Self_energy_k_real_var = Self_energy_X_k_real_var.copy()

               Self_energy_X_k_imag = np.add(Self_energy_k_imag,self_energy_k_imag)
               Self_energy_k_imag = Self_energy_X_k_imag.copy()

               Self_energy_Y_k_imag_var = np.power(self_energy_k_imag_std,2)
               Self_energy_X_k_imag_var = np.add(Self_energy_k_imag_var,Self_energy_Y_k_imag_var)
               Self_energy_k_imag_var = Self_energy_X_k_imag_var.copy()
               

       Self_energy_k_real_avg = Self_energy_k_real/r_counter

       Self_energy_k_real_std_avg = np.sqrt(Self_energy_k_real_var/(r_counter*r_counter))

       Self_energy_k_imag_avg = Self_energy_k_imag/r_counter

       Self_energy_k_imag_std_avg = np.sqrt(Self_energy_k_imag_var/(r_counter*r_counter))
       
       #=======================================Saving real part=======================================================================
       
       filename_sigma_k_real_avg = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       data_sigma_k_real_avg = Self_energy_k_real_avg
       with open(filename_sigma_k_real_avg, 'wb') as outfile:
           pickle.dump(data_sigma_k_real_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_sigma_k_real_std_avg = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       data_sigma_k_real_std_avg = Self_energy_k_real_std_avg
       with open(filename_sigma_k_real_std_avg, 'wb') as outfile:
           pickle.dump(data_sigma_k_real_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_sigma_k_real_text_avg = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       filename_sigma_k_real_variation_text_avg = '%s/Self_energy_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_sigma_k_real_text_avg,Self_energy_k_real_avg)
       np.savetxt(filename_sigma_k_real_variation_text_avg,Self_energy_k_real_std_avg)
       
       #=====================================Saving imaginary part =====================================================================
       
       filename_sigma_k_imag_avg = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       data_sigma_k_imag_avg = Self_energy_k_imag_avg
       with open(filename_sigma_k_imag_avg, 'wb') as outfile:
           pickle.dump(data_sigma_k_imag_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_sigma_k_imag_std_avg = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       data_sigma_k_imag_std_avg = Self_energy_k_imag_std_avg
       with open(filename_sigma_k_imag_std_avg, 'wb') as outfile:
           pickle.dump(data_sigma_k_imag_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_sigma_k_imag_text_avg = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)
       filename_sigma_k_imag_variation_text_avg = '%s/Self_energy_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_sigma_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_sigma_k_imag_text_avg,Self_energy_k_imag_avg)
       np.savetxt(filename_sigma_k_imag_variation_text_avg,Self_energy_k_imag_std_avg)
       
       
       #======================================Saving data for analytic continuation ====================================================
       print("run total",run_counter)
       beta = float(dtau)*float(L)
       
       F_bz,k_points = k_point_grid_upper_half_bz(N)
       print("No of k points", k_points)
       print("Shape of Sigma_k_avg",Self_energy_k_real_avg.shape)
       Omega_n = np.zeros(int(int(L)/2+1))
       
       for n in range(len(Omega_n)):
           Omega_n[n] = 2*np.pi*(n+1)/beta
    
       Norm_self_energy = np.zeros(k_points)
       
       for kk in range(k_points):
           kx = F_bz[kk,1]
           ky = F_bz[kk,2]
           
           Sigma_k_real_data = Self_energy_k_real_avg[:,kk]
           Sigma_k_real_std_data = Self_energy_k_real_std_avg[:,kk]
           Sigma_k_imag_data = Self_energy_k_imag_avg[:,kk]
           Sigma_k_imag_std_data = Self_energy_k_imag_std_avg[:,kk]
           
           Sigma_data_real = np.stack((Omega_n,Sigma_k_real_data,Sigma_k_real_std_data),axis = 1)
           Sigma_data_imag = np.stack((Omega_n,Sigma_k_imag_data,Sigma_k_imag_std_data),axis = 1)
           
           filename_sigma_k_point_real = '%s/Self_energy_real_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_sigma_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_sigma_k_point_real,Sigma_data_real)
           filename_sigma_k_point_imag = '%s/Self_energy_imaginary_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_sigma_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_sigma_k_point_imag,Sigma_data_imag)
           
       print("run total",run_counter)
   else:
       print("Error, momentum space green function file not found")




def matsubara_frequency_chi_xx_correlation_function_momentum_space_average(Text_dir_chi_xx_corr_k_mf,Text_dir_chi_xx_corr_k_mf_avg_ac,N,U,Mu,dtau,L,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_chi_xx_k_real_0 = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)

   if os.path.exists(filename_chi_xx_k_real_0):

       with open(filename_chi_xx_k_real_0, 'rb') as infile:
           chi_xx_correlation_k_real_0 = pickle.load(infile)

       filename_chi_xx_k_real_variation_0 = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       with open(filename_chi_xx_k_real_variation_0, 'rb') as infile:
           chi_xx_correlation_k_real_std_0 = pickle.load(infile)

       filename_chi_xx_k_imag_0 = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       with open(filename_chi_xx_k_imag_0, 'rb') as infile:
           chi_xx_correlation_k_imag_0 = pickle.load(infile)

       filename_chi_xx_k_imag_variation_0 = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       with open(filename_chi_xx_k_imag_variation_0, 'rb') as infile:
           chi_xx_correlation_k_imag_std_0 = pickle.load(infile)


       Chi_xx_correlation_k_real = chi_xx_correlation_k_real_0.copy()
       Chi_xx_correlation_k_real_var = np.power(chi_xx_correlation_k_real_std_0,2)

       Chi_xx_correlation_k_imag = chi_xx_correlation_k_imag_0.copy()
       Chi_xx_correlation_k_imag_var = np.power(chi_xx_correlation_k_imag_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            
            filename_chi_xx_k_real = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_chi_xx_k_real):

               r_counter=r_counter+1
               with open(filename_chi_xx_k_real, 'rb') as infile:
                   chi_xx_correlation_k_real = pickle.load(infile)

               filename_chi_xx_k_real_variation = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_chi_xx_k_real_variation, 'rb') as infile:
                   chi_xx_correlation_k_real_std = pickle.load(infile)

               filename_chi_xx_k_imag = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_chi_xx_k_imag, 'rb') as infile:
                   chi_xx_correlation_k_imag = pickle.load(infile)

               filename_chi_xx_k_imag_variation = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_chi_xx_k_imag_variation, 'rb') as infile:
                   chi_xx_correlation_k_imag_std = pickle.load(infile)
                   
                   
               Chi_xx_correlation_X_k_real = np.add(Chi_xx_correlation_k_real,chi_xx_correlation_k_real)
               Chi_xx_correlation_k_real = Chi_xx_correlation_X_k_real.copy()

               Chi_xx_correlation_Y_k_real_var = np.power(chi_xx_correlation_k_real_std,2)
               Chi_xx_correlation_X_k_real_var = np.add(Chi_xx_correlation_k_real_var,Chi_xx_correlation_Y_k_real_var)
               Chi_xx_correlation_k_real_var = Chi_xx_correlation_X_k_real_var.copy()

               Chi_xx_correlation_X_k_imag = np.add(Chi_xx_correlation_k_imag,chi_xx_correlation_k_imag)
               Chi_xx_correlation_k_imag = Chi_xx_correlation_X_k_imag.copy()

               Chi_xx_correlation_Y_k_imag_var = np.power(chi_xx_correlation_k_imag_std,2)
               Chi_xx_correlation_X_k_imag_var = np.add(Chi_xx_correlation_k_imag_var,Chi_xx_correlation_Y_k_imag_var)
               Chi_xx_correlation_k_imag_var = Chi_xx_correlation_X_k_imag_var.copy()
               

       Chi_xx_correlation_k_real_avg = Chi_xx_correlation_k_real/r_counter

       Chi_xx_correlation_k_real_std_avg = np.sqrt(Chi_xx_correlation_k_real_var/(r_counter*r_counter))

       Chi_xx_correlation_k_imag_avg = Chi_xx_correlation_k_imag/r_counter

       Chi_xx_correlation_k_imag_std_avg = np.sqrt(Chi_xx_correlation_k_imag_var/(r_counter*r_counter))
       
       #=======================================Saving real part=======================================================================
       
       filename_chi_xx_k_real_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       data_chi_xx_k_real_avg = Chi_xx_correlation_k_real_avg
       with open(filename_chi_xx_k_real_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_k_real_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_chi_xx_k_real_std_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       data_chi_xx_k_real_std_avg = Chi_xx_correlation_k_real_std_avg
       with open(filename_chi_xx_k_real_std_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_k_real_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_chi_xx_k_real_text_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       filename_chi_xx_k_real_variation_text_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_chi_xx_k_real_text_avg,Chi_xx_correlation_k_real_avg)
       np.savetxt(filename_chi_xx_k_real_variation_text_avg,Chi_xx_correlation_k_real_std_avg)
       
       #=====================================Saving imaginary part =====================================================================
       
       filename_chi_xx_k_imag_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       data_chi_xx_k_imag_avg = Chi_xx_correlation_k_imag_avg
       with open(filename_chi_xx_k_imag_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_k_imag_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_chi_xx_k_imag_std_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       data_chi_xx_k_imag_std_avg = Chi_xx_correlation_k_imag_std_avg
       with open(filename_chi_xx_k_imag_std_avg, 'wb') as outfile:
           pickle.dump(data_chi_xx_k_imag_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_chi_xx_k_imag_text_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)
       filename_chi_xx_k_imag_variation_text_avg = '%s/Spin_xx_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_chi_xx_corr_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_chi_xx_k_imag_text_avg,Chi_xx_correlation_k_imag_avg)
       np.savetxt(filename_chi_xx_k_imag_variation_text_avg,Chi_xx_correlation_k_imag_std_avg)
       
       
       #======================================Saving data for analytic continuation ====================================================
       print("run total",run_counter)
       beta = float(dtau)*float(L)
       
       F_bz,k_points = k_point_grid_upper_half_bz(N)
       print("No of k points", k_points)
       print("Shape of Chi_xx_k_avg",Chi_xx_correlation_k_real_avg.shape)
       Omega_n = np.zeros(int(int(L)/2+1))
       
       for n in range(len(Omega_n)):
           Omega_n[n] = 2*np.pi*n/beta
           
       for kk in range(k_points):
           kx = F_bz[kk,1]
           ky = F_bz[kk,2]
 
           #========== Note that here we save the correlation function by removing the - sign in the definition, as it will cause problems in MaxEnt ========================
           
           Chi_xx_k_real_data = -1*Chi_xx_correlation_k_real_avg[:,kk]
           Chi_xx_k_real_std_data = Chi_xx_correlation_k_real_std_avg[:,kk]
           Chi_xx_k_imag_data = -1*Chi_xx_correlation_k_imag_avg[:,kk]
           Chi_xx_k_imag_std_data = Chi_xx_correlation_k_imag_std_avg[:,kk]
           

           Chi_xx_data_real = np.stack((Omega_n,Chi_xx_k_real_data,Chi_xx_k_real_std_data),axis = 1)
           filename_chi_xx_k_point_real = '%s/Spin_xx_correlation_function_real_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_chi_xx_corr_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_chi_xx_k_point_real,Chi_xx_data_real)
           
           Chi_xx_data_imag = np.stack((Omega_n,Chi_xx_k_imag_data,Chi_xx_k_imag_std_data),axis = 1)
           filename_chi_xx_k_point_imag = '%s/Spin_xx_correlation_function_imaginary_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_chi_xx_corr_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_chi_xx_k_point_imag,Chi_xx_data_imag)
           
           Chi_xx_data = np.stack((Omega_n,Chi_xx_k_real_data,Chi_xx_k_real_std_data,Chi_xx_k_imag_data,Chi_xx_k_imag_std_data),axis = 1)
           filename_chi_xx_k_point = '%s/Spin_xx_correlation_function_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_chi_xx_corr_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_chi_xx_k_point,Chi_xx_data)
              
       print("run total",run_counter)
   else:
       print("Error, momentum space chi xx correlation function file not found")




def matsubara_frequency_current_correlation_function_momentum_space_average(Text_dir_curr_curr_corr_k_mf,Text_dir_curr_curr_corr_k_mf_avg_ac,N,U,Mu,dtau,L,run_no):

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   filename_curr_k_real_0 = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)

   if os.path.exists(filename_curr_k_real_0):

       with open(filename_curr_k_real_0, 'rb') as infile:
           current_correlation_k_real_0 = pickle.load(infile)

       filename_current_k_real_variation_0 = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       with open(filename_current_k_real_variation_0, 'rb') as infile:
           current_correlation_k_real_std_0 = pickle.load(infile)

       filename_curr_k_imag_0 = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       with open(filename_curr_k_imag_0, 'rb') as infile:
           current_correlation_k_imag_0 = pickle.load(infile)

       filename_current_k_imag_variation_0 = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       with open(filename_current_k_imag_variation_0, 'rb') as infile:
           current_correlation_k_imag_std_0 = pickle.load(infile)


       Current_correlation_k_real = current_correlation_k_real_0.copy()
       Current_correlation_k_real_var = np.power(current_correlation_k_real_std_0,2)

       Current_correlation_k_imag = current_correlation_k_imag_0.copy()
       Current_correlation_k_imag_var = np.power(current_correlation_k_imag_std_0,2)


       run_counter = 1
       r_counter = 1
       while(run_counter<run_no):
            realization = str(run_counter)
       
            run_counter=run_counter+1
            
            filename_curr_k_real = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_curr_k_real):

               r_counter=r_counter+1
               with open(filename_curr_k_real, 'rb') as infile:
                   current_correlation_k_real = pickle.load(infile)

               filename_current_k_real_variation = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_current_k_real_variation, 'rb') as infile:
                   current_correlation_k_real_std = pickle.load(infile)

               filename_curr_k_imag = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_curr_k_imag, 'rb') as infile:
                   current_correlation_k_imag = pickle.load(infile)

               filename_current_k_imag_variation = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L,realization)
               with open(filename_current_k_imag_variation, 'rb') as infile:
                   current_correlation_k_imag_std = pickle.load(infile)
                   
                   
               Current_correlation_X_k_real = np.add(Current_correlation_k_real,current_correlation_k_real)
               Current_correlation_k_real = Current_correlation_X_k_real.copy()

               Current_correlation_Y_k_real_var = np.power(current_correlation_k_real_std,2)
               Current_correlation_X_k_real_var = np.add(Current_correlation_k_real_var,Current_correlation_Y_k_real_var)
               Current_correlation_k_real_var = Current_correlation_X_k_real_var.copy()

               Current_correlation_X_k_imag = np.add(Current_correlation_k_imag,current_correlation_k_imag)
               Current_correlation_k_imag = Current_correlation_X_k_imag.copy()

               Current_correlation_Y_k_imag_var = np.power(current_correlation_k_imag_std,2)
               Current_correlation_X_k_imag_var = np.add(Current_correlation_k_imag_var,Current_correlation_Y_k_imag_var)
               Current_correlation_k_imag_var = Current_correlation_X_k_imag_var.copy()
               

       Current_correlation_k_real_avg = Current_correlation_k_real/r_counter

       Current_correlation_k_real_std_avg = np.sqrt(Current_correlation_k_real_var/(r_counter*r_counter))

       Current_correlation_k_imag_avg = Current_correlation_k_imag/r_counter

       Current_correlation_k_imag_std_avg = np.sqrt(Current_correlation_k_imag_var/(r_counter*r_counter))
       
       #=======================================Saving real part=======================================================================
       
       filename_current_k_real_avg = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       data_current_k_real_avg = Current_correlation_k_real_avg
       with open(filename_current_k_real_avg, 'wb') as outfile:
           pickle.dump(data_current_k_real_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_current_k_real_std_avg = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       data_current_k_real_std_avg = Current_correlation_k_real_std_avg
       with open(filename_current_k_real_std_avg, 'wb') as outfile:
           pickle.dump(data_current_k_real_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_current_k_real_text_avg = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       filename_current_k_real_variation_text_avg = '%s/Current_current_correlation_function_matsubara_frequency_real_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_current_k_real_text_avg,Current_correlation_k_real_avg)
       np.savetxt(filename_current_k_real_variation_text_avg,Current_correlation_k_real_std_avg)
       
       #=====================================Saving imaginary part =====================================================================
       
       filename_current_k_imag_avg = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       data_current_k_imag_avg = Current_correlation_k_imag_avg
       with open(filename_current_k_imag_avg, 'wb') as outfile:
           pickle.dump(data_current_k_imag_avg, outfile, pickle.HIGHEST_PROTOCOL)


       filename_current_k_imag_std_avg = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       data_current_k_imag_std_avg = Current_correlation_k_imag_std_avg
       with open(filename_current_k_imag_std_avg, 'wb') as outfile:
           pickle.dump(data_current_k_imag_std_avg, outfile, pickle.HIGHEST_PROTOCOL)

       filename_current_k_imag_text_avg = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)
       filename_current_k_imag_variation_text_avg = '%s/Current_current_correlation_function_matsubara_frequency_imaginary_part_momentum_space_standard_deviation_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_curr_curr_corr_k_mf,N,U,Mu,dtau,L)

       np.savetxt(filename_current_k_imag_text_avg,Current_correlation_k_imag_avg)
       np.savetxt(filename_current_k_imag_variation_text_avg,Current_correlation_k_imag_std_avg)
       
       
       #======================================Saving data for analytic continuation ====================================================
       print("run total",run_counter)
       beta = float(dtau)*float(L)
       
       F_bz,k_points = k_point_grid_full_bz(N)
       print("No of k points", k_points)
       print("Shape of Curr_k_avg",Current_correlation_k_real_avg.shape)
       Omega_n = np.zeros(int(int(L)/2+1))
       
       for n in range(len(Omega_n)):
           Omega_n[n] = 2*np.pi*n/beta
           
       for kk in range(k_points):
           kx = F_bz[kk,1]
           ky = F_bz[kk,2]
 
           #========== Note that here we save the correlation function by removing the - sign in the definition, as it will cause problems in MaxEnt ========================
           
           Curr_k_real_data = -1*Current_correlation_k_real_avg[:,kk]
           Curr_k_real_std_data = Current_correlation_k_real_std_avg[:,kk]
           Curr_k_imag_data = -1*Current_correlation_k_imag_avg[:,kk]
           Curr_k_imag_std_data = Current_correlation_k_imag_std_avg[:,kk]
           

           Curr_data_real = np.stack((Omega_n,Curr_k_real_data,Curr_k_real_std_data),axis = 1)
           filename_curr_k_point_real = '%s/Current_current_correlation_function_real_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_curr_curr_corr_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_curr_k_point_real,Curr_data_real)
           
           Curr_data_imag = np.stack((Omega_n,Curr_k_imag_data,Curr_k_imag_std_data),axis = 1)
           filename_curr_k_point_imag = '%s/Current_current_correlation_function_imaginary_part_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_curr_curr_corr_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_curr_k_point_imag,Curr_data_imag)
           
           Curr_data = np.stack((Omega_n,Curr_k_real_data,Curr_k_real_std_data,Curr_k_imag_data,Curr_k_imag_std_data),axis = 1)
           filename_curr_k_point = '%s/Current_current_correlation_function_matsubara_frequency_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat' %(Text_dir_curr_curr_corr_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk))
           np.savetxt(filename_curr_k_point,Curr_data)
              
       print("U,mu,L,run total: ",U,Mu,L,run_counter)
   else:
       print("Error, momentum space current correlation function file not found")




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



    #==========================Diretories for real space correlation averages =============================================

    Text_dir_gf_r = '../../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_den_corr_r = '../../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)
    
    Text_dir_chi_xx_corr_r = '../../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_xx_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_chi_zz_corr_r = '../../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_zz_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_curr_curr_corr_r = '../../../Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)



    #==================== Directories for momentum space correlation averages ===================================================

    Text_dir_gf_k = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)
    Text_dir_gf_k_mf = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space_matsubara_frequency'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_sigma_k_mf = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Self_energy_momentum_space_matsubara_frequency'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_den_corr_k = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_chi_xx_corr_k = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_xx_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)
    Text_dir_chi_xx_corr_k_mf = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_xx_correlation_functions_momentum_space_matsubara_frequency'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_chi_zz_corr_k = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_zz_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_curr_curr_corr_k = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_curr_curr_corr_k_mf = '../../../Text_files/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space_matsubara_frequency'%(N,N,U,Dtau,mu,Dtau,L)




    #===========================Directories for correlations for analytic continuation averages, imaginary time space =========================================

    Text_dir_gf_r_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_real_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_r_avg_ac):
        os.makedirs(Text_dir_gf_r_avg_ac)

    Text_dir_gf_k_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_k_avg_ac):
        os.makedirs(Text_dir_gf_k_avg_ac)
  
    Text_dir_den_k_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation_density_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_correlation_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_den_k_avg_ac):
        os.makedirs(Text_dir_den_k_avg_ac)
   
    Text_dir_chi_xx_k_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation_spin_xx_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_xx_correlation_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_chi_xx_k_avg_ac):
        os.makedirs(Text_dir_chi_xx_k_avg_ac)
        
    Text_dir_chi_zz_k_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation_spin_zz_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_zz_correlation_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_chi_zz_k_avg_ac):
        os.makedirs(Text_dir_chi_zz_k_avg_ac)
        
    Text_dir_curr_k_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation_current_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_curr_k_avg_ac):
        os.makedirs(Text_dir_curr_k_avg_ac)
 
 
     #===========================Directories for correlations for analytic continuation averages, matsubara frequency space =========================================

    Text_dir_gf_k_mf_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_k_mf_avg_ac):
        os.makedirs(Text_dir_gf_k_mf_avg_ac)

    Text_dir_sigma_k_mf_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation_self_energy/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Self_energy_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_sigma_k_mf_avg_ac):
        os.makedirs(Text_dir_sigma_k_mf_avg_ac)

    Text_dir_chi_xx_k_mf_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation_spin_xx_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_xx_correlation_functions_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_chi_xx_k_mf_avg_ac):
        os.makedirs(Text_dir_chi_xx_k_mf_avg_ac)
        
    Text_dir_curr_k_mf_avg_ac = '../../../Text_files/Text_files_N_%s_analytic_continuation_current_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Current_current_correlation_functions_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_curr_k_mf_avg_ac):
        os.makedirs(Text_dir_curr_k_mf_avg_ac)
        
        
    Text_dir_eqm = "../../../Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements"%(N,N,U,Dtau,mu,Dtau,L)
    
    filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U,mu,Dtau,L)
    with open(filename_eqm_avg, 'rb') as infile:
         sys_measure_avg = pickle.load(infile)

    Nden = sys_measure_avg['Density averaged']
    print("avg sign, density, mu", Nden, mu)
    
    
    # ====================================== Unequal time data averaging ==================================================================================
    
    unequal_time_greens_function_real_space_average(Text_dir_gf_r,Text_dir_gf_r_avg_ac,N,U,mu,Dtau,L,run_no)
    unequal_time_greens_function_momentum_space_average(Text_dir_gf_k,Text_dir_gf_k_avg_ac,N,U,mu,Dtau,L,Nden,run_no)

    #unequal_time_density_correlation_function_real_space_average(Text_dir_den_corr_r,Text_dir_den_corr_r_avg,N,U,mu,Dtau,L)
    unequal_time_density_correlation_function_momentum_space_average(Text_dir_den_corr_k,Text_dir_den_k_avg_ac,N,U,mu,Dtau,L,run_no)

    #unequal_time_chi_xx_correlation_function_real_space_average(Text_dir_chi_xx_corr_r,N,U,mu,Dtau,L,run_no)
    unequal_time_chi_xx_correlation_function_momentum_space_average(Text_dir_chi_xx_corr_k,Text_dir_chi_xx_k_avg_ac,N,U,mu,Dtau,L,run_no)

    #unequal_time_chi_zz_correlation_function_real_space_average(Text_dir_chi_zz_corr_r,N,U,mu,Dtau,L,run_no)
    unequal_time_chi_zz_correlation_function_momentum_space_average(Text_dir_chi_zz_corr_k,Text_dir_chi_zz_k_avg_ac,N,U,mu,Dtau,L,run_no)

    #unequal_time_current_correlation_function_real_space_average(Text_dir_curr_curr_corr_r,Text_dir_curr_curr_corr_r_avg,N,U,mu,Dtau,L)
    unequal_time_current_correlation_function_momentum_space_average(Text_dir_curr_curr_corr_k,Text_dir_curr_k_avg_ac,N,U,mu,Dtau,L,run_no)

    # ================================================= Matsubara frequency space data averaging   ================================================================
    
    matsubara_frequency_green_function_momentum_space_average(Text_dir_gf_k_mf,Text_dir_gf_k_mf_avg_ac,N,U,mu,Dtau,L,Nden,run_no)
    matsubara_frequency_self_energy_momentum_space_average(Text_dir_sigma_k_mf,Text_dir_sigma_k_mf_avg_ac,N,U,mu,Dtau,L,run_no)
    matsubara_frequency_chi_xx_correlation_function_momentum_space_average(Text_dir_chi_xx_corr_k_mf,Text_dir_chi_xx_k_mf_avg_ac,N,U,mu,Dtau,L,run_no)
    matsubara_frequency_current_correlation_function_momentum_space_average(Text_dir_curr_curr_corr_k_mf,Text_dir_curr_k_mf_avg_ac,N,U,mu,Dtau,L,run_no)
    print("===================================================================================")


if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)
