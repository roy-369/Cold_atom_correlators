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

def neighbor_locator(x_cord,y_cord,rad_sq):

    idx = []
    for i in range(len(x_cord)):
        rdius = x_cord[i]*x_cord[i]+y_cord[i]*y_cord[i]
        if(round(rdius,2) == round(rad_sq,2)):
          idx.append(i)

    return idx


def neighbor_average(r_rel_grid, corr_array,corr_arr_std):

    R_cord = []
    
    rx = r_rel_grid[:,0]
    ry = r_rel_grid[:,1]

    r_2 = np.add(np.power(rx,2),np.power(ry,2))

    r_2_unq = np.sort(np.unique(r_2))

    Corr_arr_nbr_avg = np.zeros(len(r_2_unq))
    Corr_arr_nbr_std_avg = np.zeros(len(r_2_unq))

    for i in range(len(r_2_unq)):
        idx = neighbor_locator(rx,ry,r_2_unq[i])
        corr_nbr = 0
        corr_nbr_var = 0
        for j in range(len(idx)):
            corr_arr_idx = idx[j]
            corr_nbr = corr_nbr+corr_array[corr_arr_idx]
            corr_nbr_var = corr_nbr_var+(corr_arr_std[corr_arr_idx])**2
        Corr_arr_nbr_avg[i] = corr_nbr/len(idx)
        Corr_arr_nbr_std_avg[i] = np.sqrt(corr_nbr_var/(len(idx)*len(idx)))

    return r_2_unq,Corr_arr_nbr_avg, Corr_arr_nbr_std_avg
    
    
def mat_append(mat_A,mat_B):

    rank_A = np.array(mat_A).ndim

    if(rank_A == 1):
       mat_C = np.zeros((len(mat_A),2))
       mat_C[:,0] = np.copy(mat_A)
       mat_C[:,1] = np.copy(mat_B)
       return mat_C

    if(rank_A == 2):
       mat_C = np.zeros((len(mat_A[:,0]),len(mat_A[0,:])+1))
       for xx in range(len(mat_A[0,:])):
           mat_C[:,xx] = np.copy(mat_A[:,xx])
       mat_C[:,-1] = np.copy(mat_B)
       return mat_C



def normal_average(data_set,data_set_std):

   data_avg = 0
   data_avg_var = 0
   
   a = 0
   
   for i in range(len(data_set)):
       data_avg = data_avg+data_set[i]
       data_avg_var = data_avg_var+data_set_std[i]*data_set_std[i]
       a= a+1
   

   return data_avg/a, np.sqrt(data_avg_var/(a*a))



        

#=========================================================Calcualting normal average over correlations across realizations====================================================================================================


def equal_time_density_density_correlation_function_real_space_normal_average(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,run_no):


   filename_eq_du_du_r_0 = '%s/Equal_time_Density_up_Density_up_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_du_du_r_0):

       with open(filename_eq_du_du_r_0, 'rb') as infile:
            eq_du_du_r_0 = pickle.load(infile)

       filename_eq_du_du_r_variation_0 = "%s/Equal_time_Density_up_Density_up_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl"%(Text_dir_eq_den_den_r,N,U,Mu,dtau,L)
       with open(filename_eq_du_du_r_variation_0, 'rb') as infile:
            eq_du_du_r_std_0 = pickle.load(infile)

       filename_eq_du_dn_r_0 = '%s/Equal_time_Density_up_Density_down_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L)
       with open(filename_eq_du_dn_r_0, 'rb') as infile:
            eq_du_dn_r_0 = pickle.load(infile)

       filename_eq_du_dn_r_variation_0 = '%s/Equal_time_Density_up_Density_down_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L)
       with open(filename_eq_du_dn_r_variation_0, 'rb') as infile:
            eq_du_dn_r_std_0 = pickle.load(infile)


       eq_den_den_r_0 = 2*np.add(eq_du_du_r_0,eq_du_dn_r_0)
       eq_den_den_r_std_0 = 2*np.sqrt(np.add(np.power(eq_du_du_r_std_0,2),np.power(eq_du_dn_r_std_0,2)))

       r_grid, r_rel_grid_size = r_point_upper_half_grid(N)

       r_sep_2,Eq_den_den_r_nbr,Eq_den_den_r_nbr_std = neighbor_average(r_grid, eq_den_den_r_0, eq_den_den_r_std_0)

       run_counter = 0
       r_c = 1

       while(run_counter<int(run_no)):

            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_du_du_r = '%s/Equal_time_Density_up_Density_up_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_du_du_r):

               with open(filename_eq_du_du_r, 'rb') as infile:
                    eq_du_du_r = pickle.load(infile)

               filename_eq_du_du_r_variation = "%s/Equal_time_Density_up_Density_up_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_du_du_r_variation, 'rb') as infile:
                    eq_du_du_r_std = pickle.load(infile)

               filename_eq_du_dn_r = '%s/Equal_time_Density_up_Density_down_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_du_dn_r, 'rb') as infile:
                    eq_du_dn_r = pickle.load(infile)

               filename_eq_du_dn_r_variation = "%s/Equal_time_Density_up_Density_down_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_du_dn_r_variation, 'rb') as infile:
                    eq_du_dn_r_std = pickle.load(infile)

               eq_den_den_r = 2*np.add(eq_du_du_r,eq_du_dn_r)
               eq_den_den_r_std = 2*np.sqrt(np.add(np.power(eq_du_du_r_std,2),np.power(eq_du_dn_r_std,2)))

               r_sep_2, eq_den_den_r_nbr, eq_den_den_r_nbr_std = neighbor_average(r_grid, eq_den_den_r, eq_den_den_r_std)

               Eq_X_Den_Den_r_nbr = mat_append(Eq_den_den_r_nbr,eq_den_den_r_nbr)
               Eq_den_den_r_nbr = np.copy(Eq_X_Den_Den_r_nbr)

               Eq_X_Den_Den_r_nbr_std = mat_append(Eq_den_den_r_nbr_std,eq_den_den_r_nbr_std)
               Eq_den_den_r_nbr_std = np.copy(Eq_X_Den_Den_r_nbr_std)
               r_c = r_c+1

       Eq_den_den_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_den_den_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_den_den_r_nbr_avg[ii],Eq_den_den_r_nbr_std_avg[ii] = normal_average(Eq_den_den_r_nbr[ii,:],Eq_den_den_r_nbr_std[ii,:])


       filename_eq_den_den_r_nbr_avg = '%s/Equal_time_Density_Density_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L)
       data_eq_den_den_r_nbr_avg = np.stack((r_sep_2,Eq_den_den_r_nbr_avg,Eq_den_den_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_den_den_r_nbr_avg,data_eq_den_den_r_nbr_avg)
       print("run total",r_c)
   else:
       print("Error, den dbn file not found")




def equal_time_density_doublon_correlation_function_real_space_normal_average(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,run_no):


   filename_eq_du_dbn_r_0 = '%s/Equal_time_Density_up_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_du_dbn_r_0):

       with open(filename_eq_du_dbn_r_0, 'rb') as infile:
            eq_du_dbn_r_0 = pickle.load(infile)

       filename_eq_du_dbn_r_variation_0 = "%s/Equal_time_Density_up_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl"%(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L)
       with open(filename_eq_du_dbn_r_variation_0, 'rb') as infile:
            eq_du_dbn_r_std_0 = pickle.load(infile)

       filename_eq_dn_dbn_r_0 = '%s/Equal_time_Density_down_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L)
       with open(filename_eq_dn_dbn_r_0, 'rb') as infile:
            eq_dn_dbn_r_0 = pickle.load(infile)

       filename_eq_dn_dbn_r_variation_0 = '%s/Equal_time_Density_down_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L)
       with open(filename_eq_dn_dbn_r_variation_0, 'rb') as infile:
            eq_dn_dbn_r_std_0 = pickle.load(infile)


       eq_den_dbn_r_0 = np.add(eq_du_dbn_r_0,eq_dn_dbn_r_0)
       eq_den_dbn_r_std_0 = np.sqrt(np.add(np.power(eq_du_dbn_r_std_0,2),np.power(eq_dn_dbn_r_std_0,2)))

       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_den_dbn_r_nbr,Eq_den_dbn_r_nbr_std = neighbor_average(r_grid, eq_den_dbn_r_0, eq_den_dbn_r_std_0)

       run_counter = 0
       r_c = 1

       while(run_counter<int(run_no)):

            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_du_dbn_r = '%s/Equal_time_Density_up_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_du_dbn_r):

               with open(filename_eq_du_dbn_r, 'rb') as infile:
                    eq_du_dbn_r = pickle.load(infile)

               filename_eq_du_dbn_r_variation = "%s/Equal_time_Density_up_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_du_dbn_r_variation, 'rb') as infile:
                    eq_du_dbn_r_std = pickle.load(infile)

               filename_eq_dn_dbn_r = '%s/Equal_time_Density_down_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_dn_dbn_r, 'rb') as infile:
                    eq_dn_dbn_r = pickle.load(infile)

               filename_eq_dn_dbn_r_variation = "%s/Equal_time_Density_down_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_dn_dbn_r_variation, 'rb') as infile:
                    eq_dn_dbn_r_std = pickle.load(infile)

               eq_den_dbn_r = np.add(eq_du_dbn_r,eq_dn_dbn_r)
               eq_den_dbn_r_std = np.sqrt(np.add(np.power(eq_du_dbn_r_std,2),np.power(eq_dn_dbn_r_std,2)))

               r_sep_2, eq_den_dbn_r_nbr, eq_den_dbn_r_nbr_std = neighbor_average(r_grid, eq_den_dbn_r, eq_den_dbn_r_std)

               Eq_X_Den_Dbn_r_nbr = mat_append(Eq_den_dbn_r_nbr,eq_den_dbn_r_nbr)
               Eq_den_dbn_r_nbr = np.copy(Eq_X_Den_Dbn_r_nbr)

               Eq_X_Den_Dbn_r_nbr_std = mat_append(Eq_den_dbn_r_nbr_std,eq_den_dbn_r_nbr_std)
               Eq_den_dbn_r_nbr_std = np.copy(Eq_X_Den_Dbn_r_nbr_std)
               r_c = r_c+1

       Eq_den_dbn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_den_dbn_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_den_dbn_r_nbr_avg[ii],Eq_den_dbn_r_nbr_std_avg[ii] = normal_average(Eq_den_dbn_r_nbr[ii,:],Eq_den_dbn_r_nbr_std[ii,:])


       filename_eq_den_dbn_r_nbr_avg = '%s/Equal_time_Density_Doublon_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L)
       data_eq_den_dbn_r_nbr_avg = np.stack((r_sep_2,Eq_den_dbn_r_nbr_avg,Eq_den_dbn_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_den_dbn_r_nbr_avg,data_eq_den_dbn_r_nbr_avg)
       print("run total",r_c)
   else:
       print("Error, den dbn file not found")




def equal_time_doublon_doublon_correlation_function_real_space_normal_average(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L,run_no):


   filename_eq_dbn_dbn_r_0 = '%s/Equal_time_Doublon_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_dbn_dbn_r_0):

       with open(filename_eq_dbn_dbn_r_0, 'rb') as infile:
            eq_dbn_dbn_r_0 = pickle.load(infile)

       filename_eq_dbn_dbn_r_variation_0 = "%s/Equal_time_Doublon_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl"%(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L)
       with open(filename_eq_dbn_dbn_r_variation_0, 'rb') as infile:
            eq_dbn_dbn_r_std_0 = pickle.load(infile)


       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_dbn_dbn_r_nbr,Eq_dbn_dbn_r_nbr_std = neighbor_average(r_grid, eq_dbn_dbn_r_0, eq_dbn_dbn_r_std_0)       

       run_counter = 0
       r_c = 1
       while(run_counter<int(run_no)):

            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_dbn_dbn_r = '%s/Equal_time_Doublon_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_dbn_dbn_r):

               with open(filename_eq_dbn_dbn_r, 'rb') as infile:
                    eq_dbn_dbn_r = pickle.load(infile)

               filename_eq_dbn_dbn_r_variation = "%s/Equal_time_Doublon_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_dbn_dbn_r_variation, 'rb') as infile:
                    eq_dbn_dbn_r_std = pickle.load(infile)
            
               r_sep_2, eq_dbn_dbn_r_nbr, eq_dbn_dbn_r_nbr_std = neighbor_average(r_grid, eq_dbn_dbn_r, eq_dbn_dbn_r_std)

               Eq_X_Dbn_Dbn_r_nbr = mat_append(Eq_dbn_dbn_r_nbr,eq_dbn_dbn_r_nbr)
               Eq_dbn_dbn_r_nbr = np.copy(Eq_X_Dbn_Dbn_r_nbr)

               Eq_X_Dbn_Dbn_r_nbr_std = mat_append(Eq_dbn_dbn_r_nbr_std,eq_dbn_dbn_r_nbr_std)
               Eq_dbn_dbn_r_nbr_std = np.copy(Eq_X_Dbn_Dbn_r_nbr_std)

               r_c = r_c+1

       Eq_dbn_dbn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_dbn_dbn_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_dbn_dbn_r_nbr_avg[ii],Eq_dbn_dbn_r_nbr_std_avg[ii] = normal_average(Eq_dbn_dbn_r_nbr[ii,:],Eq_dbn_dbn_r_nbr_std[ii,:])


       filename_eq_dbn_dbn_r_nbr_avg = '%s/Equal_time_Doublon_Doublon_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L)
       data_eq_dbn_dbn_r_nbr_avg = np.stack((r_sep_2,Eq_dbn_dbn_r_nbr_avg,Eq_dbn_dbn_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_dbn_dbn_r_nbr_avg,data_eq_dbn_dbn_r_nbr_avg)
       print("run total",r_c)
   else:
       print("Error, dbn dbn file not found")




def equal_time_moment_moment_correlation_function_real_space_normal_average(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L,run_no):


   filename_eq_m2_m2_r_0 = '%s/Equal_time_Moment_Moment_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_m2_m2_r_0):

       with open(filename_eq_m2_m2_r_0, 'rb') as infile:
            eq_m2_m2_r_0 = pickle.load(infile)

       filename_eq_m2_m2_r_variation_0 = "%s/Equal_time_Moment_Moment_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl"%(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L)
       with open(filename_eq_m2_m2_r_variation_0, 'rb') as infile:
            eq_m2_m2_r_std_0 = pickle.load(infile)


       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_m2_m2_r_nbr,Eq_m2_m2_r_nbr_std = neighbor_average(r_grid, eq_m2_m2_r_0, eq_m2_m2_r_std_0)

       run_counter = 0
       r_c = 1

       while(run_counter<int(run_no)):

            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_m2_m2_r = '%s/Equal_time_Moment_Moment_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L,realization)
            if os.path.exists(filename_eq_m2_m2_r):
               with open(filename_eq_m2_m2_r, 'rb') as infile:
                    eq_m2_m2_r = pickle.load(infile)

               filename_eq_m2_m2_r_variation = "%s/Equal_time_Moment_Moment_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_m2_m2_r_variation, 'rb') as infile:
                    eq_m2_m2_r_std = pickle.load(infile)

               r_sep_2, eq_m2_m2_r_nbr, eq_m2_m2_r_nbr_std = neighbor_average(r_grid, eq_m2_m2_r, eq_m2_m2_r_std)

               Eq_X_M2_M2_r_nbr = mat_append(Eq_m2_m2_r_nbr,eq_m2_m2_r_nbr)
               Eq_m2_m2_r_nbr = np.copy(Eq_X_M2_M2_r_nbr)

               Eq_X_M2_M2_r_nbr_std = mat_append(Eq_m2_m2_r_nbr_std,eq_m2_m2_r_nbr_std)
               Eq_m2_m2_r_nbr_std = np.copy(Eq_X_M2_M2_r_nbr_std)

               r_c = r_c +1

       Eq_m2_m2_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_m2_m2_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_m2_m2_r_nbr_avg[ii],Eq_m2_m2_r_nbr_std_avg[ii] = normal_average(Eq_m2_m2_r_nbr[ii,:],Eq_m2_m2_r_nbr_std[ii,:])


       filename_eq_m2_m2_r_nbr_avg = '%s/Equal_time_Moment_Moment_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L)
       data_eq_m2_m2_r_nbr_avg = np.stack((r_sep_2,Eq_m2_m2_r_nbr_avg,Eq_m2_m2_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_m2_m2_r_nbr_avg,data_eq_m2_m2_r_nbr_avg)
       print("run total",r_c)
   else:
       print("Error, m2 m2 file not found")


def equal_time_spin_spin_xx_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,run_no):


   filename_eq_sxsx_r_0 = '%s/Equal_time_SxSx_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_sxsx_r_0):

       with open(filename_eq_sxsx_r_0, 'rb') as infile:
            eq_sxsx_r_0 = pickle.load(infile)

       filename_eq_sxsx_r_variation_0 = "%s/Equal_time_SxSx_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl"%(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L)
       with open(filename_eq_sxsx_r_variation_0, 'rb') as infile:
            eq_sxsx_r_std_0 = pickle.load(infile)


       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_sxsx_r_nbr,Eq_sxsx_r_nbr_std = neighbor_average(r_grid, eq_sxsx_r_0, eq_sxsx_r_std_0)

       run_counter = 0
       r_c=1

       while(run_counter<int(run_no)):

            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_sxsx_r = '%s/Equal_time_SxSx_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_sxsx_r):
               with open(filename_eq_sxsx_r, 'rb') as infile:
                    eq_sxsx_r = pickle.load(infile)

               filename_eq_sxsx_r_variation = "%s/Equal_time_SxSx_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_sxsx_r_variation, 'rb') as infile:
                    eq_sxsx_r_std = pickle.load(infile)

               if os.path.exists(filename_eq_sxsx_r):

                  r_sep_2, eq_sxsx_r_nbr, eq_sxsx_r_nbr_std = neighbor_average(r_grid, eq_sxsx_r, eq_sxsx_r_std)

                  Eq_X_SxSx_r_nbr = mat_append(Eq_sxsx_r_nbr,eq_sxsx_r_nbr)
                  Eq_sxsx_r_nbr = np.copy(Eq_X_SxSx_r_nbr)

                  Eq_X_SxSx_r_nbr_std = mat_append(Eq_sxsx_r_nbr_std,eq_sxsx_r_nbr_std)
                  Eq_sxsx_r_nbr_std = np.copy(Eq_X_SxSx_r_nbr_std)
                  r_c = r_c+1

       Eq_sxsx_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_sxsx_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_sxsx_r_nbr_avg[ii],Eq_sxsx_r_nbr_std_avg[ii] = normal_average(Eq_sxsx_r_nbr[ii,:],Eq_sxsx_r_nbr_std[ii,:])


       filename_eq_sxsx_r_nbr_avg = '%s/Equal_time_SxSx_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L)
       data_eq_sxsx_r_nbr_avg = np.stack((r_sep_2,Eq_sxsx_r_nbr_avg,Eq_sxsx_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_sxsx_r_nbr_avg,data_eq_sxsx_r_nbr_avg)
       print("run total",r_c)
   else:
       print("Error, sx sx file not found")       



def equal_time_spin_spin_zz_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,run_no):


   filename_eq_szsz_r_0 = '%s/Equal_time_SzSz_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_szsz_r_0):

       with open(filename_eq_szsz_r_0, 'rb') as infile:
            eq_szsz_r_0 = pickle.load(infile)

       filename_eq_szsz_r_variation_0 = "%s/Equal_time_SzSz_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl"%(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L)
       with open(filename_eq_szsz_r_variation_0, 'rb') as infile:
            eq_szsz_r_std_0 = pickle.load(infile)


       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_szsz_r_nbr,Eq_szsz_r_nbr_std = neighbor_average(r_grid, eq_szsz_r_0, eq_szsz_r_std_0)

       run_counter = 0
       r_c=1

       while(run_counter<int(run_no)):

            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_szsz_r = '%s/Equal_time_SzSz_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_szsz_r):
               with open(filename_eq_szsz_r, 'rb') as infile:
                    eq_szsz_r = pickle.load(infile)

               filename_eq_szsz_r_variation = "%s/Equal_time_SzSz_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,realization)
               with open(filename_eq_szsz_r_variation, 'rb') as infile:
                    eq_szsz_r_std = pickle.load(infile)

               if os.path.exists(filename_eq_szsz_r):

                  r_sep_2, eq_szsz_r_nbr, eq_szsz_r_nbr_std = neighbor_average(r_grid, eq_szsz_r, eq_szsz_r_std)

                  Eq_X_SzSz_r_nbr = mat_append(Eq_szsz_r_nbr,eq_szsz_r_nbr)
                  Eq_szsz_r_nbr = np.copy(Eq_X_SzSz_r_nbr)

                  Eq_X_SzSz_r_nbr_std = mat_append(Eq_szsz_r_nbr_std,eq_szsz_r_nbr_std)
                  Eq_szsz_r_nbr_std = np.copy(Eq_X_SzSz_r_nbr_std)
                  r_c=r_c+1 
  
       Eq_szsz_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_szsz_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_szsz_r_nbr_avg[ii],Eq_szsz_r_nbr_std_avg[ii] = normal_average(Eq_szsz_r_nbr[ii,:],Eq_szsz_r_nbr_std[ii,:])


       filename_eq_szsz_r_nbr_avg = '%s/Equal_time_SzSz_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L)
       data_eq_szsz_r_nbr_avg = np.stack((r_sep_2,Eq_szsz_r_nbr_avg,Eq_szsz_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_szsz_r_nbr_avg,data_eq_szsz_r_nbr_avg)
       print("run total",r_c)
   else:
       print("Error, sz sz file not found")


def main(total,cmdargs):
    if(total!=7):
        raise ValueError('missing args')

    N = cmdargs[1]
    U = cmdargs[2]
    mu = cmdargs[3]
    L = cmdargs[4]
    Dtau = cmdargs[5]
    N_runs = int(cmdargs[6])
    print(U,"U")
    print(L,"L")
    Beta = str((float(L))*(float(Dtau)))
    print(Beta,"Beta")


#===============================================================Location of data files=============================================================================================================================================
    
    Text_dir_eqm = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(N,N,U,Dtau,mu,Dtau,L)
       
    Text_dir_eq_den_den_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_density_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_eq_dbn_dbn_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Doublon_doublon_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_eq_den_dbn_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_doublon_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_eq_m2_m2_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Moment_moment_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_eq_spin_spin_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_spin_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)



#==========================================================================================


    print("avg_running")
    equal_time_density_density_correlation_function_real_space_normal_average(Text_dir_eq_den_den_r,N,U,mu,Dtau,L,N_runs)
    equal_time_doublon_doublon_correlation_function_real_space_normal_average(Text_dir_eq_dbn_dbn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_moment_moment_correlation_function_real_space_normal_average(Text_dir_eq_m2_m2_r,N,U,mu,Dtau,L,N_runs)
    equal_time_density_doublon_correlation_function_real_space_normal_average(Text_dir_eq_den_dbn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_spin_spin_xx_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_r,N,U,mu,Dtau,L,N_runs)
    equal_time_spin_spin_zz_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_r,N,U,mu,Dtau,L,N_runs)
    print("done,", U,mu,L)

if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)
