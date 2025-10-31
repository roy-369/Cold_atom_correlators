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
#===========================================================================================Calculating connected correlators for different realizations===================================================================================



def equal_time_density_density_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_den_den_r,Text_dir_eq_den_den_conn_r,N,U,Mu,dtau,L,run_no):


    run_counter = 0
    while(run_counter<int(run_no)):

        Realization = str(run_counter)

        filename_eqm = '%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eqm,N,U,Mu,dtau,L,Realization)
        filename_eq_dudu_r = '%s/Equal_time_Density_up_Density_up_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
 
        if os.path.exists(filename_eq_dudu_r): # and os.path.exists(filename_eqm):

            with open(filename_eqm, 'rb') as infile:
                sys_measure = pickle.load(infile)

            with open(filename_eq_dudu_r, 'rb') as infile:
                eq_dudu_r = pickle.load(infile)
           
            if bool(sys_measure):
               density_st = sys_measure['Average density']
               density = density_st.split(' ')
               num_den = float(density[0].strip(' '))
               num_den_dev = float(density[8].strip(' '))

            #print(num_den,num_den_dev,"num_den_meas")

               filename_eq_dudu_r_variation = '%s/Equal_time_Density_up_Density_up_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
               with open(filename_eq_dudu_r_variation, 'rb') as infile:
                    eq_dudu_r_std = pickle.load(infile)

               filename_eq_dudn_r = '%s/Equal_time_Density_up_Density_down_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
               with open(filename_eq_dudn_r, 'rb') as infile:
                    eq_dudn_r = pickle.load(infile)

               filename_eq_dudn_r_variation = '%s/Equal_time_Density_up_Density_down_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_den_r,N,U,Mu,dtau,L,Realization)
               with open(filename_eq_dudn_r_variation, 'rb') as infile:
                    eq_dudn_r_std = pickle.load(infile)

               eq_den_den_r = 2*np.add(eq_dudu_r,eq_dudn_r)
               eq_den_den_r_var = 4*np.add(np.power(eq_dudu_r_std,2),np.power(eq_dudn_r_std,2))
       
               eq_den_den_r_connected = np.zeros(len(eq_den_den_r))
               eq_den_den_r_connected_std = np.zeros(len(eq_den_den_r))

               for ii in range(len(eq_den_den_r)):
                   eq_den_den_r_connected[ii] = eq_den_den_r[ii]-num_den*num_den
                   eq_den_den_r_connected_var = eq_den_den_r_var[ii]+(2*num_den*num_den_dev)**2 
                   eq_den_den_r_connected_std[ii] = np.sqrt(eq_den_den_r_connected_var)

               filename_eq_den_den_r_conn = "%s/Equal_time_Density_Density_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_den_den_conn_r,N,U,Mu,dtau,L,Realization)
               filename_eq_den_den_r_conn_txt = "%s/Equal_time_Density_Density_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(Text_dir_eq_den_den_conn_r,N,U,Mu,dtau,L,Realization)
               Eq_den_den_conn_corr_data =  np.stack((eq_den_den_r_connected,eq_den_den_r_connected_std),axis = 1)

               with open(filename_eq_den_den_r_conn, 'wb') as outfile:
                  pickle.dump(Eq_den_den_conn_corr_data, outfile, pickle.HIGHEST_PROTOCOL)

               np.savetxt(filename_eq_den_den_r_conn_txt,Eq_den_den_conn_corr_data)
               print("uff kaka")   
        run_counter=run_counter+1



def equal_time_doublon_doublon_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_dbn_dbn_r,Text_dir_eq_dbn_dbn_conn_r,N,U,Mu,dtau,L,run_no):
    
    run_counter = 0
    while(run_counter<int(run_no)):
         
        Realization = str(run_counter)
        
        filename_eqm = "%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eqm,N,U,Mu,dtau,L,Realization)
        filename_eq_dbn_dbn_r = "%s/Equal_time_Doublon_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L,Realization)

        if os.path.exists(filename_eq_dbn_dbn_r): # and os.path.exists(filename_eqm):

             with open(filename_eqm, 'rb') as infile:
                  sys_measure = pickle.load(infile)

             with open(filename_eq_dbn_dbn_r, 'rb') as infile:
                  eq_dbn_dbn_r = pickle.load(infile)
             
             if bool(sys_measure):
                
                nup_ndn_st = sys_measure['Average Nup*Ndn']
                nup_ndn = nup_ndn_st.split(' ')
                Dbn_occ = float(nup_ndn[0].strip(' '))
                Dbn_occ_dev = float(nup_ndn[8].strip(' '))
             #print(Dbn_occ,Dbn_occ_dev,"Doublon_no")

                filename_eq_dbn_dbn_r_variation = "%s/Equal_time_Doublon_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_dbn_dbn_r,N,U,Mu,dtau,L,Realization)
                with open(filename_eq_dbn_dbn_r_variation, 'rb') as infile:
                     eq_dbn_dbn_r_std = pickle.load(infile)

                eq_dbn_dbn_r_var = np.power(eq_dbn_dbn_r_std,2)

                eq_dbn_dbn_r_connected = np.zeros(len(eq_dbn_dbn_r))
                eq_dbn_dbn_r_connected_std = np.zeros(len(eq_dbn_dbn_r))

                for ii in range(len(eq_dbn_dbn_r)):
                    
                    eq_dbn_dbn_r_connected[ii] = eq_dbn_dbn_r[ii]-Dbn_occ*Dbn_occ
                    eq_dbn_dbn_r_connected_var = eq_dbn_dbn_r_var[ii]+(2*Dbn_occ*Dbn_occ_dev)**2
                    eq_dbn_dbn_r_connected_std[ii] = np.sqrt(eq_dbn_dbn_r_connected_var)


                filename_eq_dbn_dbn_r_conn = "%s/Equal_time_Doublon_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_dbn_dbn_conn_r,N,U,Mu,dtau,L,Realization)            
                filename_eq_dbn_dbn_r_conn_txt = "%s/Equal_time_Doublon_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(Text_dir_eq_dbn_dbn_conn_r,N,U,Mu,dtau,L,Realization)

                Eq_dbn_dbn_conn_corr_data =  np.stack((eq_dbn_dbn_r_connected,eq_dbn_dbn_r_connected_std),axis = 1)

                with open(filename_eq_dbn_dbn_r_conn, 'wb') as outfile:
                   pickle.dump(Eq_dbn_dbn_conn_corr_data, outfile, pickle.HIGHEST_PROTOCOL)
               
                np.savetxt(filename_eq_dbn_dbn_r_conn_txt,Eq_dbn_dbn_conn_corr_data)
                print("uff kaka")

        run_counter=run_counter+1


def equal_time_density_doublon_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_den_dbn_r,Text_dir_eq_den_dbn_conn_r,N,U,Mu,dtau,L,run_no):


    run_counter = 0
    while(run_counter<int(run_no)):

        Realization = str(run_counter)

        filename_eqm = '%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eqm,N,U,Mu,dtau,L,Realization)
        filename_eq_du_dbn_r = '%s/Equal_time_Density_up_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)

        if os.path.exists(filename_eq_du_dbn_r) and os.path.exists(filename_eqm):

            with open(filename_eqm, 'rb') as infile:
                sys_measure = pickle.load(infile)

            with open(filename_eq_du_dbn_r, 'rb') as infile:
                eq_du_dbn_r = pickle.load(infile)

            if bool(sys_measure):

               density_st = sys_measure['Average density']
               density = density_st.split(' ')
               num_den = float(density[0].strip(' '))
               num_den_dev = float(density[8].strip(' '))
               nup_ndn_st = sys_measure['Average Nup*Ndn']
               nup_ndn = nup_ndn_st.split(' ')
               Dbn_occ = float(nup_ndn[0].strip(' '))
               Dbn_occ_dev = float(nup_ndn[8].strip(' '))
            #print(Dbn_occ,Dbn_occ_dev,"Doublon_no")
            #print(num_den,num_den_dev,"num_den_meas")


               filename_eq_du_dbn_r_variation = '%s/Equal_time_Density_up_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)
               with open(filename_eq_du_dbn_r_variation, 'rb') as infile:
                   eq_du_dbn_r_std = pickle.load(infile)

               filename_eq_dn_dbn_r = '%s/Equal_time_Density_down_Doublon_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)
               with open(filename_eq_dn_dbn_r, 'rb') as infile:
                    eq_dn_dbn_r = pickle.load(infile)

               filename_eq_dn_dbn_r_variation = '%s/Equal_time_Density_down_Doublon_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_den_dbn_r,N,U,Mu,dtau,L,Realization)
               with open(filename_eq_du_dbn_r_variation, 'rb') as infile:
                    eq_dn_dbn_r_std = pickle.load(infile)


               eq_den_dbn_r = np.add(eq_du_dbn_r,eq_dn_dbn_r)
               eq_den_dbn_r_var = np.add(np.power(eq_du_dbn_r_std,2),np.power(eq_dn_dbn_r_std,2))

               eq_den_dbn_r_connected = np.zeros(len(eq_den_dbn_r))
               eq_den_dbn_r_connected_std = np.zeros(len(eq_den_dbn_r))

               for ii in range(len(eq_den_dbn_r)):
                   eq_den_dbn_r_connected[ii] = eq_den_dbn_r[ii]-num_den*Dbn_occ
                   eq_den_dbn_r_connected_var = eq_den_dbn_r_var[ii]+(num_den*Dbn_occ_dev)**2+(Dbn_occ*num_den_dev)**2
                   eq_den_dbn_r_connected_std[ii] = np.sqrt(eq_den_dbn_r_connected_var)
 
               filename_eq_den_dbn_r_conn = "%s/Equal_time_Density_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_den_dbn_conn_r,N,U,Mu,dtau,L,Realization)
               filename_eq_den_dbn_r_conn_txt = "%s/Equal_time_Density_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(Text_dir_eq_den_dbn_conn_r,N,U,Mu,dtau,L,Realization)

               Eq_den_dbn_conn_corr_data =  np.stack((eq_den_dbn_r_connected,eq_den_dbn_r_connected_std),axis = 1)
            
               with open(filename_eq_den_dbn_r_conn, 'wb') as outfile:
                  pickle.dump(Eq_den_dbn_conn_corr_data, outfile, pickle.HIGHEST_PROTOCOL)
               
               np.savetxt(filename_eq_den_dbn_r_conn_txt,Eq_den_dbn_conn_corr_data) 
               print("Uff kaka")
            
        run_counter=run_counter+1




def equal_time_moment_moment_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_m2_m2_r,Text_dir_eq_m2_m2_conn_r,N,U,Mu,dtau,L,run_no):


    run_counter = 0
    while(run_counter<int(run_no)):

        Realization = str(run_counter)

        filename_eqm = '%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eqm,N,U,Mu,dtau,L,Realization)
        filename_eq_m2_m2_r = '%s/Equal_time_Moment_Moment_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L,Realization)

        if os.path.exists(filename_eq_m2_m2_r) and os.path.exists(filename_eqm):

            with open(filename_eqm, 'rb') as infile:
                sys_measure = pickle.load(infile)

            with open(filename_eq_m2_m2_r, 'rb') as infile:
               eq_m2_m2_r = pickle.load(infile)
          
            if bool(sys_measure):
               density_st = sys_measure['Average density']
               density = density_st.split(' ')
               num_den = float(density[0].strip(' '))
               num_den_dev = float(density[8].strip(' '))
               nup_ndn_st = sys_measure['Average Nup*Ndn']
               nup_ndn = nup_ndn_st.split(' ')
               Dbn_occ = float(nup_ndn[0].strip(' '))
               Dbn_occ_dev = float(nup_ndn[8].strip(' '))
            #print(Dbn_occ,Dbn_occ_dev,"Doublon_no")
            #print(num_den,num_den_dev,"num_den_meas")


               filename_eq_m2_m2_r_variation = '%s/Equal_time_Moment_Moment_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_m2_m2_r,N,U,Mu,dtau,L,Realization)
               with open(filename_eq_m2_m2_r_variation, 'rb') as infile:
                   eq_m2_m2_r_std = pickle.load(infile)

               eq_m2_m2_r_var = np.power(eq_m2_m2_r_std,2)

               eq_m2_m2_r_connected = np.zeros(len(eq_m2_m2_r))
               eq_m2_m2_r_connected_std = np.zeros(len(eq_m2_m2_r))

               for ii in range(len(eq_m2_m2_r)):
                   eq_m2_m2_r_connected[ii] = eq_m2_m2_r[ii]-(num_den-2*Dbn_occ)**2
                   eq_m2_m2_r_connected_var = eq_m2_m2_r_var[ii]+4*(num_den-2*Dbn_occ)*(num_den-2*Dbn_occ)*(num_den_dev*num_den_dev+4*Dbn_occ_dev*Dbn_occ_dev)
                   eq_m2_m2_r_connected_std[ii] = np.sqrt(eq_m2_m2_r_connected_var)

               filename_eq_m2_m2_r_conn = "%s/Equal_time_Moment_Moment_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_m2_m2_conn_r,N,U,Mu,dtau,L,Realization)
               filename_eq_m2_m2_r_conn_txt = "%s/Equal_time_Moment_Moment_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(Text_dir_eq_m2_m2_conn_r,N,U,Mu,dtau,L,Realization)
            
               Eq_m2_m2_conn_corr_data =  np.stack((eq_m2_m2_r_connected,eq_m2_m2_r_connected_std),axis = 1)

               with open(filename_eq_m2_m2_r_conn, 'wb') as outfile:
                   pickle.dump(Eq_m2_m2_conn_corr_data, outfile, pickle.HIGHEST_PROTOCOL)
                         
               np.savetxt(filename_eq_m2_m2_r_conn_txt,Eq_m2_m2_conn_corr_data)

               print("Uff kaka")
        run_counter=run_counter+1




def equal_time_spin_spin_xx_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_spin_spin_r,Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,run_no):
    
    
    run_counter = 0
    while(run_counter<int(run_no)):

        Realization = str(run_counter)
        
        filename_eqm = '%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eqm,N,U,Mu,dtau,L,Realization)
        filename_eq_sxsx_r = '%s/Equal_time_SxSx_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)

        if os.path.exists(filename_eq_sxsx_r) and os.path.exists(filename_eqm):
                
            with open(filename_eqm, 'rb') as infile:
                sys_measure = pickle.load(infile)
               
            with open(filename_eq_sxsx_r, 'rb') as infile:
               eq_sxsx_r = pickle.load(infile)
            
       
            filename_eq_sxsx_r_variation = '%s/Equal_time_SxSx_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)
            with open(filename_eq_sxsx_r_variation, 'rb') as infile:
               eq_sxsx_r_std = pickle.load(infile)

            if bool(sys_measure):
               density_st = sys_measure['Average density']
               density = density_st.split(' ')
               num_den = float(density[0].strip(' '))
               num_den_dev = float(density[8].strip(' '))

               nup_ndn_st = sys_measure['Average Nup*Ndn']
               nup_ndn = nup_ndn_st.split(' ')
               Dbn_occ = float(nup_ndn[0].strip(' '))
               Dbn_occ_dev = float(nup_ndn[8].strip(' '))

               Up_occ_st = sys_measure['Average up occupancy']
               Dn_occ_st = sys_measure['Average dn occupancy']
               Up_occ = Up_occ_st.split(' ')
               Dn_occ = Dn_occ_st.split(' ')

               Up_den = float(Up_occ[0].strip(' '))
               Dn_den = float(Dn_occ[0].strip(' '))
            
               Up_den_dev = float(Up_occ[8].strip(' '))
               Dn_den_dev = float(Dn_occ[8].strip(' '))
 
            #AF_corr_func_xx_st = sys_measure['AF correlation function (xx)'] 
            #AF_corr_func_zz_st = sys_measure['AF correlation function (zz)'] 
            #Ferro_corr_func_xx_st = sys_measure['Ferro corr. func. (xx)']
            #Ferro_corr_func_zz_st = sys_measure['Ferro corr. func. (zz)'] 

            #print(Dbn_occ,Dbn_occ_dev,"Doublon_no")
            #print(num_den,num_den_dev,"num_den_meas")
            #print(Up_den,Dn_den,Up_den_dev,Dn_den_dev,"up,down occupations")            
                        
               eq_sxsx_r_var = np.power(eq_sxsx_r_std,2)
            
               eq_sxsx_r_connected = np.zeros(len(eq_sxsx_r))
               eq_sxsx_r_connected_std = np.zeros(len(eq_sxsx_r))
                
               for ii in range(len(eq_sxsx_r)):
                   eq_sxsx_r_connected[ii] = eq_sxsx_r[ii]-0.5*(Up_den-Dn_den)
                   eq_sxsx_r_connected_var = eq_sxsx_r_var[ii]+0.25*(Up_den_dev*Up_den_dev+Dn_den_dev*Dn_den_dev)
                   eq_sxsx_r_connected_std[ii] = np.sqrt(eq_sxsx_r_connected_var)
           
               filename_eq_sxsx_r_conn = "%s/Equal_time_SxSx_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,Realization)
               filename_eq_sxsx_r_conn_txt = "%s/Equal_time_SxSx_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,Realization)

               Eq_sxsx_conn_corr_data =  np.stack((eq_sxsx_r_connected,eq_sxsx_r_connected_std),axis = 1)

               with open(filename_eq_sxsx_r_conn, 'wb') as outfile:
                  pickle.dump(Eq_sxsx_conn_corr_data, outfile, pickle.HIGHEST_PROTOCOL)
   
               np.savetxt(filename_eq_sxsx_r_conn_txt,Eq_sxsx_conn_corr_data)

               print("Uff kaka")

        run_counter=run_counter+1


        
def equal_time_spin_spin_zz_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_spin_spin_r,Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,run_no):
    
    
    run_counter = 0
    while(run_counter<int(run_no)):

        Realization = str(run_counter)
        
        filename_eqm = '%s/Thermodynamic_measurements_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eqm,N,U,Mu,dtau,L,Realization)
        filename_eq_szsz_r = '%s/Equal_time_SzSz_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)

        if os.path.exists(filename_eq_szsz_r) and os.path.exists(filename_eqm):
                
            with open(filename_eqm, 'rb') as infile:
                sys_measure = pickle.load(infile)
               
            with open(filename_eq_szsz_r, 'rb') as infile:
               eq_szsz_r = pickle.load(infile)
            
            
            filename_eq_szsz_r_variation = '%s/Equal_time_SzSz_correlation_real_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl' %(Text_dir_eq_spin_spin_r,N,U,Mu,dtau,L,Realization)
            with open(filename_eq_szsz_r_variation, 'rb') as infile:
               eq_szsz_r_std = pickle.load(infile)

            if bool(sys_measure):
               
               density_st = sys_measure['Average density']
               density = density_st.split(' ')
               num_den = float(density[0].strip(' '))
               num_den_dev = float(density[8].strip(' '))

               nup_ndn_st = sys_measure['Average Nup*Ndn']
               nup_ndn = nup_ndn_st.split(' ')
               Dbn_occ = float(nup_ndn[0].strip(' '))
               Dbn_occ_dev = float(nup_ndn[8].strip(' '))

               Up_occ_st = sys_measure['Average up occupancy']
               Dn_occ_st = sys_measure['Average dn occupancy']
               Up_occ = Up_occ_st.split(' ')
               Dn_occ = Dn_occ_st.split(' ')

               Up_den = float(Up_occ[0].strip(' '))
               Dn_den = float(Dn_occ[0].strip(' '))
            
               Up_den_dev = float(Up_occ[8].strip(' '))
               Dn_den_dev = float(Dn_occ[8].strip(' '))
 
            #AF_corr_func_xx_st = sys_measure['AF correlation function (xx)'] 
            #AF_corr_func_zz_st = sys_measure['AF correlation function (zz)'] 
            #Ferro_corr_func_xx_st = sys_measure['Ferro corr. func. (xx)']
            #Ferro_corr_func_zz_st = sys_measure['Ferro corr. func. (zz)'] 

            #print(Dbn_occ,Dbn_occ_dev,"Doublon_no")
            #print(num_den,num_den_dev,"num_den_meas")
            #print(Up_den,Dn_den,Up_den_dev,Dn_den_dev,"up,down occupations")            
                        
               eq_szsz_r_var = np.power(eq_szsz_r_std,2)
            
               eq_szsz_r_connected = np.zeros(len(eq_szsz_r))
               eq_szsz_r_connected_std = np.zeros(len(eq_szsz_r))
                
               for ii in range(len(eq_szsz_r)):
                   eq_szsz_r_connected[ii] = eq_szsz_r[ii]-0.5*(Up_den-Dn_den)
                   eq_szsz_r_connected_var = eq_szsz_r_var[ii]+0.25*(Up_den_dev*Up_den_dev+Dn_den_dev*Dn_den_dev)
                   eq_szsz_r_connected_std[ii] = np.sqrt(eq_szsz_r_connected_var)
           
               filename_eq_szsz_r_conn = "%s/Equal_time_SzSz_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.pkl"%(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,Realization)
               filename_eq_szsz_r_conn_txt = "%s/Equal_time_SzSz_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,Realization)
            
               Eq_szsz_conn_corr_data =  np.stack((eq_szsz_r_connected,eq_szsz_r_connected_std),axis = 1)

               with open(filename_eq_szsz_r_conn, 'wb') as outfile:
                   pickle.dump(Eq_szsz_conn_corr_data, outfile, pickle.HIGHEST_PROTOCOL)
             
               np.savetxt(filename_eq_szsz_r_conn_txt,Eq_szsz_conn_corr_data)

               print("Uff kaka")
            
        run_counter=run_counter+1



        

#=========================================================Calcualting weighted average over connected correlations across realizations====================================================================================================



def equal_time_density_density_connected_correlation_function_real_space_normal_average(Text_dir_eq_den_den_conn_r,N,U,Mu,dtau,L,run_no):


   filename_eq_den_den_conn_r_0 = '%s/Equal_time_Density_Density_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat' %(Text_dir_eq_den_den_conn_r,N,U,Mu,dtau,L)


   if os.path.exists(filename_eq_den_den_conn_r_0):


       eq_den_den_conn_r_0,eq_den_den_conn_r_std_0 = np.loadtxt(filename_eq_den_den_conn_r_0,unpack = 'True',usecols = [0,1])

       r_grid, r_rel_grid_size = r_point_upper_half_grid(N)
       r_sep_2,Eq_den_den_conn_r_nbr,Eq_den_den_conn_r_nbr_std = neighbor_average(r_grid, eq_den_den_conn_r_0, eq_den_den_conn_r_std_0)
       

       run_counter = 0
       while(run_counter<int(run_no)):
            run_counter=run_counter+1
            realization = str(run_counter)


            filename_eq_den_den_conn_r = '%s/Equal_time_Density_Density_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_eq_den_den_conn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_den_den_conn_r):
        
               eq_den_den_conn_r,eq_den_den_conn_r_std = np.loadtxt(filename_eq_den_den_conn_r,unpack = 'True',usecols = [0,1])
           
               r_sep_2, eq_den_den_conn_r_nbr, eq_den_den_conn_r_nbr_std = neighbor_average(r_grid, eq_den_den_conn_r, eq_den_den_conn_r_std)
               
               Eq_X_Den_Den_conn_r_nbr = mat_append(Eq_den_den_conn_r_nbr,eq_den_den_conn_r_nbr)
               Eq_den_den_conn_r_nbr = np.copy(Eq_X_Den_Den_conn_r_nbr)


               Eq_X_Den_Den_conn_r_nbr_std = mat_append(Eq_den_den_conn_r_nbr_std,eq_den_den_conn_r_nbr_std)
               Eq_den_den_conn_r_nbr_std = np.copy(Eq_X_Den_Den_conn_r_nbr_std)


       Eq_den_den_conn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_den_den_conn_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
            Eq_den_den_conn_r_nbr_avg[ii],Eq_den_den_conn_r_nbr_std_avg[ii] = normal_average(Eq_den_den_conn_r_nbr[ii,:],Eq_den_den_conn_r_nbr_std[ii,:])

       
       filename_eq_den_den_conn_r_nbr_avg = '%s/Equal_time_Density_Density_connected_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_den_den_conn_r,N,U,Mu,dtau,L)
       data_eq_den_den_conn_r_nbr_avg = np.stack((r_sep_2,Eq_den_den_conn_r_nbr_avg,Eq_den_den_conn_r_nbr_std_avg),axis = 1)
       np.savetxt(filename_eq_den_den_conn_r_nbr_avg,data_eq_den_den_conn_r_nbr_avg)


       #print("Density up correlation", Eq_dudu_r_avg)
       print("run total",run_counter)
   else:
       print("Error, den den file not found")




def equal_time_doublon_doublon_connected_correlation_function_real_space_normal_average(Text_dir_eq_dbn_dbn_conn_r,N,U,Mu,dtau,L,run_no):



   filename_eq_dbn_dbn_conn_r_0 = '%s/Equal_time_Doublon_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat' %(Text_dir_eq_dbn_dbn_conn_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_dbn_dbn_conn_r_0):

       eq_dbn_dbn_conn_r_0,eq_dbn_dbn_conn_r_std_0 = np.loadtxt(filename_eq_dbn_dbn_conn_r_0,unpack = 'True',usecols = [0,1])

       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_dbn_dbn_conn_r_nbr,Eq_dbn_dbn_conn_r_nbr_std = neighbor_average(r_grid, eq_dbn_dbn_conn_r_0, eq_dbn_dbn_conn_r_std_0)       

       run_counter = 0
       while(run_counter<int(run_no)):
            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_dbn_dbn_conn_r = '%s/Equal_time_Doublon_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_eq_dbn_dbn_conn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_dbn_dbn_conn_r):


               eq_dbn_dbn_conn_r,eq_dbn_dbn_conn_r_std = np.loadtxt(filename_eq_dbn_dbn_conn_r,unpack = 'True',usecols = [0,1])
               
               r_sep_2, eq_dbn_dbn_conn_r_nbr, eq_dbn_dbn_conn_r_nbr_std = neighbor_average(r_grid, eq_dbn_dbn_conn_r, eq_dbn_dbn_conn_r_std)

               Eq_X_Dbn_Dbn_conn_r_nbr = mat_append(Eq_dbn_dbn_conn_r_nbr,eq_dbn_dbn_conn_r_nbr)
               Eq_dbn_dbn_conn_r_nbr = np.copy(Eq_X_Dbn_Dbn_conn_r_nbr)

               Eq_X_Dbn_Dbn_conn_r_nbr_std = mat_append(Eq_dbn_dbn_conn_r_nbr_std,eq_dbn_dbn_conn_r_nbr_std)
               Eq_dbn_dbn_conn_r_nbr_std = np.copy(Eq_X_Dbn_Dbn_conn_r_nbr_std)


       Eq_dbn_dbn_conn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_dbn_dbn_conn_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_dbn_dbn_conn_r_nbr_avg[ii],Eq_dbn_dbn_conn_r_nbr_std_avg[ii] = normal_average(Eq_dbn_dbn_conn_r_nbr[ii,:],Eq_dbn_dbn_conn_r_nbr_std[ii,:])


       filename_eq_dbn_dbn_conn_r_nbr_avg = '%s/Equal_time_Doublon_Doublon_connected_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_dbn_dbn_conn_r,N,U,Mu,dtau,L)
       data_eq_dbn_dbn_conn_r_nbr_avg = np.stack((r_sep_2,Eq_dbn_dbn_conn_r_nbr_avg,Eq_dbn_dbn_conn_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_dbn_dbn_conn_r_nbr_avg,data_eq_dbn_dbn_conn_r_nbr_avg)
       print("run total",run_counter)
   else:
       print("Error, dbn dbn file not found")


       
def equal_time_moment_moment_connected_correlation_function_real_space_normal_average(Text_dir_eq_m2_m2_conn_r,N,U,Mu,dtau,L,run_no):


   filename_eq_m2_m2_conn_r_0 = '%s/Equal_time_Moment_Moment_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat' %(Text_dir_eq_m2_m2_conn_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_m2_m2_conn_r_0):
       
       eq_m2_m2_conn_r_0,eq_m2_m2_conn_r_std_0 = np.loadtxt(filename_eq_m2_m2_conn_r_0,unpack = 'True', usecols = [0,1])

       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_m2_m2_conn_r_nbr,Eq_m2_m2_conn_r_nbr_std = neighbor_average(r_grid, eq_m2_m2_conn_r_0, eq_m2_m2_conn_r_std_0)
         

       run_counter = 0
       while(run_counter<int(run_no)):
            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_m2_m2_conn_r = '%s/Equal_time_Moment_Moment_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_eq_m2_m2_conn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_m2_m2_conn_r):


               eq_m2_m2_conn_r,eq_m2_m2_conn_r_std = np.loadtxt(filename_eq_m2_m2_conn_r,unpack = 'True',usecols = [0,1])  
               r_sep_2,eq_m2_m2_conn_r_nbr,eq_m2_m2_conn_r_nbr_std = neighbor_average(r_grid,eq_m2_m2_conn_r,eq_m2_m2_conn_r_std)
               
               Eq_X_m2_m2_conn_r_nbr = mat_append(Eq_m2_m2_conn_r_nbr,eq_m2_m2_conn_r_nbr)
               Eq_m2_m2_conn_r_nbr = np.copy(Eq_X_m2_m2_conn_r_nbr)

               Eq_X_m2_m2_conn_r_nbr_std = mat_append(Eq_m2_m2_conn_r_nbr_std,eq_m2_m2_conn_r_nbr_std)
               Eq_m2_m2_conn_r_nbr_std = np.copy(Eq_X_m2_m2_conn_r_nbr_std)


       Eq_m2_m2_conn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_m2_m2_conn_r_nbr_std_avg = np.zeros(len(r_sep_2))

       for ii in range(len(r_sep_2)):
           Eq_m2_m2_conn_r_nbr_avg[ii],Eq_m2_m2_conn_r_nbr_std_avg[ii] = normal_average(Eq_m2_m2_conn_r_nbr[ii,:],Eq_m2_m2_conn_r_nbr_std[ii,:])

       filename_eq_m2_m2_conn_r_nbr_avg = '%s/Equal_time_Moment_Moment_connected_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_m2_m2_conn_r,N,U,Mu,dtau,L)
       data_eq_m2_m2_conn_r_nbr_avg = np.stack((r_sep_2,Eq_m2_m2_conn_r_nbr_avg,Eq_m2_m2_conn_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_m2_m2_conn_r_nbr_avg,data_eq_m2_m2_conn_r_nbr_avg)
       print("run total",run_counter)
   else:
       print("Error, m2 m2 file not found")
       


def equal_time_density_doublon_connected_correlation_function_real_space_normal_average(Text_dir_eq_den_dbn_conn_r,N,U,Mu,dtau,L,run_no):



   filename_eq_den_dbn_conn_r_0 = '%s/Equal_time_Density_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat' %(Text_dir_eq_den_dbn_conn_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_den_dbn_conn_r_0):

       eq_den_dbn_conn_r_0,eq_den_dbn_conn_r_std_0 = np.loadtxt(filename_eq_den_dbn_conn_r_0,unpack = 'True',usecols = [0,1])

       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_den_dbn_conn_r_nbr,Eq_den_dbn_conn_r_nbr_std = neighbor_average(r_grid, eq_den_dbn_conn_r_0, eq_den_dbn_conn_r_std_0)       

       run_counter = 0
       while(run_counter<int(run_no)):
            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_den_dbn_conn_r = '%s/Equal_time_Density_Doublon_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_eq_den_dbn_conn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_den_dbn_conn_r):


               eq_den_dbn_conn_r,eq_den_dbn_conn_r_std = np.loadtxt(filename_eq_den_dbn_conn_r,unpack = 'True',usecols = [0,1])
               
               r_sep_2, eq_den_dbn_conn_r_nbr, eq_den_dbn_conn_r_nbr_std = neighbor_average(r_grid, eq_den_dbn_conn_r, eq_den_dbn_conn_r_std)

               Eq_X_Den_Dbn_conn_r_nbr = mat_append(Eq_den_dbn_conn_r_nbr,eq_den_dbn_conn_r_nbr)
               Eq_den_dbn_conn_r_nbr = np.copy(Eq_X_Den_Dbn_conn_r_nbr)

               Eq_X_Den_Dbn_conn_r_nbr_std = mat_append(Eq_den_dbn_conn_r_nbr_std,eq_den_dbn_conn_r_nbr_std)
               Eq_den_dbn_conn_r_nbr_std = np.copy(Eq_X_Den_Dbn_conn_r_nbr_std)


       Eq_den_dbn_conn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_den_dbn_conn_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_den_dbn_conn_r_nbr_avg[ii],Eq_den_dbn_conn_r_nbr_std_avg[ii] = normal_average(Eq_den_dbn_conn_r_nbr[ii,:],Eq_den_dbn_conn_r_nbr_std[ii,:])


       filename_eq_den_dbn_conn_r_nbr_avg = '%s/Equal_time_Density_Doublon_connected_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_den_dbn_conn_r,N,U,Mu,dtau,L)
       data_eq_den_dbn_conn_r_nbr_avg = np.stack((r_sep_2,Eq_den_dbn_conn_r_nbr_avg,Eq_den_dbn_conn_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_den_dbn_conn_r_nbr_avg,data_eq_den_dbn_conn_r_nbr_avg)
       print("run total",run_counter)
   else:
       print("Error, den dbn file not found")


def equal_time_spin_spin_xx_connected_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,run_no):


   filename_eq_sxsx_conn_r_0 = '%s/Equal_time_SxSx_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat' %(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_sxsx_conn_r_0):
       print("file lol")
       eq_sxsx_conn_r_0,eq_sxsx_conn_r_std_0 = np.loadtxt(filename_eq_sxsx_conn_r_0,unpack = 'True',usecols = [0,1])

       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_sxsx_conn_r_nbr,Eq_sxsx_conn_r_nbr_std = neighbor_average(r_grid, eq_sxsx_conn_r_0, eq_sxsx_conn_r_std_0)       

       run_counter = 0
       while(run_counter<int(run_no)):
            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_sxsx_conn_r = '%s/Equal_time_SxSx_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_sxsx_conn_r):


               eq_sxsx_conn_r,eq_sxsx_conn_r_std = np.loadtxt(filename_eq_sxsx_conn_r,unpack = 'True',usecols = [0,1])
               
               r_sep_2, eq_sxsx_conn_r_nbr, eq_sxsx_conn_r_nbr_std = neighbor_average(r_grid, eq_sxsx_conn_r, eq_sxsx_conn_r_std)

               Eq_X_SxSx_conn_r_nbr = mat_append(Eq_sxsx_conn_r_nbr,eq_sxsx_conn_r_nbr)
               Eq_sxsx_conn_r_nbr = np.copy(Eq_X_SxSx_conn_r_nbr)

               Eq_X_SxSx_conn_r_nbr_std = mat_append(Eq_sxsx_conn_r_nbr_std,eq_sxsx_conn_r_nbr_std)
               Eq_sxsx_conn_r_nbr_std = np.copy(Eq_X_SxSx_conn_r_nbr_std)


       Eq_sxsx_conn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_sxsx_conn_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_sxsx_conn_r_nbr_avg[ii],Eq_sxsx_conn_r_nbr_std_avg[ii] = normal_average(Eq_sxsx_conn_r_nbr[ii,:],Eq_sxsx_conn_r_nbr_std[ii,:])


       filename_eq_sxsx_conn_r_nbr_avg = '%s/Equal_time_SxSx_connected_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L)
       data_eq_sxsx_conn_r_nbr_avg = np.stack((r_sep_2,Eq_sxsx_conn_r_nbr_avg,Eq_sxsx_conn_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_sxsx_conn_r_nbr_avg,data_eq_sxsx_conn_r_nbr_avg)
       print("run total",run_counter)
   else:
       print("Error, sx sx file not found")


       
def equal_time_spin_spin_zz_connected_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,run_no):


   filename_eq_szsz_conn_r_0 = '%s/Equal_time_SzSz_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat' %(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L)

   if os.path.exists(filename_eq_szsz_conn_r_0):
       print("file lol")
       eq_szsz_conn_r_0,eq_szsz_conn_r_std_0 = np.loadtxt(filename_eq_szsz_conn_r_0,unpack = 'True',usecols = [0,1])

       r_grid, r_rel_grid_size = r_point_full_grid(N)

       r_sep_2,Eq_szsz_conn_r_nbr,Eq_szsz_conn_r_nbr_std = neighbor_average(r_grid, eq_szsz_conn_r_0, eq_szsz_conn_r_std_0)       

       run_counter = 0
       while(run_counter<int(run_no)):
            run_counter=run_counter+1
            realization = str(run_counter)

            filename_eq_szsz_conn_r = '%s/Equal_time_SzSz_connected_correlation_real_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat' %(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L,realization)

            if os.path.exists(filename_eq_szsz_conn_r):


               eq_szsz_conn_r,eq_szsz_conn_r_std = np.loadtxt(filename_eq_szsz_conn_r,unpack = 'True',usecols = [0,1])
               
               r_sep_2, eq_szsz_conn_r_nbr, eq_szsz_conn_r_nbr_std = neighbor_average(r_grid, eq_szsz_conn_r, eq_szsz_conn_r_std)

               Eq_X_SzSz_conn_r_nbr = mat_append(Eq_szsz_conn_r_nbr,eq_szsz_conn_r_nbr)
               Eq_szsz_conn_r_nbr = np.copy(Eq_X_SzSz_conn_r_nbr)

               Eq_X_SzSz_conn_r_nbr_std = mat_append(Eq_szsz_conn_r_nbr_std,eq_szsz_conn_r_nbr_std)
               Eq_szsz_conn_r_nbr_std = np.copy(Eq_X_SzSz_conn_r_nbr_std)


       Eq_szsz_conn_r_nbr_avg = np.zeros(len(r_sep_2))
       Eq_szsz_conn_r_nbr_std_avg = np.zeros(len(r_sep_2))


       for ii in range(len(r_sep_2)):
           Eq_szsz_conn_r_nbr_avg[ii],Eq_szsz_conn_r_nbr_std_avg[ii] = normal_average(Eq_szsz_conn_r_nbr[ii,:],Eq_szsz_conn_r_nbr_std[ii,:])


       filename_eq_szsz_conn_r_nbr_avg = '%s/Equal_time_SzSz_connected_correlation_real_space_neighbor_averaged_normal_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat' %(Text_dir_eq_spin_spin_conn_r,N,U,Mu,dtau,L)
       data_eq_szsz_conn_r_nbr_avg = np.stack((r_sep_2,Eq_szsz_conn_r_nbr_avg,Eq_szsz_conn_r_nbr_std_avg),axis =1)
       np.savetxt(filename_eq_szsz_conn_r_nbr_avg,data_eq_szsz_conn_r_nbr_avg)
       print("run total",run_counter)
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


#==================================================================Directory for saving neighbor averaged connected correlators for each run========================================================================================

    
    Text_dir_eq_den_den_conn_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_density_connected_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_den_den_conn_r):
       os.makedirs(Text_dir_eq_den_den_conn_r)

    Text_dir_eq_dbn_dbn_conn_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Doublon_doublon_connected_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_dbn_dbn_conn_r):
       os.makedirs(Text_dir_eq_dbn_dbn_conn_r)

    Text_dir_eq_den_dbn_conn_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_doublon_connected_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_den_dbn_conn_r):
       os.makedirs(Text_dir_eq_den_dbn_conn_r)
       
    Text_dir_eq_m2_m2_conn_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Moment_moment_connected_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_m2_m2_conn_r):
       os.makedirs(Text_dir_eq_m2_m2_conn_r)
       
    Text_dir_eq_spin_spin_conn_r = '/fs/byo/trivedi/roy.369/Research/Determinant_Quantum_Monte_Carlo/UC_code/legacy-qmc/Cold_atom_correlators/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_spin_connected_correlation_functions'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_eq_spin_spin_conn_r):
       os.makedirs(Text_dir_eq_spin_spin_conn_r)

#=============================================================================Redundant directories, may be useful for something else================================================================================================       
       
    #Text_dir_eq_den_den_r_avg = '/Users/sayantanroy/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_density_correlation_functions_neighbor_averaged_weighted_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    #if not os.path.exists(Text_dir_eq_den_den_r_avg):
    #   os.makedirs(Text_dir_eq_den_den_r_avg)

    #Text_dir_eq_dbn_dbn_r_avg = '/Users/sayantanroy/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Doublon_doublon_correlation_functions_neighbor_averaged_weighted_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    #if not os.path.exists(Text_dir_eq_dbn_dbn_r_avg):
    #   os.makedirs(Text_dir_eq_dbn_dbn_r_avg)

    #Text_dir_eq_den_dbn_r_avg = '/Users/sayantanroy/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Density_doublon_correlation_functions_neighbor_averaged_weighted_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    #if not os.path.exists(Text_dir_eq_den_dbn_r_avg):
    #   os.makedirs(Text_dir_eq_den_dbn_r_avg)

    #Text_dir_eq_m2_m2_r_avg = '/Users/sayantanroy/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Moment_moment_correlation_functions_neighbor_averaged_weighted_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    #if not os.path.exists(Text_dir_eq_m2_m2_r_avg):
    #   os.makedirs(Text_dir_eq_m2_m2_r_avg)

    #Text_dir_eq_spin_spin_r_avg = '/Users/sayantanroy/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Text_files/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spin_spin_correlation_functions_neighbor_averaged_weighted_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    #if not os.path.exists(Text_dir_eq_spin_spin_r_avg):
    #   os.makedirs(Text_dir_eq_spin_spin_r_avg)

#==========================================================================================

    equal_time_density_density_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_den_den_r,Text_dir_eq_den_den_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_doublon_doublon_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_dbn_dbn_r,Text_dir_eq_dbn_dbn_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_moment_moment_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_m2_m2_r,Text_dir_eq_m2_m2_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_density_doublon_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_den_dbn_r,Text_dir_eq_den_dbn_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_spin_spin_xx_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_spin_spin_r,Text_dir_eq_spin_spin_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_spin_spin_zz_connected_correlation_function_real_space(Text_dir_eqm,Text_dir_eq_spin_spin_r,Text_dir_eq_spin_spin_conn_r,N,U,mu,Dtau,L,N_runs)

    print("avg_running")
    equal_time_density_density_connected_correlation_function_real_space_normal_average(Text_dir_eq_den_den_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_doublon_doublon_connected_correlation_function_real_space_normal_average(Text_dir_eq_dbn_dbn_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_moment_moment_connected_correlation_function_real_space_normal_average(Text_dir_eq_m2_m2_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_density_doublon_connected_correlation_function_real_space_normal_average(Text_dir_eq_den_dbn_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_spin_spin_xx_connected_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_conn_r,N,U,mu,Dtau,L,N_runs)
    equal_time_spin_spin_zz_connected_correlation_function_real_space_normal_average(Text_dir_eq_spin_spin_conn_r,N,U,mu,Dtau,L,N_runs)
    print("done,", U,mu,L)

if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)
