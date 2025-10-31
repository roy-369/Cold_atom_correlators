#/usr/bin/env python3
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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colorbar
from matplotlib import rc
import matplotlib as mpl
mpl.rcParams['axes.linewidth'] = 3
import matplotlib as mpl
mpl.rc('font',family='Times New Roman')
from matplotlib.pyplot import cm
from scipy import integrate



def arr_loc(val,X):
   ind = 0
   for i in range(len(X)):
      if(round(float(X[i]),2) == round(val,2)):
         ind = i
         break
   return ind


def normal_average(data_set,data_set_var):

   data_avg = 0
   data_avg_var = 0
   
   a = 0
   
   for i in range(len(data_set)):
       data_avg = data_avg+data_set[i]
       data_avg_var = data_avg_var+data_set_var[i]
       a= a+1
   

   return data_avg/a, np.sqrt(data_avg_var/(a*a))
   
def numerical_derivative_calc(x_array,y_array,y_array_std):

    dy_dx = np.zeros(len(x_array))
    dy_dx_std = np.zeros(len(x_array))

    for i in range(len(x_array)):

        if i == 0:
           dx = x_array[i+1]-x_array[i]+10**(-6)
           dy_dx[i] = (y_array[i+1]-y_array[i])/dx
           dy_dx_std[i] = np.sqrt((y_array_std[i+1]*y_array_std[i+1]+y_array_std[i]*y_array_std[i])/(dx*dx))
  
        if i > 0 and i < len(x_array)-1:
           dx = x_array[i+1]-x_array[i-1]+10**(-6)
           dy_dx[i] = (y_array[i+1]-y_array[i-1])/dx
           dy_dx_std[i] = np.sqrt((y_array_std[i+1]*y_array_std[i+1]+y_array_std[i-1]*y_array_std[i-1])/(dx*dx))

        if i == len(x_array)-1:
           dx = x_array[i]-x_array[i-1]+10**(-6)
           dy_dx[i] = (y_array[i]-y_array[i-1])/dx
           dy_dx_std[i] = np.sqrt((y_array_std[i]*y_array_std[i]+y_array_std[i-1]*y_array_std[i-1])/(dx*dx))

    return dy_dx,dy_dx_std

def numerical_integrate_trapz(x_array,y_array,y_array_std):
    Intg_trapz = 0
    Intg_trapz_var = 0

    if len(x_array) == 2:
       dx = x_array[1]-x_array[0]
       Intg_trapz = Intg_trapz+0.5*(y_array[0]+y_array[1])*dx
       Intg_trapz_var = Intg_trapz_var+0.25*(y_array_std[1]*y_array_std[1]+y_array_std[0]*y_array_std[0])*dx*dx

    if len(x_array) >2:
       for i in range(1,len(x_array)):

           dx = x_array[i]-x_array[i-1]
           Intg_trapz = Intg_trapz+0.5*(y_array[i-1]+y_array[i])*dx
           Intg_trapz_var = Intg_trapz_var+0.25*(y_array_std[i]*y_array_std[i]+y_array_std[i-1]*y_array_std[i-1])*dx*dx

    return Intg_trapz,np.sqrt(Intg_trapz_var)


     
def Bz_label_coordinate_parser(K,k_class):
   k_label = K[:,0]
   kx = K[:,1]
   ky = K[:,2]
   b = 0
   kx_cord = kx[0]
   ky_cord = ky[0]
#   print(kx[12],ky[12],"Max_point")
   for i in range(len(k_label)):
       if b == 1:
          break
       if(k_label[i] == k_class+1):
          kx_cord = kx[i]
          ky_cord = ky[i]
          b = 1
   return kx_cord,ky_cord



def plot_func_multiple_beta(N,dtau,trot,plot_title,plot_x_label,plot_y_label,x_val,y_val,y_val_var,plot_name,plot_legend_name,x_cord,x_min,x_max):

   T_val = np.zeros(len(trot))
   
   color_1 = iter(cm.gnuplot(np.linspace(0, 1, len(trot)+1)))
   plt.figure(figsize = (25,20))
   plt.xticks(x_cord,fontsize = 80)
   plt.yticks(fontsize = 80)
   #ax = plt.gca()
   
   for jj in range(len(trot)-1,-1,-1):
       c_1 = next(color_1)         
       beta = (float(dtau))*(float(trot[jj]))
       T_val[jj] = round(1/beta,2)
       plt.errorbar(x_val[:,jj],y_val[:,jj],yerr = y_val_var[:,jj],c=c_1,marker = 'o',markersize = 20, linewidth = 5, elinewidth = 3, capsize = 5,label = r"T = %s"%str(T_val[jj]))
   plt.grid()
   plt.xlim(x_min,x_max)
   #plt.axvline(x = 1,color="black", linestyle=":")
   plt.axhline(y = 0,color="black", linestyle=":")
   plt.savefig(plot_name)
   plt.close()

   color_1 = iter(cm.gnuplot(np.linspace(0, 1, len(trot)+1)))

   fig = plt.figure(figsize = (30,30))
   ax = plt.subplot(111)
   for jj in range(len(trot)-1,-1,-1):
       c_1 = next(color_1)         
       ax.errorbar(x_val[:,jj],y_val[:,jj],yerr = y_val_var[:,jj],c=c_1,marker = 'o',markersize = 20, linewidth = 5, elinewidth = 3, capsize = 5,label = r"T = %s"%str(T_val[jj]))
   box = ax.get_position()
   ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
   ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=False, ncol=1,fontsize=40)
   #plt.legend(loc = 'upper right', fontsize = 40)
   plt.savefig(plot_legend_name)
   plt.close() 


def plot_func_multiple_u(N,u,plot_title,plot_x_label,plot_y_label,x_val,y_val,y_val_var,plot_name,plot_legend_name,x_cord,x_min,x_max):

   color_2 = iter(cm.tab10(np.linspace(0, 1, len(u)+1)))

   plt.figure(figsize = (25,20))
   plt.xticks(x_cord,fontsize = 80)
   plt.yticks(fontsize = 80)
   for jj in range(len(u)):
       c_2 = next(color_2)
       plt.errorbar(x_val[:,jj],y_val[:,jj],yerr = y_val_var[:,jj],c=c_2,marker = 'o',markersize = 15, linewidth = 5,elinewidth = 3, capsize = 5,label = r"U=%s"%str(u[jj]))

   plt.grid(True)
   plt.xlim(x_min,x_max)
   #plt.axvline(x = 1,color="black", linestyle=":")
   plt.axhline(y = 0,color="black", linestyle=":")
   #plt.legend(fontsize = 60, loc = 'best')
   plt.savefig(plot_name)
   plt.close()

   color_2 = iter(cm.tab10(np.linspace(0, 1, len(u)+1)))

   fig = plt.figure(figsize = (35,20))
   ax = plt.subplot(111)
   for jj in range(len(u)):
       c_2 = next(color_2)         
       ax.errorbar(x_val[:,jj],y_val[:,jj],yerr = y_val_var[:,jj],c=c_2,marker = 'o',markersize = 20, linewidth = 5, elinewidth = 3, capsize = 5,label = r"U = %s"%str(u[jj]))
   box = ax.get_position()
   ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
   ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=False, ncol=2,fontsize=40)
   #plt.legend(loc = 'upper right', fontsize = 40)
   plt.savefig(plot_legend_name)
   plt.close()


   
def plot_func_multiple_beta_half_filling(Graph_dir_hf,N,dtau,trot,plot_title,plot_x_label,plot_y_label,x_val,y_val,y_val_var,plot_name,plot_legend_name,x_cord,x_min,x_max):

   color_3 = iter(cm.gnuplot(np.linspace(0, 1, len(trot)+1)))
   plt.figure(figsize = (25,20))
   plt.xticks(x_cord,fontsize = 80)
   plt.yticks(fontsize = 80)
   for jj in range(len(trot)-1,-1,-1):
       c_3 = next(color_3)
       beta = (float(dtau))*(float(trot[jj]))
       T = round(1/beta,2)
       plt.errorbar(x_val[:,jj],y_val[:,jj],yerr = y_val_var[:,jj],c=c_3,marker = 'o',markersize = 20, linewidth = 5, elinewidth = 3, capsize = 5,label = r"T = %s"%str(T))
   plt.grid()
   plt.xlim(x_min,x_max)
   if plot_name == "%s/Magnetization_N_%s_dtau_%s_half_filling.png"%(Graph_dir_hf,N,dtau):
      plt.ylim(0.5,1.0)
   #plt.axvline(x = 1,color="black", linestyle=":")
   plt.axhline(y = 0,color="black", linestyle=":")
   if plot_name == "%s/Magnetization_N_%s_dtau_%s_half_filling.png"%(Graph_dir_hf,N,dtau):
       plt.yscale("log")
   #plt.legend(fontsize = 60, loc = 'best')
   plt.savefig(plot_name)
   plt.close()

   color_3 = iter(cm.gnuplot(np.linspace(0, 1, len(trot)+1)))

   fig = plt.figure(figsize = (30,20))
   ax = plt.subplot(111)
   for jj in range(len(trot)-1,-1,-1):
       c_3 = next(color_3)
       beta = (float(dtau))*(float(trot[jj]))
       T = round(1/beta,2)
       ax.errorbar(x_val[:,jj],y_val[:,jj],yerr = y_val_var[:,jj],c=c_3,marker = 'o',markersize = 20, linewidth = 5, elinewidth = 3, capsize = 5,label = r"T = %s"%str(T))
   box = ax.get_position()
   #ax.set_position([box.x0, box.y0 + box.height * 0.1, box.width, box.height * 0.9])
   ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05),fancybox=True, shadow=False, ncol=len(trot),fontsize=40)
   box = ax.get_position()
   ax.set_position([box.x0, box.y0, box.width * 0.8, box.height])

# Put a legend to the right of the current axis
   ax.legend(loc='center left', bbox_to_anchor=(1, 0.5),fancybox=True, shadow=False, ncol=1,fontsize=40)
   #plt.legend(loc = 'upper right', fontsize = 40)
   plt.savefig(plot_legend_name)
   plt.close()

   
   
def brillouin_zone_labels(Text_dir,Graph_dir,N):
   filename_k = '%s/Brillouin_zone_co-oridinates_N_%s.pkl' %(Text_dir,N)
   with open(filename_k, 'rb') as infile:
        K = pickle.load(infile)
  

#   print(K,"K_grid")
#   print(K.shape)
   k_label = K[:,0]
   kx = K[:,1]
   ky = K[:,2]

   cm = plt.cm.get_cmap('viridis')
   plt.figure(figsize = (20,20))
   plt.title(r"Brillouin_zone, $N = %sx%s$"%(N,N),fontsize = 40)
   ax = plt.gca()
   plt.xlabel('Kx',fontsize = 40)
   plt.ylabel('Ky',fontsize = 40)
   scat = plt.scatter(kx,ky,c=k_label,s =1000,cmap=cm)
   plt.xticks(fontsize = 40)
   plt.yticks(fontsize = 40)
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
  # plt.colorbar(scat, cax=cax)
   cbar = plt.colorbar(scat, cax= cax)
   cbar.ax.tick_params(labelsize = 30)
   plt.legend(loc = 'upper right')
   for i, txt in enumerate(k_label):
      ax.annotate(txt, (kx[i], ky[i]))
   plt.savefig('%s/Brillouin_zone_%d_sites.png'%(Graph_dir,int(N)))
   plt.close()
   return K


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
       print("no realization found at U, mu ,L", u,mu,L)


def thermodynamic_entropy_seebeck_calc(Text_dir_main,N,U,Mu,dtau,Trot):

   Density = np.zeros((len(U),len(Mu),len(Trot)))
   dN_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
   dN_dT = np.zeros((len(U),len(Mu),len(Trot)))

   Density_dev = np.zeros((len(U),len(Mu),len(Trot)))
   dN_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
   dN_dT_std = np.zeros((len(U),len(Mu),len(Trot)))

   Entropy = np.zeros((len(U),len(Mu),len(Trot)))
   S_k = np.zeros((len(U),len(Mu),len(Trot)))
   Entropy_std = np.zeros((len(U),len(Mu),len(Trot)))
   S_k_std = np.zeros((len(U),len(Mu),len(Trot)))    
   Entropy_hf = np.zeros((len(U),len(Trot)))
   
   for i in range(len(U)):
       for j in range(len(Mu)):
           for k in range(len(Trot)):
               Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],dtau,Mu[j],dtau,Trot[k])
               filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
               with open(filename_eqm_avg, 'rb') as infile:
                    sys_measure_avg = pickle.load(infile)

               filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
               with open(filename_eqm_std, 'rb') as infile:
                    sys_measure_std = pickle.load(infile)

               Density[i][j][k] = sys_measure_avg['Density averaged']
               Density_dev[i][j][k] = sys_measure_std['Density standard deviation']
       
   Mu_val = np.zeros(len(Mu))
   for j in range(len(Mu)):
       Mu_val[j] = float(Mu[j])
    
   Beta = np.zeros(len(Trot))
   for k in range(len(Trot)):
       Beta[k] = float(dtau)*float(Trot[k])
   
   U_val = np.zeros(len(U))
   for i in range(len(U)):
       U_val[i] = float(U[i])
    
   for i in range(len(U)):
       for j in range(len(Mu)):
           
           dN_dBeta[i,j,:],dN_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,Density[i,j,:],Density_dev[i,j,:])
           
   for k in range(len(Trot)):
       beta = float(dtau)*float(Trot[k])
       T = 1/beta
       dN_dT[:,:,k] = (-1/(T*T))*dN_dBeta[:,:,k]
       dN_dT_std[:,:,k] = (1/(T*T))*dN_dBeta_std[:,:,k]

   
   for i in range(len(U)):
       for j in range(len(Mu)):
           for k in range(len(Trot)):
               Entropy[i,j,k],Entropy_std[i,j,k] = numerical_integrate_trapz(Mu_val[:j+1],dN_dT[i,:j+1,k],dN_dT_std[i,:j+1,k])
   
   for i in range(len(U)):
       for k in range(len(Trot)):
           S_k[i,:,k],S_k_std[i,:,k] = numerical_derivative_calc(Density[i,:,k],Entropy[i,:,k],Entropy_std[i,:,k])

   for i in range(len(U)):
       
       Text_dir_thermo = "%s/Text_files_N_%s_thermodynamic_entropy_seebeck/Text_files_N_%s_U_%s_dtau_%s"%(Text_dir_main,N,N,U[i],dtau)
       if not os.path.exists(Text_dir_thermo):
          os.makedirs(Text_dir_thermo)

       for k in range(len(Trot)):

           Entropy_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],Entropy[i,:,k],Entropy_std[i,:,k]),axis = 1)
           S_kelvin_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],S_k[i,:,k],S_k_std[i,:,k]),axis = 1)

           filename_entropy = "%s/Thermodynamic_entropy_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
           filename_seebeck = "%s/Seebeck_kelvin_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])

           np.savetxt(filename_entropy,Entropy_data)
           np.savetxt(filename_seebeck,S_kelvin_data) 

#=========================================================================================================================================


def doublon_derivative_calc(Text_dir_main,N,U,Mu,dtau,Trot):
 
    Density = np.zeros((len(U),len(Mu),len(Trot)))
    Doublon = np.zeros((len(U),len(Mu),len(Trot)))
    Doublon_corr = np.zeros((len(U),len(Mu),len(Trot)))
    dD_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
    dD_corr_dBeta = np.zeros((len(U),len(Mu),len(Trot)))

    dD_dT = np.zeros((len(U),len(Mu),len(Trot)))
    dD_corr_dT = np.zeros((len(U),len(Mu),len(Trot)))
    Doublon_hf = np.zeros((len(U),len(Trot)))
    dD_dBeta_hf = np.zeros((len(U),len(Trot)))
    dD_dT_hf = np.zeros((len(U),len(Trot)))

    Density_dev = np.zeros((len(U),len(Mu),len(Trot)))
    Doublon_dev = np.zeros((len(U),len(Mu),len(Trot)))
    dD_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
    dD_dT_std = np.zeros((len(U),len(Mu),len(Trot)))
    dD_corr_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
    dD_corr_dT_std = np.zeros((len(U),len(Mu),len(Trot)))
    Doublon_hf_std = np.zeros((len(U),len(Trot)))
    dD_dBeta_hf_std = np.zeros((len(U),len(Trot)))
    dD_dT_hf_std = np.zeros((len(U),len(Trot)))

 
    for i in range(len(U)):
        for j in range(len(Mu)):
            for k in range(len(Trot)):
                Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],dtau,Mu[j],dtau,Trot[k])
                filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_avg, 'rb') as infile:
                     sys_measure_avg = pickle.load(infile)
 
                filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_std, 'rb') as infile:
                     sys_measure_std = pickle.load(infile)
 
                Density[i][j][k] = sys_measure_avg['Density averaged']
                Density_dev[i][j][k] = sys_measure_std['Density standard deviation']
                Doublon[i][j][k] = sys_measure_avg['Doublon number averaged']
                Doublon_dev[i][j][k] = sys_measure_std['Doublon number standard deviation']
                Doublon_corr[i][j][k] = Doublon[i][j][k]-0.25*Density[i][j][k]*Density[i][j][k]

    Mu_val = np.zeros(len(Mu))
    for j in range(len(Mu)):
        Mu_val[j] = float(Mu[j])
 
    Beta = np.zeros(len(Trot))
    T_val = np.zeros(len(Trot))

    for k in range(len(Trot)):
        Beta[k] = float(dtau)*float(Trot[k])
 
    U_val = np.zeros(len(U))
    for i in range(len(U)):
        U_val[i] = float(U[i])

    zero_ind = arr_loc(0.00,Mu_val)
    
    for i in range(len(U)):
        for j in range(len(Mu)):
 
            dD_dBeta[i,j,:],dD_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,Doublon[i,j,:],Doublon_dev[i,j,:])
            dD_corr_dBeta[i,j,:],dD_corr_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,Doublon_corr[i,j,:],Doublon_dev[i,j,:])

        dD_dBeta_hf[i,:],dD_dBeta_hf_std[i,:] = numerical_derivative_calc(Beta,Doublon[i,zero_ind,:],Doublon_dev[i,zero_ind,:])


    for k in range(len(Trot)):
        beta = float(dtau)*float(Trot[k])
        T = 1/beta
        T_val[k] = T
        dD_dT[:,:,k] = (-1/(T*T))*dD_dBeta[:,:,k]
        dD_dT_std[:,:,k] = (1/(T*T))*dD_dBeta_std[:,:,k]
        dD_corr_dT[:,:,k] = (-1/(T*T))*dD_corr_dBeta[:,:,k]
        dD_corr_dT_std[:,:,k] = (1/(T*T))*dD_corr_dBeta_std[:,:,k]

        dD_dT_hf[:,k] = (-1/(T*T))*dD_dBeta_hf[:,k]
        dD_dT_hf_std[:,k] = (1/(T*T))*dD_dBeta_hf_std[:,k] 
        Doublon_hf[:,k] = np.copy(Doublon[:,zero_ind,k])
        Doublon_hf_std[:,k] = np.copy(Doublon_dev[:,zero_ind,k])
 
    Text_dir_thermo_hf = "%s/Text_files_N_%s_doublon_number_half_filling/Text_files_N_%s_dtau_%s"%(Text_dir_main,N,N,dtau)
    if not os.path.exists(Text_dir_thermo_hf):
       os.makedirs(Text_dir_thermo_hf)
 
    for i in range(len(U)):
 
        Text_dir_thermo = "%s/Text_files_N_%s_doublon_number/Text_files_N_%s_U_%s_dtau_%s"%(Text_dir_main,N,N,U[i],dtau)
        if not os.path.exists(Text_dir_thermo):
           os.makedirs(Text_dir_thermo) 

        for k in range(len(Trot)):
 
            Doublon_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],Doublon[i,:,k],Doublon_dev[i,:,k],dD_dT[i,:,k],dD_dT_std[i,:,k],dD_corr_dT[i,:,k],dD_corr_dT_std[i,:,k]),axis = 1)
            filename_doublon = "%s/Doublon_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
            np.savetxt(filename_doublon,Doublon_data)

        Doublon_hf_data = np.stack((T_val,Doublon_hf[i,:],Doublon_hf_std[i,:],dD_dT_hf[i,:],dD_dT_hf_std[i,:]),axis = 1)
        filename_doublon_hf = "%s/Doublon_data_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
        np.savetxt(filename_doublon_hf,Doublon_hf_data)
 
  #=========================================================================================================================================

def kinetic_energy_derivative_calc(Text_dir_main,N,U,Mu,dtau,Trot):

    Density = np.zeros((len(U),len(Mu),len(Trot)))
    KE = np.zeros((len(U),len(Mu),len(Trot)))
    dKE_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
    dKE_dT = np.zeros((len(U),len(Mu),len(Trot)))

    Density_dev = np.zeros((len(U),len(Mu),len(Trot)))
    KE_dev = np.zeros((len(U),len(Mu),len(Trot)))
    dKE_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
    dKE_dT_std = np.zeros((len(U),len(Mu),len(Trot)))


    for i in range(len(U)):
        for j in range(len(Mu)):
            for k in range(len(Trot)):
                Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],dtau,Mu[j],dtau,Trot[k])
                filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_avg, 'rb') as infile:
                     sys_measure_avg = pickle.load(infile)

                filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_std, 'rb') as infile:
                     sys_measure_std = pickle.load(infile)

                Density[i][j][k] = sys_measure_avg['Density averaged']
                Density_dev[i][j][k] = sys_measure_std['Density standard deviation']

                KE[i][j][k] = sys_measure_avg['Kinetic energy averaged']
                KE_dev[i][j][k] = sys_measure_std['Kinetic energy standard deviation']


    Mu_val = np.zeros(len(Mu))
    for j in range(len(Mu)):
        Mu_val[j] = float(Mu[j])

    Beta = np.zeros(len(Trot))
    T_val = np.zeros(len(Trot))

    for k in range(len(Trot)):
        Beta[k] = float(dtau)*float(Trot[k])

    U_val = np.zeros(len(U))
    for i in range(len(U)):
        U_val[i] = float(U[i])

    zero_ind = 1 #arr_loc(0.00,Mu_val)

    for i in range(len(U)):
        for j in range(len(Mu)):

            dKE_dBeta[i,j,:],dKE_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,KE[i,j,:],KE_dev[i,j,:])

    for k in range(len(Trot)):
        beta = float(dtau)*float(Trot[k])
        T = 1/beta
        T_val[k] = T
        dKE_dT[:,:,k] = (-1/(T*T))*dKE_dBeta[:,:,k]
        dKE_dT_std[:,:,k] = (1/(T*T))*dKE_dBeta_std[:,:,k]


    for i in range(len(U)):

        Text_dir_thermo = "%s/Text_files_N_%s_kinetic_energy/Text_files_N_%s_U_%s_dtau_%s"%(Text_dir_main,N,N,U[i],dtau)
        if not os.path.exists(Text_dir_thermo):
           os.makedirs(Text_dir_thermo)

        for k in range(len(Trot)):

            KE_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],KE[i,:,k],KE_dev[i,:,k],dKE_dT[i,:,k],dKE_dT_std[i,:,k]),axis = 1)
            filename_KE = "%s/Kinetic_energy_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
            np.savetxt(filename_KE,KE_data)


  #=========================================================================================================================================



 
def local_moment_derivative_calc(Text_dir_main,N,U,Mu,dtau,Trot):

    Density = np.zeros((len(U),len(Mu),len(Trot)))
    Moment = np.zeros((len(U),len(Mu),len(Trot)))
    dM_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
    dM_dT = np.zeros((len(U),len(Mu),len(Trot)))
    Moment_hf = np.zeros((len(U),len(Trot)))
    dM_dBeta_hf = np.zeros((len(U),len(Trot)))
    dM_dT_hf = np.zeros((len(U),len(Trot)))

    Density_dev = np.zeros((len(U),len(Mu),len(Trot)))
    Moment_dev = np.zeros((len(U),len(Mu),len(Trot)))
    dM_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
    dM_dT_std = np.zeros((len(U),len(Mu),len(Trot)))
    Moment_hf_std = np.zeros((len(U),len(Trot)))
    dM_dBeta_hf_std = np.zeros((len(U),len(Trot)))
    dM_dT_hf_std = np.zeros((len(U),len(Trot)))


    for i in range(len(U)):
        for j in range(len(Mu)):
            for k in range(len(Trot)):
                Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],dtau,Mu[j],dtau,Trot[k])
                filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_avg, 'rb') as infile:
                     sys_measure_avg = pickle.load(infile)

                filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_std, 'rb') as infile:
                     sys_measure_std = pickle.load(infile)

                Density[i][j][k] = sys_measure_avg['Density averaged']
                Density_dev[i][j][k] = sys_measure_std['Density standard deviation']
                doublon = sys_measure_avg['Doublon number averaged']
                doublon_dev = sys_measure_std['Doublon number standard deviation']
                Moment[i][j][k] = Density[i][j][k]-2*doublon
                Moment_dev[i][j][k] = np.sqrt(Density_dev[i][j][k]*Density_dev[i][j][k]+4*doublon_dev*doublon_dev)


    Mu_val = np.zeros(len(Mu))
    for j in range(len(Mu)):
        Mu_val[j] = float(Mu[j])

    Beta = np.zeros(len(Trot))
    T_val = np.zeros(len(Trot))

    for k in range(len(Trot)):
        Beta[k] = float(dtau)*float(Trot[k])

    U_val = np.zeros(len(U))
    for i in range(len(U)):
        U_val[i] = float(U[i])

    zero_ind = arr_loc(0.00,Mu_val)

    for i in range(len(U)):
        for j in range(len(Mu)):

            dM_dBeta[i,j,:],dM_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,Moment[i,j,:],Moment_dev[i,j,:])
                                                                                                                                                                           
        dM_dBeta_hf[i,:],dM_dBeta_hf_std[i,:] = numerical_derivative_calc(Beta,Moment[i,zero_ind,:],Moment_dev[i,zero_ind,:])


    for k in range(len(Trot)):
        beta = float(dtau)*float(Trot[k])
        T = 1/beta
        T_val[k] = T
        dM_dT[:,:,k] = (-1/(T*T))*dM_dBeta[:,:,k]
        dM_dT_std[:,:,k] = (1/(T*T))*dM_dBeta_std[:,:,k]
        dM_dT_hf[:,k] = (-1/(T*T))*dM_dBeta_hf[:,k]
        dM_dT_hf_std[:,k] = (1/(T*T))*dM_dBeta_hf_std[:,k]
        Moment_hf[:,k] = np.copy(Moment[:,zero_ind,k])
        Moment_hf_std[:,k] = np.copy(Moment_dev[:,zero_ind,k])

    Text_dir_thermo_hf = "%s/Text_files_N_%s_local_moment_half_filling/Text_files_N_%s_dtau_%s"%(Text_dir_main,N,N,dtau)
    if not os.path.exists(Text_dir_thermo_hf):
       os.makedirs(Text_dir_thermo_hf)

    for i in range(len(U)):

        Text_dir_thermo = "%s/Text_files_N_%s_local_moment/Text_files_N_%s_U_%s_dtau_%s"%(Text_dir_main,N,N,U[i],dtau)
        if not os.path.exists(Text_dir_thermo):
           os.makedirs(Text_dir_thermo)

        for k in range(len(Trot)):

            Moment_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],Moment[i,:,k],Moment_dev[i,:,k],dM_dT[i,:,k],dM_dT_std[i,:,k]),axis = 1)
            filename_moment = "%s/Moment_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
            np.savetxt(filename_moment,Moment_data)

        Moment_hf_data = np.stack((T_val,Moment_hf[i,:],Moment_hf_std[i,:],dM_dT_hf[i,:],dM_dT_hf_std[i,:]),axis = 1)
        filename_moment_hf = "%s/Moment_data_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
        np.savetxt(filename_moment_hf,Moment_hf_data)

  #=========================================================================================================================================


def energy_derivative_calc(Text_dir_main,N,U,Mu,dtau,Trot):

    Density = np.zeros((len(U),len(Mu),len(Trot)))
    KE = np.zeros((len(U),len(Mu),len(Trot)))
    PE = np.zeros((len(U),len(Mu),len(Trot)))
    En = np.zeros((len(U),len(Mu),len(Trot)))


    Density_dev = np.zeros((len(U),len(Mu),len(Trot)))
    KE_dev = np.zeros((len(U),len(Mu),len(Trot)))
    PE_dev = np.zeros((len(U),len(Mu),len(Trot)))
    En_dev = np.zeros((len(U),len(Mu),len(Trot)))


    dKE_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
    dKE_dT = np.zeros((len(U),len(Mu),len(Trot)))
    KE_hf = np.zeros((len(U),len(Trot)))
    dKE_dBeta_hf = np.zeros((len(U),len(Trot)))
    dKE_dT_hf = np.zeros((len(U),len(Trot)))

    dPE_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
    dPE_dT = np.zeros((len(U),len(Mu),len(Trot)))
    PE_hf = np.zeros((len(U),len(Trot)))
    dPE_dBeta_hf = np.zeros((len(U),len(Trot)))
    dPE_dT_hf = np.zeros((len(U),len(Trot)))

    dEn_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
    dEn_dT = np.zeros((len(U),len(Mu),len(Trot)))
    En_hf = np.zeros((len(U),len(Trot)))
    dEn_dBeta_hf = np.zeros((len(U),len(Trot)))
    dEn_dT_hf = np.zeros((len(U),len(Trot)))

    dKE_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
    dKE_dT_std = np.zeros((len(U),len(Mu),len(Trot)))
    KE_hf_std = np.zeros((len(U),len(Trot)))
    dKE_dBeta_hf_std = np.zeros((len(U),len(Trot)))
    dKE_dT_hf_std = np.zeros((len(U),len(Trot)))

    dPE_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
    dPE_dT_std = np.zeros((len(U),len(Mu),len(Trot)))
    PE_hf_std = np.zeros((len(U),len(Trot)))
    dPE_dBeta_hf_std = np.zeros((len(U),len(Trot)))
    dPE_dT_hf_std = np.zeros((len(U),len(Trot)))

    dEn_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
    dEn_dT_std = np.zeros((len(U),len(Mu),len(Trot)))
    En_hf_std = np.zeros((len(U),len(Trot)))
    dEn_dBeta_hf_std = np.zeros((len(U),len(Trot)))
    dEn_dT_hf_std = np.zeros((len(U),len(Trot)))


    for i in range(len(U)):
        for j in range(len(Mu)):
            for k in range(len(Trot)):
                Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],dtau,Mu[j],dtau,Trot[k])
                filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_avg, 'rb') as infile:
                     sys_measure_avg = pickle.load(infile)

                filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
                with open(filename_eqm_std, 'rb') as infile:
                     sys_measure_std = pickle.load(infile)

                Density[i][j][k] = sys_measure_avg['Density averaged']
                Density_dev[i][j][k] = sys_measure_std['Density standard deviation']
                En[i][j][k] = sys_measure_avg['Total energy averaged']-float(Mu[j])*sys_measure_avg['Density averaged'] 
                KE[i][j][k] = sys_measure_avg['Kinetic energy averaged']-float(Mu[j])*sys_measure_avg['Density averaged'] 
                PE[i][j][k] = float(U[i])*sys_measure_avg['Doublon number averaged']

                En_dev[i][j][k] = np.sqrt(sys_measure_std['Total energy standard deviation']*sys_measure_std['Total energy standard deviation']+float(Mu[j])*float(Mu[j])*sys_measure_std['Density standard deviation']*sys_measure_std['Density standard deviation']) 
                KE_dev[i][j][k] = np.sqrt(sys_measure_std['Kinetic energy standard deviation']*sys_measure_std['Kinetic energy standard deviation']+float(Mu[j])*float(Mu[j])*sys_measure_std['Density standard deviation']*sys_measure_std['Density standard deviation']) 
                PE_dev[i][j][k] = float(U[i])*sys_measure_std['Doublon number standard deviation']



    Mu_val = np.zeros(len(Mu))
    for j in range(len(Mu)):
        Mu_val[j] = float(Mu[j])

    Beta = np.zeros(len(Trot))
    T_val = np.zeros(len(Trot))

    for k in range(len(Trot)):
        Beta[k] = float(dtau)*float(Trot[k])

    U_val = np.zeros(len(U))
    for i in range(len(U)):
        U_val[i] = float(U[i])

    zero_ind = arr_loc(0.00,Mu_val) 


    for i in range(len(U)):
        for j in range(len(Mu)):

            dKE_dBeta[i,j,:],dKE_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,KE[i,j,:],KE_dev[i,j,:])
            dPE_dBeta[i,j,:],dPE_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,PE[i,j,:],PE_dev[i,j,:])
            dEn_dBeta[i,j,:],dEn_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,En[i,j,:],En_dev[i,j,:])

        dKE_dBeta_hf[i,:],dKE_dBeta_hf_std[i,:] = numerical_derivative_calc(Beta,KE[i,zero_ind,:],KE_dev[i,zero_ind,:])
        dPE_dBeta_hf[i,:],dPE_dBeta_hf_std[i,:] = numerical_derivative_calc(Beta,PE[i,zero_ind,:],PE_dev[i,zero_ind,:])
        dEn_dBeta_hf[i,:],dEn_dBeta_hf_std[i,:] = numerical_derivative_calc(Beta,En[i,zero_ind,:],En_dev[i,zero_ind,:])

    for k in range(len(Trot)):
        beta = float(dtau)*float(Trot[k])
        T = 1/beta
        T_val[k] = T

        dKE_dT[:,:,k] = (-1/(T*T))*dKE_dBeta[:,:,k]
        dKE_dT_std[:,:,k] = (1/(T*T))*dKE_dBeta_std[:,:,k]
        dKE_dT_hf[:,k] = (-1/(T*T))*dKE_dBeta_hf[:,k]
        dKE_dT_hf_std[:,k] = (1/(T*T))*dKE_dBeta_hf_std[:,k]
        KE_hf[:,k] = np.copy(KE[:,zero_ind,k])
        KE_hf_std[:,k] = np.copy(KE_dev[:,zero_ind,k])

        dPE_dT[:,:,k] = (-1/(T*T))*dPE_dBeta[:,:,k]
        dPE_dT_std[:,:,k] = (1/(T*T))*dPE_dBeta_std[:,:,k]
        dPE_dT_hf[:,k] = (-1/(T*T))*dPE_dBeta_hf[:,k]
        dPE_dT_hf_std[:,k] = (1/(T*T))*dPE_dBeta_hf_std[:,k]
        PE_hf[:,k] = np.copy(PE[:,zero_ind,k])
        PE_hf_std[:,k] = np.copy(PE_dev[:,zero_ind,k])

        dEn_dT[:,:,k] = (-1/(T*T))*dEn_dBeta[:,:,k]
        dEn_dT_std[:,:,k] = (1/(T*T))*dEn_dBeta_std[:,:,k]
        dEn_dT_hf[:,k] = (-1/(T*T))*dEn_dBeta_hf[:,k]
        dEn_dT_hf_std[:,k] = (1/(T*T))*dEn_dBeta_hf_std[:,k]
        En_hf[:,k] = np.copy(En[:,zero_ind,k])
        En_hf_std[:,k] = np.copy(En_dev[:,zero_ind,k])


    Text_dir_thermo_hf = "%s/Text_files_N_%s_energies_half_filling/Text_files_N_%s_dtau_%s"%(Text_dir_main,N,N,dtau)
    if not os.path.exists(Text_dir_thermo_hf):
       os.makedirs(Text_dir_thermo_hf)

    for i in range(len(U)):

        Text_dir_thermo = "%s/Text_files_N_%s_energies/Text_files_N_%s_U_%s_dtau_%s"%(Text_dir_main,N,N,U[i],dtau)
        if not os.path.exists(Text_dir_thermo):
           os.makedirs(Text_dir_thermo)

        for k in range(len(Trot)):

            KE_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],KE[i,:,k],KE_dev[i,:,k],dKE_dT[i,:,k],dKE_dT_std[i,:,k]),axis = 1)
            filename_KE = "%s/Kinetic_Energy_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
            np.savetxt(filename_KE,KE_data)

            PE_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],PE[i,:,k],PE_dev[i,:,k],dPE_dT[i,:,k],dPE_dT_std[i,:,k]),axis = 1)
            filename_PE = "%s/Potential_Energy_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
            np.savetxt(filename_PE,PE_data)

            En_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],En[i,:,k],En_dev[i,:,k],dEn_dT[i,:,k],dEn_dT_std[i,:,k]),axis = 1)
            filename_En = "%s/Total_Energy_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
            np.savetxt(filename_En,En_data)


        KE_hf_data = np.stack((T_val,KE_hf[i,:],KE_hf_std[i,:],dKE_dT_hf[i,:],dKE_dT_hf_std[i,:]),axis = 1)
        filename_KE_hf = "%s/Kinetic_Energy_data_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
        np.savetxt(filename_KE_hf,KE_hf_data)

        PE_hf_data = np.stack((T_val,PE_hf[i,:],PE_hf_std[i,:],dPE_dT_hf[i,:],dPE_dT_hf_std[i,:]),axis = 1)
        filename_PE_hf = "%s/Potential_Energy_data_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
        np.savetxt(filename_PE_hf,PE_hf_data)

        En_hf_data = np.stack((T_val,En_hf[i,:],En_hf_std[i,:],dEn_dT_hf[i,:],dEn_dT_hf_std[i,:]),axis = 1)
        filename_En_hf = "%s/Total_Energy_data_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
        np.savetxt(filename_En_hf,En_hf_data)
  #=========================================================================================================================================





def compressibility_charge_gap_calc(Text_dir_main,N,U,Mu,dtau,Trot):

   Density = np.zeros((len(U),len(Mu),len(Trot)))
   dN_dmu = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dBeta = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dT = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dBeta_hf = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf = np.zeros((len(U),len(Trot)))

   Density_dev = np.zeros((len(U),len(Mu),len(Trot)))
   dN_dmu_std = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dBeta_std = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dT_std = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dBeta_hf_std = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf_std = np.zeros((len(U),len(Trot)))

   dN_dmu_num_drv = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dBeta_hf_num_drv = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf_num_drv = np.zeros((len(U),len(Trot)))

   dN_dmu_std_num_drv = np.zeros((len(U),len(Mu),len(Trot)))
   dKappa_dBeta_hf_std_num_drv = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf_std_num_drv = np.zeros((len(U),len(Trot)))

   for i in range(len(U)):
       for j in range(len(Mu)):
           for k in range(len(Trot)):
               Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],dtau,Mu[j],dtau,Trot[k])
               filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
               with open(filename_eqm_avg, 'rb') as infile:
                    sys_measure_avg = pickle.load(infile)

               filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
               with open(filename_eqm_std, 'rb') as infile:
                    sys_measure_std = pickle.load(infile)

               Density[i][j][k] = sys_measure_avg['Density averaged']
               Density_dev[i][j][k] = sys_measure_std['Density standard deviation']


   Mu_val = np.zeros(len(Mu))
   for j in range(len(Mu)):
       Mu_val[j] = float(Mu[j])

   zero_ind = arr_loc(0.00,Mu_val)

   Beta = np.zeros(len(Trot))
   for k in range(len(Trot)):
       Beta[k] = float(dtau)*float(Trot[k])

   U_val = np.zeros(len(U))
   for i in range(len(U)):
       U_val[i] = float(U[i])

   for i in range(len(U)):
       for k in range(len(Trot)):

           dN_dmu[i,:,k],dN_dmu_std[i,:,k] = numerical_derivative_calc(Mu_val,Density[i,:,k],Density_dev[i,:,k])
           dN_dmu_num_drv[i,:,k] = np.gradient(Density[i,:,k],Mu_val)
           dN_dmu_std_num_drv[i,:,k] = np.copy(Density_dev[i,:,k])

       dKappa_dBeta_hf[i,:],dKappa_dBeta_hf_std[i,:] = numerical_derivative_calc(Beta,dN_dmu[i,zero_ind,:],dN_dmu_std[i,zero_ind,:])
       dKappa_dBeta_hf_num_drv[i,:] = np.gradient(dN_dmu_num_drv[i,zero_ind,:],Beta)
       dKappa_dBeta_hf_std_num_drv[i,:] = np.copy(Density_dev[i,zero_ind,:])


   for i in range(len(U)):
       for j in range(len(Mu)):

           dKappa_dBeta[i,j,:],dKappa_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,dN_dmu[i,j,:],dN_dmu_std[i,j,:])

   T_val = np.zeros(len(Beta))
   for k in range(len(Trot)):
       beta = float(dtau)*float(Trot[k])
       T = 1/beta
       T_val[k] = T
       dKappa_dT_hf[:,k] = (-1/(T*T))*dKappa_dBeta_hf[:,k]
       dKappa_dT_hf_std[:,k] = (1/(T*T))*dKappa_dBeta_hf_std[:,k]
       dKappa_dT_hf_num_drv[:,k] = (-1/(T*T))*dKappa_dBeta_hf_num_drv[:,k]
       dKappa_dT_hf_std_num_drv[:,k] = (1/(T*T))*dKappa_dBeta_hf_std_num_drv[:,k]
       dKappa_dT[:,:,k] = (-1/(T*T))*dKappa_dBeta[:,:,k] 
       dKappa_dT_std[:,:,k] = (1/(T*T))*dKappa_dBeta_std[:,:,k]


   Text_dir_thermo_hf = "%s/Text_files_N_%s_charge_gap_half_filling/Text_files_N_%s_dtau_%s"%(Text_dir_main,N,N,dtau)
   if not os.path.exists(Text_dir_thermo_hf):
      os.makedirs(Text_dir_thermo_hf)

   for i in range(len(U)):
       
       Text_dir_thermo = "%s/Text_files_N_%s_charge_gap/Text_files_N_%s_U_%s_dtau_%s"%(Text_dir_main,N,N,U[i],dtau)
       if not os.path.exists(Text_dir_thermo):
          os.makedirs(Text_dir_thermo)

       for k in range(len(Trot)):

           Kappa_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],dN_dmu[i,:,k],dN_dmu_std[i,:,k]),axis = 1)
           filename_kappa = "%s/Charge_gap_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
           np.savetxt(filename_kappa,Kappa_data)

           dKappa_dT_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],dKappa_dT[i,:,k],dKappa_dT_std[i,:,k]),axis = 1)
           filename_dkappa_dT = "%s/dkappa_dT_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
           np.savetxt(filename_dkappa_dT,dKappa_dT_data)


       Kappa_hf_data = np.stack((T_val,dN_dmu[i,zero_ind,:],dN_dmu_std[i,zero_ind,:],dKappa_dT_hf[i,:],dKappa_dT_hf_std[i,:]),axis = 1)
       filename_kappa_hf = "%s/Charge_gap_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
       np.savetxt(filename_kappa_hf,Kappa_hf_data)

       Kappa_hf_data_num_drv = np.stack((T_val,dKappa_dT_hf_num_drv[i,:],dKappa_dT_hf_std_num_drv[i,:]),axis = 1)
       filename_kappa_hf_num_drv = "%s/Charge_gap_numpy_gradient_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
       np.savetxt(filename_kappa_hf_num_drv,Kappa_hf_data_num_drv)

#=========================================================================================================================================


def compressibility_charge_gap_calc_half_filling(Text_dir_main,N,U,Mu,dtau,Trot):

   Density = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dN_dmu = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dBeta = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dT = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dBeta_hf = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf = np.zeros((len(U),len(Trot)))

   Density_dev = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dN_dmu_std = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dBeta_std = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dT_std = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dBeta_hf_std = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf_std = np.zeros((len(U),len(Trot)))

   dN_dmu_num_drv = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dBeta_hf_num_drv = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf_num_drv = np.zeros((len(U),len(Trot)))

   dN_dmu_std_num_drv = np.zeros((len(U),len(Mu)+1,len(Trot)))
   dKappa_dBeta_hf_std_num_drv = np.zeros((len(U),len(Trot)))
   dKappa_dT_hf_std_num_drv = np.zeros((len(U),len(Trot)))


   for i in range(len(U)):
       for k in range(len(Trot)):
           for j in range(len(Mu)):
               Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],dtau,Mu[j],dtau,Trot[k])
               filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
               with open(filename_eqm_avg, 'rb') as infile:
                    sys_measure_avg = pickle.load(infile)

               filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U[i],Mu[j],dtau,Trot[k])
               with open(filename_eqm_std, 'rb') as infile:
                    sys_measure_std = pickle.load(infile)

               Density[i][j][k] = sys_measure_avg['Density averaged']
               Density_dev[i][j][k] = sys_measure_std['Density standard deviation']
           Density[i][-1][k] = 2-Density[i][0][k]
           Density_dev[i][-1][k] = Density_dev[i][0][k]


   Mu_val = np.zeros(len(Mu)+1)
   for j in range(len(Mu)):
       Mu_val[j] = float(Mu[j])
   Mu_val[-1] = -1*Mu_val[0]

   zero_ind = arr_loc(0.00,Mu_val)

   Beta = np.zeros(len(Trot))
   for k in range(len(Trot)):
       Beta[k] = float(dtau)*float(Trot[k])

   U_val = np.zeros(len(U))
   for i in range(len(U)):
       U_val[i] = float(U[i])

   for i in range(len(U)):
       for k in range(len(Trot)):

           dN_dmu[i,:,k],dN_dmu_std[i,:,k] = numerical_derivative_calc(Mu_val,Density[i,:,k],Density_dev[i,:,k])
           dN_dmu_num_drv[i,:,k] = np.gradient(Density[i,:,k],Mu_val)
           dN_dmu_std_num_drv[i,:,k] = np.copy(Density_dev[i,:,k])

       dKappa_dBeta_hf[i,:],dKappa_dBeta_hf_std[i,:] = numerical_derivative_calc(Beta,dN_dmu[i,zero_ind,:],dN_dmu_std[i,zero_ind,:])
       dKappa_dBeta_hf_num_drv[i,:] = np.gradient(dN_dmu_num_drv[i,zero_ind,:],Beta)
       dKappa_dBeta_hf_std_num_drv[i,:] = np.copy(Density_dev[i,zero_ind,:])


   for i in range(len(U)):
       for j in range(len(Mu)+1):

           dKappa_dBeta[i,j,:],dKappa_dBeta_std[i,j,:] = numerical_derivative_calc(Beta,dN_dmu[i,j,:],dN_dmu_std[i,j,:])

   T_val = np.zeros(len(Beta))
   for k in range(len(Trot)):
       beta = float(dtau)*float(Trot[k])
       T = 1/beta
       T_val[k] = T
       dKappa_dT_hf[:,k] = (-1/(T*T))*dKappa_dBeta_hf[:,k]
       dKappa_dT_hf_std[:,k] = (1/(T*T))*dKappa_dBeta_hf_std[:,k]
       dKappa_dT_hf_num_drv[:,k] = (-1/(T*T))*dKappa_dBeta_hf_num_drv[:,k]
       dKappa_dT_hf_std_num_drv[:,k] = (1/(T*T))*dKappa_dBeta_hf_std_num_drv[:,k]
       dKappa_dT[:,:,k] = (-1/(T*T))*dKappa_dBeta[:,:,k]
       dKappa_dT_std[:,:,k] = (1/(T*T))*dKappa_dBeta_std[:,:,k]


   Text_dir_thermo_hf = "%s/Text_files_N_%s_charge_gap_half_filling/Text_files_N_%s_dtau_%s"%(Text_dir_main,N,N,dtau)
   if not os.path.exists(Text_dir_thermo_hf):
      os.makedirs(Text_dir_thermo_hf)


   for i in range(len(U)):

      # Text_dir_thermo = "%s/Text_files_N_%s_charge_gap/Text_files_N_%s_U_%s_dtau_%s"%(Text_dir_main,N,N,U[i],dtau)
      # if not os.path.exists(Text_dir_thermo):
      #    os.makedirs(Text_dir_thermo)

      # for k in range(len(Trot)):

      #     Kappa_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],dN_dmu[i,:,k],dN_dmu_std[i,:,k]),axis = 1)
      #     filename_kappa = "%s/Charge_gap_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
      #     np.savetxt(filename_kappa,Kappa_data)

      #     dKappa_dT_data = np.stack((Mu_val,Density[i,:,k],Density_dev[i,:,k],dKappa_dT[i,:,k],dKappa_dT_std[i,:,k]),axis = 1)
      #     filename_dkappa_dT = "%s/dkappa_dT_data_N_%s_U_%s_dtau_%s_L_%s.dat"%(Text_dir_thermo,N,U[i],dtau,Trot[k])
      #     np.savetxt(filename_dkappa_dT,dKappa_dT_data)


       Kappa_hf_data = np.stack((T_val,dN_dmu[i,zero_ind,:],dN_dmu_std[i,zero_ind,:],dKappa_dT_hf[i,:],dKappa_dT_hf_std[i,:]),axis = 1)
       filename_kappa_hf = "%s/Charge_gap_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
       np.savetxt(filename_kappa_hf,Kappa_hf_data)

       Kappa_hf_data_num_drv = np.stack((T_val,dKappa_dT_hf_num_drv[i,:],dKappa_dT_hf_std_num_drv[i,:]),axis = 1)
       filename_kappa_hf_num_drv = "%s/Charge_gap_numpy_gradient_half_filling_data_N_%s_U_%s_dtau_%s.dat"%(Text_dir_thermo_hf,N,U[i],dtau)
       np.savetxt(filename_kappa_hf_num_drv,Kappa_hf_data_num_drv)




def main(total,cmdargs):
    if(total!=4):
        raise ValueError('missing args')

    N = cmdargs[1]
    Dtau = cmdargs[2]
    Run_no=int(cmdargs[3])
    U = ["3.10","3.20","3.30","3.40","3.50","3.60","3.70","3.80","3.90","4.00","4.10","4.20","4.30","4.40","4.50","5.00"] #"3.00","3.50","4.00","4.50","5.00","5.50","6.00"] #["1.0","2.0","3.0","5.0","7.0","9.0"] #"3.00","3.50","4.00","4.50","5.00","5.50","6.00","6.50","7.00","7.50","8.00","8.50","9.00","9.50","10.00"] #,"5.0"] #,"5.50"] #["1.0","2.0","3.0","3.50","4.00","4.50","5.00","5.50","6.0","7.0","8.00","9.0"] #["3.20","3.40","3.50","3.60","3.70","3.80","3.90","4.00","4.10","4.20","4.30","4.40","4.50","4.60","4.70","4.80","4.90","5.00"] #,"5.50","6.00","6.50","7.00","7.50","8.00","8.50","9.00","9.50","10.00"] 
    #Mu = ["-0.20","0.00","0.20","0.40","0.60","0.80","1.00","1.20","1.40","1.60","1.80","2.00","2.20"]
    #Mu = ["-10.00","-9.50","-9.00","-8.50","-8.00","-7.50","-7.00","-6.50","-6.00","-5.50","-5.00","-4.50","-4.00","-3.50","-3.00","-2.50","-2.00","-1.50","-1.00","-0.50","0.00"] #,"0.50","1.00","1.50","2.00","2.50","3.00","3.50","4.00","4.50","5.00","5.50","6.00","6.50","7.00","7.50","8.00","8.50","9.00","9.50","10.00"]
    Mu = ["-0.20","0.00"] #,"0.20","0.40","0.60","0.80","1.00","1.20","1.40","1.60","1.80","2.00","2.20"] #,"2.40","2.60","2.80","3.00","3.20","3.40","3.60","3.80","4.00","4.50","5.00","5.50","6.00","6.50","7.00","7.50","8.00","8.50","9.00","9.50","10.00"]
    Trot = ["18","20","22","24","26","28","30","32","34","36","38","40","42"] #,"42","44","46","48","50"] #,"42","44","46","48","50","52","54","60"] #,"24","26","28","30","32","34","36","38","40","42","44","46","48","50","52","54","56","58","60","62"]

    Text_dir_main = "../../Text_files"
 


   
#=====================================================Averaging over multiple realizations =================================

    for i in range(len(U)):
        for j in range(len(Mu)):
            for k in range(len(Trot)):
                Text_dir_eqm = '%s/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements'%(Text_dir_main,N,N,U[i],Dtau,Mu[j],Dtau,Trot[k])
                thermodynamic_measurements_average(Text_dir_eqm,N,U[i],Mu[j],Dtau,Trot[k],Run_no)

    #thermodynamic_entropy_seebeck_calc(Text_dir_main,N,U,Mu,Dtau,Trot)
    compressibility_charge_gap_calc_half_filling(Text_dir_main,N,U,Mu,Dtau,Trot)
    #doublon_derivative_calc(Text_dir_main,N,U,Mu,Dtau,Trot)
    #local_moment_derivative_calc(Text_dir_main,N,U,Mu,Dtau,Trot)
   # energy_derivative_calc(Text_dir_main,N,U,Mu,Dtau,Trot)

if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)












