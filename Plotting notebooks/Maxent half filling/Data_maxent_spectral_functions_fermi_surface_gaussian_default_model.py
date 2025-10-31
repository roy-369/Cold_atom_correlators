#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 30 10:11:44 2022

@author: roy.369
"""


import numpy as np
import pickle
import matplotlib as mpl
import matplotlib.pyplot as plt
import os
import sys
from scipy.interpolate import griddata
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colorbar
from matplotlib import rc
from scipy.optimize import fsolve
from scipy.integrate import simpson
from numpy import trapz
from numpy.polynomial.polynomial import polyfit
from scipy import integrate


def fermi_distribution(en,mu,beta):
    fd = 0
    if beta>500:
       if en<=mu:
          fd = 1
       if en>mu:
          fd  = 0

    else:
       fd = 1.0/(np.exp(beta*(en-mu))+1.0)
    
    return fd

def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx,array[idx]

def radius_theta_cons(kx_mesh,ky_mesh):
   
  
   r = np.zeros(kx_mesh.shape)
   theta = np.zeros(ky_mesh.shape)
   for i in range(len(r[:,0])):
       for j in range(len(r[0,:])):
           r[i][j] = np.abs(np.sqrt((kx_mesh[i][j])**2+(ky_mesh[i][j])**2))
           if(r[i][j]>0):
              if(kx_mesh[i][j] == 0):
                 if(ky_mesh[i][j]>0):
                    theta[i][j] = np.arctan(np.inf)
                 if(ky_mesh[i][j]<0):
                    theta[i][j] = np.arctan(-1*np.inf)
              if(kx_mesh[i][j]>0):
                    theta[i][j] = np.arctan(ky_mesh[i][j]/kx_mesh[i][j])
              if(kx_mesh[i][j]<0):
                    theta[i][j] = np.pi+np.arctan(ky_mesh[i][j]/kx_mesh[i][j])

   return r,theta
   
   
def number_density_calc(chemical_potential,beta):

   initial_x = 0
   initial_y = 0
   final_x = np.pi
   final_y = np.pi

   dx = np.pi/50
   dy = np.pi/50

   x_grid_cord = np.arange(-np.pi,np.pi,dx)
   y_grid_cord = np.arange(-np.pi,np.pi,dy)

   Ek = 0 #np.zeros((len(x_grid_cord),len(y_grid_cord)))
   nk_noninteracting = 0 #np.zeros((len(x_grid_cord),len(y_grid_cord)))

   t = 1
   bz_area = 0

   for i in range(len(x_grid_cord)):
       for j in range(len(y_grid_cord)):
           Ek = -2*t*(np.cos(x_grid_cord[i])+np.cos(y_grid_cord[j]))
           nk_noninteracting = nk_noninteracting+2*fermi_distribution(Ek,chemical_potential,beta)
           bz_area = bz_area+1

   Nk = nk_noninteracting/bz_area
   return Nk
   
   
def brillouin_zone_labels(Text_dir,N):
   filename_k = '%s/Brillouin_zone_co-oridinates_N_%s.pkl' %(Text_dir,N)
   with open(filename_k, 'rb') as infile:
        K = pickle.load(infile)
  
   kx = K[:,0]
   ky = K[:,1]
   K_grid = np.zeros((3,len(kx)))
   
   for j in range(len(kx)):
       K_grid[0][j] = j
       K_grid[1][j] = kx[j]
       K_grid[2][j] = ky[j]
   print(K_grid,"BZ")
   return K_grid


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


def FS_locator(klabel,kxgrid,kygrid,Kx_FS_cord,Ky_FS_cord):
 
    k_label_FS = 0   
    a = 0
    for i in range(len(kxgrid)):
        if a == 1:
           break
        if kxgrid[i] == Kx_FS_cord and kygrid[i] == Ky_FS_cord :
           k_label_FS = klabel[i]
           a = 1
    return k_label_FS


def number_density_func(mu_val,num_den,beta):

   y = 0

   initial_x = -np.pi
   initial_y= -np.pi
   final_x = np.pi
   final_y = np.pi

   dx = 2*np.pi/50
   dy = 2*np.pi/50

   x_grid_cord = np.arange(-np.pi,np.pi,dx)
   y_grid_cord = np.arange(-np.pi,np.pi,dy)

   Ek = 0 #np.zeros((len(x_grid_cord),len(y_grid_cord)))
   nk_noninteracting = 0 #np.zeros((len(x_grid_cord),len(y_grid_cord)))

   t = 1
   bz_area = 0

   for i in range(len(x_grid_cord)):
       for j in range(len(y_grid_cord)):
           Ek = -2*t*(np.cos(x_grid_cord[i])+np.cos(y_grid_cord[j]))
           nk_noninteracting = nk_noninteracting+2*fermi_distribution(Ek,mu_val,beta)
           bz_area = bz_area+1
   
   Nk = nk_noninteracting/bz_area
   y = Nk-num_den

   return y

   
def fermi_momenta_spectral_functions_plot(Text_dir_gf_k_avg_ac_main,Graph_dir_main,n_den,freq_no,sigma,N,u,mu,trot_slice,dtau,omega_max,omega_min):

   beta = float(dtau)*float(trot_slice)
   T = 1/beta
   
   bz,k_pairs = k_point_grid_upper_half_bz(N)
   k_label = bz[:,0]
   Kx_grid = bz[:,1]
   Ky_grid = bz[:,2]

   #Kx_grid = (2*np.pi/int(N))*Kx_grid
   #Ky_grid = (2*np.pi/int(N))*Ky_grid

   K_grid = np.stack((Kx_grid,Ky_grid),axis = 1)
   
   
   T = 1/beta
   t = 1
   Omega_gaus = np.zeros((len(k_label),freq_no))
   A_k_gaus = np.zeros((len(k_label),freq_no))
   FS_gaus = np.zeros((len(k_label)))

   for kk in range(len(k_label)):

       k_point = str(int(k_label[kk]))
       filename_a_k_avg_ac_gaus = '%s/k_point_%s/Spectral_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_gaussian_sigma_%s_nfreq_%s_omega_max_%s_omega_min_%s.dat' %(Text_dir_gf_k_avg_ac_main,k_point,N,u,mu,dtau,trot_slice,k_point,sigma,str(freq_no),omega_max,omega_min)
       omega_gaus,a_k_gaus = np.loadtxt(filename_a_k_avg_ac_gaus,unpack = 'True',usecols = [0,1])

       Omega_gaus[kk,:] = np.copy(omega_gaus)
       A_k_gaus[kk,:] = np.copy(a_k_gaus)
       zero_ind_gaus,zero_val_gaus = find_nearest(omega_gaus,0.00000001)
       FS_gaus[kk] = (a_k_gaus[zero_ind_gaus]+a_k_gaus[zero_ind_gaus-1])/2

   # this holds the co-ordinates of the k points along pi,pi direction
   FS_gaus_pi_pi = []
   k_label_pi_pi = []
   for kk in range(len(k_label)):
       if Kx_grid[kk] == Ky_grid[kk]:
          FS_gaus_pi_pi.append(FS_gaus[kk])
          k_label_pi_pi.append(k_label[kk])
   
   #print(k_label_pi_pi,"k label pi pi")
   k_pi_max = k_label_pi_pi[np.argmax(FS_gaus_pi_pi)]
   
   #print(k_pi_max, "k_pi_max")
   if k_pi_max not in k_label_pi_pi:
      print("error")
      
   Graph_dir_gaus_fs_kf = "%s/Gaussian_sigma_%s_spectral_functions_fermi_momenta"%(Graph_dir_main,sigma)
   if not os.path.exists(Graph_dir_gaus_fs_kf):
      os.makedirs(Graph_dir_gaus_fs_kf)
      
   plt.figure(figsize = (25,20))
   plt.title(r"$A(k_F,\omega)$ along $(\pi,\pi), N = %s \times %s,  U = %s, \mu = %s, n = %s, T = %s$"%(N,N,u,mu,str(round(n_den,3)),str(round(T,3))),fontsize = 40)
   plt.xticks(fontsize = 40)
   plt.yticks(fontsize = 40)
   plt.xlabel(r"$\omega$",fontsize = 40)
   plt.ylabel(r"$A(k_F,\omega)$",fontsize = 40)
   plt.plot(Omega_gaus[k_pi_max],A_k_gaus[k_pi_max],label = r"$k_x=\frac{2\pi}{N} %s, k_y = \frac{2\pi}{N} %s$"%(str(Kx_grid[k_pi_max]),str(Ky_grid[k_pi_max])))
   plt.grid('True',which='both')
   plt.legend(loc = 'best',fontsize = 40)
   plt.savefig("%s/Gaussian_sigma_%s_spectral_functions_fermi_momenta_N_%s_U_%s_mu_%s_dtau_%s_L_%s_omega_max_%s_omega_min_%s.png"%(Graph_dir_gaus_fs_kf,sigma,N,u,mu,dtau,trot_slice,omega_max,omega_min))
          
   plt.close('all')

def fermi_surface_spectral_functions_full_bz_plot_v2(Text_dir_gf_k_avg_ac_main,Graph_dir_main,n_den,freq_no,sigma,N,u,mu,trot_slice,dtau,omega_max,omega_min,x_grid_size,y_grid_size):


   beta = float(dtau)*float(trot_slice)
   T = 1/beta
   
   mu_ni_actual = fsolve(number_density_func,[0.0],args = (n_den,beta))
   nd_ni_actual = number_density_calc(mu_ni_actual[0],beta)

   print(nd_ni_actual,n_den,"Non_interacting, interacting density")
   bz,k_pairs = k_point_grid_upper_half_bz(N)
   k_label_0 = bz[:,0]
   Kx_grid_0 = bz[:,1]
   Ky_grid_0 = bz[:,2]

   k_label_1 = []
   Kx_grid_1 = []
   Ky_grid_1 = []


   for ii in range(len(k_label_0)):

       if(Kx_grid_0[ii] != Ky_grid_0[ii]):
          kx = Ky_grid_0[ii]
          ky = Kx_grid_0[ii]
          Kx_grid_1.append(kx)
          Ky_grid_1.append(ky)
          k_label_1.append(k_label_0[ii])


   k_label = np.append(k_label_0,np.array(k_label_1))
   Kx_grid_2 = np.append(Kx_grid_0,np.array(Kx_grid_1))
   Ky_grid_2 = np.append(Ky_grid_0,np.array(Ky_grid_1))

   Kx_grid = (2*np.pi/int(N))*Kx_grid_2
   Ky_grid = (2*np.pi/int(N))*Ky_grid_2

   K_grid = np.stack((Kx_grid,Ky_grid),axis = 1)
   
   initial_x = np.nanmin(Kx_grid)
   initial_y= np.nanmin(Ky_grid)
   final_x = np.nanmax(Kx_grid)
   final_y = np.nanmax(Ky_grid)

   x_grid_cord = np.linspace(initial_x,final_x,num = x_grid_size)
   y_grid_cord = np.linspace(initial_y,final_y,num = y_grid_size)
   Kx_mesh,Ky_mesh = np.meshgrid(x_grid_cord,y_grid_cord,indexing = 'xy')
   R,Phi = radius_theta_cons(Kx_mesh,Ky_mesh)
   
   dx = Kx_mesh[0,:]
   dy = Ky_mesh[:,0]

   x_grid_cord_ni = np.linspace(0,np.pi,50)
   y_grid_cord_ni = np.linspace(0,np.pi,50)

   T = 1/beta
   t = 1
   Kx_mesh_ni,Ky_mesh_ni = np.meshgrid(x_grid_cord_ni,y_grid_cord_ni,indexing = 'xy')
   A_k_NI = np.zeros(Kx_mesh_ni.shape)
   Ek_NI =  np.zeros(Kx_mesh_ni.shape)
   for i in range(len(Ky_mesh_ni[:,0])):
       for j in range(len(Kx_mesh_ni[0,:])):
           Ek_NI[j][i] = -2*t*(np.cos(Kx_mesh_ni[j][i])+np.cos(Ky_mesh_ni[j][i]))
           Epsilon_kk = Ek_NI[j][i]-mu_ni_actual[0]
           A_k_NI[j][i] = 0.5/(np.cosh((beta*Epsilon_kk)/2))
      
   n_k = np.zeros(len(k_label))
   G_k_b2 = np.zeros(len(k_label))
   Omega_gaus = np.zeros((len(k_label),freq_no))
   A_k_gaus = np.zeros((len(k_label),freq_no))
   FS_gaus = np.zeros((len(k_label)))

   for kk in range(len(k_label)):

       k_point = str(int(k_label[kk]))
       filename_gf_k_avg = '%s/k_point_%s/Retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_ac.dat' %(Text_dir_gf_k_avg_ac_main,k_point,N,u,mu,dtau,trot_slice,k_point)
       tau,g_k,g_k_std = np.loadtxt(filename_gf_k_avg,unpack = 'True',usecols = [0,1,2])
       n_k[kk] = 1-g_k[0]
       lp = int(int(trot_slice)/2)
       G_k_b2[kk] = -beta*g_k[lp]/np.pi
       
       filename_a_k_avg_ac_gaus = '%s/k_point_%s/Spectral_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_gaussian_sigma_%s_nfreq_%s_omega_max_%s_omega_min_%s.dat' %(Text_dir_gf_k_avg_ac_main,k_point,N,u,mu,dtau,trot_slice,k_point,sigma,str(freq_no),omega_max,omega_min)
       omega_gaus,a_k_gaus = np.loadtxt(filename_a_k_avg_ac_gaus,unpack = 'True',usecols = [0,1])

       Omega_gaus[kk,:] = np.copy(omega_gaus)
       A_k_gaus[kk,:] = np.copy(a_k_gaus)
       zero_ind_gaus,zero_val_gaus = find_nearest(omega_gaus,0.00000001)
       FS_gaus[kk] = (a_k_gaus[zero_ind_gaus]+a_k_gaus[zero_ind_gaus-1])/2


   grid_n_k_c = griddata(K_grid,n_k,(Kx_mesh,Ky_mesh),method = 'cubic')
   grid_g_k_c = griddata(K_grid,G_k_b2,(Kx_mesh,Ky_mesh),method = 'cubic')
   mdf_actual,ind_min = find_nearest(n_k,0.5)
   
   dx_ni = np.pi/50
   dy_ni = np.pi/50

   del_x_grid_A_k_NI = np.gradient(A_k_NI,dx_ni,axis = 1)
   del_y_grid_A_k_NI = np.gradient(A_k_NI,dy_ni,axis = 0)
   del_x_grid_g_k_c = np.gradient(grid_g_k_c,dx,axis = 1)
   del_y_grid_g_k_c = np.gradient(grid_g_k_c,dy,axis = 0)
   del_x_grid_mdf_c = np.gradient(grid_n_k_c,dx,axis = 1)
   del_y_grid_mdf_c = np.gradient(grid_n_k_c,dy,axis = 0)

   radial_del_NI = np.zeros(A_k_NI.shape)
   radial_mdf_del = np.zeros(grid_n_k_c.shape)
   radial_g_k_del = np.zeros(grid_g_k_c.shape)


   for i in range(len(Ky_mesh[:,0])):
       for j in range(len(Kx_mesh[0,:])):
           radial_g_k_del[j][i] = (np.cos(Phi[j][i]))*(del_x_grid_g_k_c[j][i])+(np.sin(Phi[j][i]))*(del_y_grid_g_k_c[j][i])
           radial_mdf_del[j][i] = (np.cos(Phi[j][i]))*(del_x_grid_mdf_c[j][i])+(np.sin(Phi[j][i]))*(del_y_grid_mdf_c[j][i])
   
   R_ni,Phi_ni = radius_theta_cons(Kx_mesh_ni,Ky_mesh_ni)
   for i in range(len(Ky_mesh_ni[:,0])):
       for j in range(len(Kx_mesh_ni[0,:])):
           radial_del_NI[j][i] = (np.cos(Phi_ni[j][i]))*(del_x_grid_A_k_NI[j][i])+(np.sin(Phi_ni[j][i]))*(del_y_grid_A_k_NI[j][i])

   
   #===================Calculating FS contours for gaussian spectral functions====================================================================
   

   grid_A_k_c = griddata(K_grid,FS_gaus,(Kx_mesh,Ky_mesh),method = 'cubic')
   max_a_k = np.nanmax(FS_gaus)
   
   del_x_grid_A_k_c = np.gradient(grid_A_k_c,dx,axis = 1)
   del_y_grid_A_k_c = np.gradient(grid_A_k_c,dy,axis = 0)

   radial_del = np.zeros(grid_A_k_c.shape)
   
   for i in range(len(Ky_mesh[:,0])):
       for j in range(len(Kx_mesh[0,:])):

           radial_del[j][i] = (np.cos(Phi[j][i]))*(del_x_grid_A_k_c[j][i])+(np.sin(Phi[j][i]))*(del_y_grid_A_k_c[j][i])
             
   cm = matplotlib.colormaps.get_cmap('plasma')
   
   
   Graph_dir_gaus_fs_scat = "%s/Gaussian_sigma_%s_spectral_functions_scattered"%(Graph_dir_main,sigma)
   if not os.path.exists(Graph_dir_gaus_fs_scat):
      os.makedirs(Graph_dir_gaus_fs_scat)

   Graph_dir_gaus_fs_cont = "%s/Gaussian_sigma_%s_spectral_functions_cubic_interpolation_contour"%(Graph_dir_main,sigma)
   if not os.path.exists(Graph_dir_gaus_fs_cont):
      os.makedirs(Graph_dir_gaus_fs_cont)

   plt.figure(figsize = (20,20))
   plt.title(r"$G(k,\tau=\beta/2), N = %sx%s, U = %s, \mu = %s, n=%s, T = %s$"%(N,N,u,mu,str(round(n_den,4)),str(round(T,3))),fontsize = 40)
   ax = plt.gca()
   plt.xlabel('Kx',fontsize = 40)
   plt.ylabel('Ky',fontsize = 40)
   scat = plt.scatter(Kx_grid,Ky_grid,c=G_k_b2,marker = "s",s = 40000,cmap=cm)
   plt.xticks(fontsize = 40)
   plt.yticks(fontsize = 40)
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   cbar = plt.colorbar(scat, cax= cax)
   cbar.ax.tick_params(labelsize = 30)
   plt.savefig("%s/G_k_beta_2_fermi_surface_N_%s_U_%s_mu_%s_dtau_%s_L_%s_omega_max_%s_omega_min_%s_scatter.png"%(Graph_dir_gaus_fs_scat,N,u,mu,dtau,trot_slice,omega_max,omega_min))
   plt.close('all')
       
   plt.figure(figsize = (20,20))
   plt.title(r"$A(k,\omega=0), N = %sx%s, U = %s, \mu = %s, n=%s, T = %s$"%(N,N,u,mu,str(round(n_den,4)),str(round(T,3))),fontsize = 40)
   ax = plt.gca()
   plt.xlabel('Kx',fontsize = 40)
   plt.ylabel('Ky',fontsize = 40)
   scat = plt.scatter(Kx_grid,Ky_grid,c=FS_gaus[:],marker = "s",s = 40000,cmap=cm)
   plt.xticks(fontsize = 40)
   plt.yticks(fontsize = 40)
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   cbar = plt.colorbar(scat, cax= cax)
   cbar.ax.tick_params(labelsize = 30)
   plt.savefig("%s/Gaussian_sigma_%s_spectral_functions_fermi_surface_N_%s_U_%s_mu_%s_dtau_%s_L_%s_omega_max_%s_omega_min_%s_scatter.png"%(Graph_dir_gaus_fs_scat,sigma,N,u,mu,dtau,trot_slice,omega_max,omega_min))
   plt.close('all')
   
   plt.figure(figsize = (20,20))
   plt.title(r"$A(k,\omega=0), N = %sx%s, U = %s, \mu=%s, n =%s, T=%s$" %(N,N,u,mu,str(round(n_den,4)),str(round(T,3))),fontsize = 40)
   ax = plt.gca()
   plt.xlabel(r"$K_x$",fontsize = 40)
   plt.ylabel(r"$K_y$",fontsize = 40)
   plt.xticks(fontsize = 40)
   plt.yticks(fontsize = 40)
   im = plt.contourf(Kx_mesh,Ky_mesh,grid_A_k_c,levels=50,cmap = cm)
   CS = ax.contour(Kx_mesh, Ky_mesh, radial_del, [0],colors=['black'],linewidths = [5])
   CS_3 = ax.contour(Kx_mesh_ni, Ky_mesh_ni, Ek_NI, [mu_ni_actual], colors = ['green'],linewidths = [5])
   CS_2 = ax.contour(Kx_mesh, Ky_mesh, radial_g_k_del, [0.0], colors = ['white'], linestyles = ['dashdot'],linewidths = [5])
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   cbar = plt.colorbar(im, cax= cax)
   cbar.ax.tick_params(labelsize = 30)
   plt.grid('True')
   plt.savefig("%s/Gaussian_sigma_%s_spectral_functions_fermi_surface_N_%s_U_%s_mu_%s_dtau_%s_L_%s_omega_max_%s_omega_min_%s.png"%(Graph_dir_gaus_fs_cont,sigma,N,u,mu,dtau,trot_slice,omega_max,omega_min))
   
   plt.close('all')

      
def luttinger_surface_spectral_functions_full_bz_plot_v2(Text_dir_gf_k_avg_ac_main,Graph_dir_main,n_den,freq_no,sigma,N,u,mu,trot_slice,dtau,omega_max,omega_min,x_grid_size,y_grid_size):

   beta = float(dtau)*float(trot_slice)
   T = 1/beta
   
   mu_ni_actual = fsolve(number_density_func,[0.0],args = (n_den,beta))
   nd_ni_actual = number_density_calc(mu_ni_actual[0],beta)

   print(nd_ni_actual,n_den,"Non_interacting, interacting density")
   bz,k_pairs = k_point_grid_upper_half_bz(N)
   k_label_0 = bz[:,0]
   Kx_grid_0 = bz[:,1]
   Ky_grid_0 = bz[:,2]

   k_label_1 = []
   Kx_grid_1 = []
   Ky_grid_1 = []


   for ii in range(len(k_label_0)):

       if(Kx_grid_0[ii] != Ky_grid_0[ii]):
          kx = Ky_grid_0[ii]
          ky = Kx_grid_0[ii]
          Kx_grid_1.append(kx)
          Ky_grid_1.append(ky)
          k_label_1.append(k_label_0[ii])


   k_label = np.append(k_label_0,np.array(k_label_1))
   Kx_grid_2 = np.append(Kx_grid_0,np.array(Kx_grid_1))
   Ky_grid_2 = np.append(Ky_grid_0,np.array(Ky_grid_1))

   Kx_grid = (2*np.pi/int(N))*Kx_grid_2
   Ky_grid = (2*np.pi/int(N))*Ky_grid_2

   K_grid = np.stack((Kx_grid,Ky_grid),axis = 1)
   
   initial_x = np.nanmin(Kx_grid)
   initial_y= np.nanmin(Ky_grid)
   final_x = np.nanmax(Kx_grid)
   final_y = np.nanmax(Ky_grid)

   x_grid_cord = np.linspace(initial_x,final_x,num = x_grid_size)
   y_grid_cord = np.linspace(initial_y,final_y,num = y_grid_size)
   Kx_mesh,Ky_mesh = np.meshgrid(x_grid_cord,y_grid_cord,indexing = 'xy')
   R,Phi = radius_theta_cons(Kx_mesh,Ky_mesh)
   
   dx = Kx_mesh[0,:]
   dy = Ky_mesh[:,0]

   x_grid_cord_ni = np.linspace(0,np.pi,50)
   y_grid_cord_ni = np.linspace(0,np.pi,50)

   T = 1/beta
   t = 1
   Kx_mesh_ni,Ky_mesh_ni = np.meshgrid(x_grid_cord_ni,y_grid_cord_ni,indexing = 'xy')
   A_k_NI = np.zeros(Kx_mesh_ni.shape)
   Ek_NI =  np.zeros(Kx_mesh_ni.shape)
   for i in range(len(Ky_mesh_ni[:,0])):
       for j in range(len(Kx_mesh_ni[0,:])):
           Ek_NI[j][i] = -2*t*(np.cos(Kx_mesh_ni[j][i])+np.cos(Ky_mesh_ni[j][i]))
           Epsilon_kk = Ek_NI[j][i]-mu_ni_actual[0]
           A_k_NI[j][i] = 0.5/(np.cosh((beta*Epsilon_kk)/2))
      
   n_k = np.zeros(len(k_label))
   G_k_b2 = np.zeros(len(k_label))
   Omega_gaus = np.zeros((len(k_label),freq_no))
   A_k_gaus = np.zeros((len(k_label),freq_no))
   Luttinger_surface_gaus = np.zeros((len(k_label)))
   FS_gaus = np.zeros((len(k_label)))

   for kk in range(len(k_label)):

       k_point = str(int(k_label[kk]))
       filename_gf_k_avg = '%s/k_point_%s/Retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_ac.dat' %(Text_dir_gf_k_avg_ac_main,k_point,N,u,mu,dtau,trot_slice,k_point)
       tau,g_k,g_k_std = np.loadtxt(filename_gf_k_avg,unpack = 'True',usecols = [0,1,2])
       n_k[kk] = 1+g_k[0]
       lp = int(int(trot_slice)/2)
       G_k_b2[kk] = -beta*g_k[lp]/np.pi
       
       filename_a_k_avg_ac_gaus = '%s/k_point_%s/Spectral_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_gaussian_sigma_%s_nfreq_%s_omega_max_%s_omega_min_%s.dat' %(Text_dir_gf_k_avg_ac_main,k_point,N,u,mu,dtau,trot_slice,k_point,sigma,str(freq_no),omega_max,omega_min)
       omega_gaus,a_k_gaus = np.loadtxt(filename_a_k_avg_ac_gaus,unpack = 'True',usecols = [0,1])

       Omega_gaus[kk,:] = np.copy(omega_gaus)
       A_k_gaus[kk,:] = np.copy(a_k_gaus)
       zero_ind_gaus,zero_val_gaus = find_nearest(omega_gaus,0.00000001)
       Luttinger_surface_gaus[kk] = -1*integrate.trapz(np.divide(A_k_gaus[kk,:],Omega_gaus[kk,:]),Omega_gaus[kk,:])
       FS_gaus[kk] = (a_k_gaus[zero_ind_gaus]+a_k_gaus[zero_ind_gaus-1])/2

   grid_n_k_c = griddata(K_grid,n_k,(Kx_mesh,Ky_mesh),method = 'cubic')
   grid_g_k_c = griddata(K_grid,G_k_b2,(Kx_mesh,Ky_mesh),method = 'cubic')
   mdf_actual,ind_min = find_nearest(n_k,0.5)
   
   dx_ni = np.pi/50
   dy_ni = np.pi/50

   del_x_grid_A_k_NI = np.gradient(A_k_NI,dx_ni,axis = 1)
   del_y_grid_A_k_NI = np.gradient(A_k_NI,dy_ni,axis = 0)
   del_x_grid_g_k_c = np.gradient(grid_g_k_c,dx,axis = 1)
   del_y_grid_g_k_c = np.gradient(grid_g_k_c,dy,axis = 0)
   del_x_grid_mdf_c = np.gradient(grid_n_k_c,dx,axis = 1)
   del_y_grid_mdf_c = np.gradient(grid_n_k_c,dy,axis = 0)

   radial_del_NI = np.zeros(A_k_NI.shape)
   radial_mdf_del = np.zeros(grid_n_k_c.shape)
   radial_g_k_del = np.zeros(grid_g_k_c.shape)


   for i in range(len(Ky_mesh[:,0])):
       for j in range(len(Kx_mesh[0,:])):
           radial_g_k_del[j][i] = (np.cos(Phi[j][i]))*(del_x_grid_g_k_c[j][i])+(np.sin(Phi[j][i]))*(del_y_grid_g_k_c[j][i])
           radial_mdf_del[j][i] = (np.cos(Phi[j][i]))*(del_x_grid_mdf_c[j][i])+(np.sin(Phi[j][i]))*(del_y_grid_mdf_c[j][i])
   
   R_ni,Phi_ni = radius_theta_cons(Kx_mesh_ni,Ky_mesh_ni)
   for i in range(len(Ky_mesh_ni[:,0])):
       for j in range(len(Kx_mesh_ni[0,:])):
           radial_del_NI[j][i] = (np.cos(Phi_ni[j][i]))*(del_x_grid_A_k_NI[j][i])+(np.sin(Phi_ni[j][i]))*(del_y_grid_A_k_NI[j][i])

   
   #===================Calculating Fermi surface and Luttinger surface contours for gaussian spectral functions====================================================================
   
   grid_A_k_c = griddata(K_grid,FS_gaus,(Kx_mesh,Ky_mesh),method = 'cubic')
   max_a_k = np.nanmax(FS_gaus)
   
   del_x_grid_A_k_c = np.gradient(grid_A_k_c,dx,axis = 1)
   del_y_grid_A_k_c = np.gradient(grid_A_k_c,dy,axis = 0)

   radial_del = np.zeros(grid_A_k_c.shape)
   
   for i in range(len(Ky_mesh[:,0])):
       for j in range(len(Kx_mesh[0,:])):

           radial_del[j][i] = (np.cos(Phi[j][i]))*(del_x_grid_A_k_c[j][i])+(np.sin(Phi[j][i]))*(del_y_grid_A_k_c[j][i])
           
           
   grid_Lutt_k_c = griddata(K_grid,Luttinger_surface_gaus,(Kx_mesh,Ky_mesh),method = 'cubic')
           
   cm = matplotlib.colormaps.get_cmap('plasma')
   
   
   Graph_dir_gaus_fs_scat = "%s/Gaussian_sigma_%s_luttinger_surface_spectral_functions_scattered"%(Graph_dir_main,sigma)
   if not os.path.exists(Graph_dir_gaus_fs_scat):
      os.makedirs(Graph_dir_gaus_fs_scat)

   Graph_dir_gaus_fs_cont = "%s/Gaussian_sigma_%s_luttinger_surface_spectral_functions_cubic_interpolation_contour"%(Graph_dir_main,sigma)
   if not os.path.exists(Graph_dir_gaus_fs_cont):
      os.makedirs(Graph_dir_gaus_fs_cont)

   Graph_dir_gaus_fs_cont_data_fit = "%s/Gaussian_sigma_%s_luttinger_surface_spectral_functions_cubic_interpolation_contour/Polynomial_contour_fit"%(Graph_dir_main,sigma)
   if not os.path.exists(Graph_dir_gaus_fs_cont_data_fit):
      os.makedirs(Graph_dir_gaus_fs_cont_data_fit)

   plt.figure(figsize = (20,20))
   plt.title(r"$G(k,\tau=\beta/2), N = %sx%s, U = %s, \mu = %s, n=%s, T = %s$"%(N,N,u,mu,str(round(n_den,4)),str(round(T,3))),fontsize = 40)
   ax = plt.gca()
   plt.xlabel('Kx',fontsize = 40)
   plt.ylabel('Ky',fontsize = 40)
   scat = plt.scatter(Kx_grid,Ky_grid,c=G_k_b2,marker = "s",s = 40000,cmap=cm)
   plt.xticks(fontsize = 40)
   plt.yticks(fontsize = 40)
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   cbar = plt.colorbar(scat, cax= cax)
   cbar.ax.tick_params(labelsize = 30)
   plt.savefig("%s/G_k_beta_2_fermi_surface_N_%s_U_%s_mu_%s_dtau_%s_L_%s_omega_max_%s_omega_min_%s_scatter.png"%(Graph_dir_gaus_fs_scat,N,u,mu,dtau,trot_slice,omega_max,omega_min))
   plt.close('all')
       
   plt.figure(figsize = (20,20))
   plt.title(r"$Re G(k,\omega=0), N = %sx%s, U = %s, \mu = %s, n=%s, T = %s$"%(N,N,u,mu,str(round(n_den,4)),str(round(T,3))),fontsize = 40)
   ax = plt.gca()
   plt.xlabel('Kx',fontsize = 40)
   plt.ylabel('Ky',fontsize = 40)
   scat = plt.scatter(Kx_grid,Ky_grid,c=Luttinger_surface_gaus,marker = "s",s = 40000,cmap=cm)
   plt.xticks(fontsize = 40)
   plt.yticks(fontsize = 40)
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   cbar = plt.colorbar(scat, cax= cax)
   cbar.ax.tick_params(labelsize = 30)
   plt.savefig("%s/Gaussian_sigma_%s_spectral_functions_luttinger_surface_N_%s_U_%s_mu_%s_dtau_%s_L_%s_omega_max_%s_omega_min_%s_scatter.png"%(Graph_dir_gaus_fs_scat,sigma,N,u,mu,dtau,trot_slice,omega_max,omega_min))
   plt.close('all')
   
   plt.figure(figsize = (28,25))
   #plt.title(r"$Re G(k,\omega=0), N = %sx%s, U = %s, \mu=%s, n =%s, T=%s$" %(N,N,u,mu,str(round(n_den,4)),str(round(T,3))),fontsize = 40)
   ax = plt.gca()
   plt.xlabel(r"$K_x$",fontsize = 80)
   plt.ylabel(r"$K_y$",fontsize = 80)
   plt.xticks(fontsize = 80)
   plt.yticks(fontsize = 80)
   im = plt.contourf(Kx_mesh,Ky_mesh,grid_Lutt_k_c,levels=50,cmap = cm)
   CS = ax.contour(Kx_mesh, Ky_mesh, grid_Lutt_k_c, [0],colors=['black'],linewidths = [5])
   #CS_1 = ax.contour(Kx_mesh, Ky_mesh, radial_del, [0.0], colors = ['white'], linestyles = ['dashdot'],linewidths = [5])
   #CS_2 = ax.contour(Kx_mesh, Ky_mesh, radial_g_k_del, [0.0], colors = ['orange'], linestyles = ['dashdot'],linewidths = [5])
   CS_3 = ax.contour(Kx_mesh_ni, Ky_mesh_ni, Ek_NI, [mu_ni_actual], colors = ['white'], linestyles = ['dashdot'],linewidths = [5])
   divider = make_axes_locatable(ax)
   cax = divider.append_axes("right", size="5%", pad=0.05)
   cbar = plt.colorbar(im, cax= cax)
   cbar.ax.tick_params(labelsize = 60)
   plt.grid('True')
   plt.tight_layout()
   plt.savefig("%s/Gaussian_sigma_%s_spectral_functions_luttinger_surface_N_%s_U_%s_mu_%s_dtau_%s_L_%s_omega_max_%s_omega_min_%s.png"%(Graph_dir_gaus_fs_cont,sigma,N,u,mu,dtau,trot_slice,omega_max,omega_min))
   
   print(len(CS.allsegs),CS.allsegs)
   

   if not any(CS.allsegs):
      fs_area = 0.0
      plt.figure(figsize = (20,20))
      plt.xlabel('Kx',fontsize = 40)
      plt.ylabel('Ky',fontsize = 40)
      plt.xticks(fontsize = 40)
      plt.yticks(fontsize = 40)
      plt.xlim(0,np.pi)
      plt.ylim(0,np.pi)
      im = plt.contourf(Kx_mesh,Ky_mesh,grid_Lutt_k_c,levels=50,cmap = cm)
      plt.grid('True')
      plt.legend(loc = 'best', fontsize = 40)
      plt.savefig('%s/Luttinger_surface_contour_only_N_%s_U_%s_mu_%s_L_%s_cubic_interpolation_data_fit.png'%(Graph_dir_gaus_fs_cont_data_fit,N,u,mu,trot_slice))
      plt.close('all')
   
   else:
      
      dat0 = CS.allsegs[0][0]
      Poly_interacting = np.polynomial.polynomial.Polynomial(np.polynomial.polynomial.polyfit(dat0[:,0],dat0[:,1],8))
      Kx_fit_intr = np.linspace(np.nanmin(dat0[:,0]),np.nanmax(dat0[:,0]),num = 101)
      
      plt.figure(figsize = (20,20))
      plt.xlabel('Kx',fontsize = 40)
      plt.ylabel('Ky',fontsize = 40)
      plt.xticks(fontsize = 40)
      plt.yticks(fontsize = 40)
      plt.xlim(0,np.pi)
      plt.ylim(0,np.pi)
      im = plt.contourf(Kx_mesh,Ky_mesh,grid_Lutt_k_c,levels=50,cmap = cm)
      plt.plot(dat0[:,0],dat0[:,1],'o',markersize = 10, label = "Luttinger surface")
      plt.plot(Kx_fit_intr, Poly_interacting(Kx_fit_intr),linewidth = 5)
      plt.grid('True')
      plt.legend(loc = 'best', fontsize = 40)
      plt.savefig('%s/Luttinger_surface_contour_only_N_%s_U_%s_mu_%s_L_%s_cubic_interpolation_data_fit.png'%(Graph_dir_gaus_fs_cont_data_fit,N,u,mu,trot_slice))
      plt.close('all')
   
      if Kx_fit_intr[0] == 0:
         fs_area = integrate.trapz(Poly_interacting(Kx_fit_intr), Kx_fit_intr)
    
      else:
         fs_area = np.pi*Kx_fit_intr[0]+integrate.trapezoid(Poly_interacting(Kx_fit_intr), Kx_fit_intr)
      
   dat_ni = CS_3.allsegs[0][0]
      
   Poly_noninteracting = np.polynomial.polynomial.Polynomial(np.polynomial.polynomial.polyfit(dat_ni[:,0],dat_ni[:,1],4))
   Kx_fit_ni = np.linspace(np.nanmin(dat_ni[:,0]),np.nanmax(dat_ni[:,0]),num = 101)
      
   if Kx_fit_ni[0] == 0:
      fs_area_ni = integrate.trapz(Poly_noninteracting(Kx_fit_ni), Kx_fit_ni)
      
   else:
      fs_area_ni = np.pi*Kx_fit_ni[0]+integrate.trapezoid(Poly_noninteracting(Kx_fit_ni), Kx_fit_ni)
         
         
   return fs_area,fs_area_ni
   
   
   
def main():
    #if(total!=6):
    #    raise ValueError('missing args')

    N = 12
    U = "10.0"
    Mu = ["0.00","0.20","0.40","0.60","0.80","1.00","1.20","1.40","1.60","1.80","2.00","2.20","2.40","2.60","2.80","3.00","3.20","3.40","3.80","4.00","4.40","4.60","5.00","5.20","5.40","5.60","6.00"] #,"6.20","6.40"]
    Trot = ["10","12","14","16","18","20"] #,"30","40","50","60","70","80"]
    FS_area = np.zeros((len(Mu),len(Trot)))
    FS_area_NI = np.zeros((len(Mu),len(Trot)))
    Nden = np.zeros((len(Mu),len(Trot)))
    
    Dtau = "0.025"
    Omega_max = "25.0"
    Omega_min = "-25.0"
    Freq_no = 2500
    
    XGrid_size = 20
    YGrid_size = 20
    Sigma = ["2.0"]
    
    
    for sig in range(len(Sigma)):
        
        Graph_dir_main = "/Users/roy.369/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Graphs/Graphs_N_%s_maxent/Graphs_N_%s_U_%s_dtau_%s_maxent_spectral_function_plots/Graphs_N_%s_U_%s_dtau_%s_maxent_plots_fermi_surface/Spectral_functions_fermi_surface_gaussian_sigma_%s_default_model_omega_max_%s_omega_min_%s"%(N,N,U,Dtau,N,U,Dtau,Sigma[sig],Omega_max,Omega_min)
        
        Graph_dir_main_kf = "/Users/roy.369/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Graphs/Graphs_N_%s_maxent/Graphs_N_%s_U_%s_dtau_%s_maxent_spectral_function_plots/Graphs_N_%s_U_%s_dtau_%s_maxent_plots_fermi_surface/Spectral_functions_fermi_momenta_gaussian_sigma_%s_default_model_omega_max_%s_omega_min_%s"%(N,N,U,Dtau,N,U,Dtau,Sigma[sig],Omega_max,Omega_min)
        if not os.path.exists(Graph_dir_main_kf):
           os.makedirs(Graph_dir_main_kf)
           
        for k in range(len(Trot)):
        
            Graph_dir_ac_fs = "%s/dtau_%s_L_%s"%(Graph_dir_main,Dtau,Trot[k])
            if not os.path.exists(Graph_dir_ac_fs):
                   os.makedirs(Graph_dir_ac_fs)
            
            Graph_dir_ac_kf = "%s/dtau_%s_L_%s"%(Graph_dir_main_kf,Dtau,Trot[k])
            if not os.path.exists(Graph_dir_ac_kf):
                   os.makedirs(Graph_dir_ac_kf)
            
            for i in range(len(Mu)):
        

       
       
                Text_dir_gf_k_avg_ac_main = "/Users/roy.369/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Text_files/Text_files_N_%s_analytic_continuation_spectral_functions/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Spectral_functions_gaussian_sigma_%s_nfreq_%s"%(N,N,U,Dtau,Mu[i],Dtau,Trot[k],Sigma[sig],str(Freq_no))
                
                Text_dir_eqm = "/Users/roy.369/Documents/Luttinger_theorem/FHM_Legacy_data/FHM_Roy_data/Text_files/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements"%(N,N,U,Dtau,Mu[i],Dtau,Trot[k])
                filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U,Mu[i],Dtau,Trot[k])
                with open(filename_eqm_avg, 'rb') as infile:
                    sys_measure_avg = pickle.load(infile)

                filename_eqm_std = '%s/Thermodynamic_measurements_normal_standard_deviation_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U,Mu[i],Dtau,Trot[k])
                with open(filename_eqm_std, 'rb') as infile:
                    sys_measure_std = pickle.load(infile)

                Density = sys_measure_avg['Density averaged']
                Density_dev = sys_measure_std['Density standard deviation']
                Avg_sign = sys_measure_avg['Total sign averaged']
                #print("avg sign, density, mu", Avg_sign, Density, Mu[i])
    
                fermi_momenta_spectral_functions_plot(Text_dir_gf_k_avg_ac_main,Graph_dir_ac_kf,Density,Freq_no,Sigma[sig],N,U,Mu[i],Trot[k],Dtau,Omega_max,Omega_min)

                #fermi_surface_spectral_functions_full_bz_plot_v2(Text_dir_gf_k_avg_ac_main,Graph_dir_ac_fs,Density,Freq_no,Sigma[sig],N,U,Mu[i],Trot[k],Dtau,Omega_max,Omega_min,XGrid_size,YGrid_size)
                Nden[i][k] = Density
                FS_area[i][k], FS_area_NI[i][k] = luttinger_surface_spectral_functions_full_bz_plot_v2(Text_dir_gf_k_avg_ac_main,Graph_dir_ac_fs,Density,Freq_no,Sigma[sig],N,U,Mu[i],Trot[k],Dtau,Omega_max,Omega_min,XGrid_size,YGrid_size)

        
        nden=np.linspace(1.02,1.45,num=101)
        plt.figure(figsize = (30,20))
        plt.title(r"Luttinger surface volume, $N=%s \times %s, U = %s$"%(N,N,U),fontsize = 40)
        plt.xlabel(r"n",fontsize = 40)
        plt.ylabel(r"$\sum_{k} \Theta [Re G(k,\omega=0)]$",fontsize = 40)
        plt.xticks(fontsize = 40)
        plt.yticks(fontsize = 40)
        for k in range(len(Trot)):
            T_val = 1/(float(Dtau)*float(Trot[k]))
            plt.plot(Nden[:,k],2*FS_area[:,k]/(np.pi*np.pi),marker="o",markersize = 15,label = r"T=%s"%(str(round(T_val,3))))
        plt.plot(nden,nden,linestyle="dashed",linewidth=5)
        plt.plot(nden,nden-1,linestyle="dashed",linewidth=5)
        plt.grid('True',which = 'both')
        plt.xlim(1.00,1.25)
        plt.legend(loc='best',fontsize = 40)
        plt.savefig("%s/Gaussian_sigma_%s_default_model_Luttinger_surface_area_N_%s_U_%s_dtau_%s.png"%(Graph_dir_main,Sigma[sig],N,U,Dtau))
        
        #plt.figure(figsize = (30,20))
        #plt.title(r"Luttinger surface volume, $N=%s \times %s, U = %s$"%(N,N,U),fontsize = 40)
        #plt.xlabel(r"n",fontsize = 40)
        #plt.ylabel(r"$\sum_{k} \Theta [Re G(k,\omega=0)]$",fontsize = 40)
        #plt.xticks(fontsize = 40)
        #plt.yticks(fontsize = 40)
        #for k in range(len(Trot)):
        #    T_val = 1/(float(Dtau)*float(Trot[k]))
        #    plt.plot(Nden[:,k],2*FS_area[:,k]/(np.pi*np.pi),marker="o",markersize = 15,label = r"T=%s"%(str(round(T_val,3))))
        #    plt.plot(Nden[:,k],2*FS_area_NI[:,k]/(np.pi*np.pi),linestyle = 'dashed',linewidth = 3) #,label = r"T=%s"%(str(round(T_val,3))))
        #plt.plot(nden,nden,linestyle="dashed",linewidth=5)
        #plt.plot(nden,nden-1,linestyle="dashed",linewidth=5)
        #plt.grid('True',which = 'both')
        #plt.xlim(1.02,1.45)
        #plt.legend(loc='best',fontsize = 40)
        #plt.savefig("%s/Gaussian_sigma_%s_default_model_Luttinger_surface_non_interacting_volume_compare_area_N_%s_U_%s_dtau_%s.png"%(Graph_dir_main,Sigma[sig],N,U,Dtau))
main()        
#if __name__ == '__main__':
#    sys.argv
#    total = len(sys.argv)
#    print("No of sys arguments",total)
#    cmdargs = sys.argv
#    main(total,cmdargs)








