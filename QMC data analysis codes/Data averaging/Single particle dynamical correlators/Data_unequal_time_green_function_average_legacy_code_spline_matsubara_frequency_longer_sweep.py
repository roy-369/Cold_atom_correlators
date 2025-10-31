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
from scipy import integrate
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
from scipy.interpolate import griddata
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.colorbar
from matplotlib import rc
mpl.rcParams['axes.linewidth'] = 5
from scipy import integrate


def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return array[idx],idx



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

############################################## Jackknife statistics of mean and error estimation ##############################################################################

def jackknife_sampling(data):
    """
    Perform jackknife resampling to estimate the mean and variance of a dataset.
    
    Parameters:
        data (array-like): Input array of data points.
    
    Returns:
        tuple: (jackknife mean, jackknife variance)
    """
    n = len(data)
    if n <= 1:
        raise ValueError("Data array must have more than one element.")
    
    # Create jackknife samples
    jackknife_samples = [np.delete(data, i) for i in range(n)]
    
    # Calculate statistics for each jackknife sample
    sample_means = np.array([np.mean(sample) for sample in jackknife_samples])
    
    # Calculate jackknife mean
    jackknife_mean = np.mean(sample_means)
    
    # Calculate jackknife variance, using formula from Wikipedia, var(jackknife) = (n-1)/n \sum_{i}(sample_mean(i)-jackknife_mean)^2
    jackknife_variance = (n - 1) * np.var(sample_means, ddof=0)
    
    return jackknife_mean, np.sqrt(jackknife_variance)



############################################## Calculating Matsubara frequency space data from averages by performing spline interpolation#####################################


def scipy_spline_fourier_transform(dtau,L,Tau,G_tau,G_tau_std,n_max,filename_fit):
    
    beta = float(dtau)*int(L)

    delta = beta/L
    ##G_tau[-1] = -(1+G_tau[0])
    Tau_fit = np.linspace(0,beta, num = 101)
    cs = CubicSpline(Tau, G_tau)
    G_fit = cs(Tau_fit)
    
    fitdata = np.stack((Tau_fit,G_fit),axis = 1)
    np.savetxt(filename_fit,fitdata)

    coefficient = cs.c
                
    a = np.zeros(L)
    b = np.zeros(L)
    c = np.zeros(L)
    d = np.zeros(L)

    a_var = np.zeros(L)
    b_var = np.zeros(L)
    c_var = np.zeros(L)
    d_var = np.zeros(L)

    for j in range(L):
        
        a[j] = coefficient[3][j]
        b[j] = coefficient[2][j]
        c[j] = coefficient[1][j]
        d[j] = coefficient[0][j]
    
    Omega = np.zeros(n_max+1)
    G_omega = np.zeros(n_max+1,dtype = np.complex128)
    G_omega_real_var = np.zeros(n_max+1)
    G_omega_imag_var = np.zeros(n_max+1)

    
    for n in range(n_max+1):
        om = (2*n+1)*np.pi/beta 
        Omega[n] = om
        for j in range(L):

            arg_j = om*delta*j
            arg_jp = om*delta*(j+1)
            
            exp_j = np.exp(1j*om*delta*j)
            exp_jplus = np.exp(1j*om*delta*(j+1))
            
            term_jplus = exp_jplus*(-6*d[j]/(om**4)+(2*1j*c[j]+6*1j*d[j]*delta)/(om**3) \
                                    +(b[j]+2*c[j]*delta+3*d[j]*delta*delta)/(om**2)-1*(1j*a[j]+1j*b[j]*delta+1j*c[j]*(delta**2)+1j*d[j]*(delta**3))/om)
            term_j = exp_j*(1j*a[j]/om-b[j]/(om**2)-1j*2*c[j]/(om**3)+6*d[j]/(om**4))
            
            G_omega[n] = G_omega[n]+term_jplus+term_j

            g_real_err = np.cos(arg_jp)*np.cos(arg_jp)*(36*d_var[j]/(om**8)+(b_var[j]+4*c_var[j]*(delta**2)+9*d_var[j]*(delta**4))/(om**4)) \
                        +np.sin(arg_jp)*np.sin(arg_jp)*((a_var[j]+b_var[j]*(delta**2)+c_var[j]*(delta**4)+d_var[j]*(delta**6))/(om**2)+(4*c_var[j]+36*d_var[j]*(delta**2))/(om**6)) \
                        +np.cos(arg_j)*np.cos(arg_j)*(36*d_var[j]/(om**8)+b_var[j]/(om**4))+np.sin(arg_j)*np.sin(arg_j)*(4*c_var[j]/(om**6)+a_var[j]/(om**2))

            g_imag_err = np.sin(arg_jp)*np.sin(arg_jp)*(36*d_var[j]/(om**8)+(b_var[j]+4*c_var[j]*(delta**2)+9*d_var[j]*(delta**4))/(om**4)) \
                        +np.cos(arg_jp)*np.cos(arg_jp)*((a_var[j]+b_var[j]*(delta**2)+c_var[j]*(delta**4)+d_var[j]*(delta**6))/(om**2)+(4*c_var[j]+36*d_var[j]*(delta**2))/(om**6)) \
                        +np.sin(arg_j)*np.sin(arg_j)*(36*d_var[j]/(om**8)+b_var[j]/(om**4))+np.cos(arg_j)*np.cos(arg_j)*(4*c_var[j]/(om**6)+a_var[j]/(om**2))

            G_omega_real_var[n] = G_omega_real_var[n]+g_real_err
            G_omega_imag_var[n] = G_omega_imag_var[n]+g_imag_err

    return Omega,np.real(G_omega), np.imag(G_omega), (1e-6)*np.ones(len(G_omega)), (1e-6)*np.ones(len(G_omega)) #np.sqrt(G_omega_real_var), np.sqrt(G_omega_imag_var)
   

#=====================================================================================================================================================================

#============================== The following function computes the spline interpolation and fourier transform of G(tau) to get G(iomega_n), then compute average over all realizations.===================================

#===============================Self energy Sigma(i\omega_n) is computed from the realization averaged G(iomega_n)=============================================

#===================================================================================================================================================


def green_function_spline_FT_matsubara_averaging_momentum_space_v1(Text_dir_gf_k,Text_dir_gf_k_mf,Text_dir_gf_k_avg_ac,Text_dir_gf_k_mf_avg_ac,Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,n_den,run_no,n_omeg):

   bz, k_pair = k_point_grid_upper_half_bz(N)
   
   k_lab = bz[:,0]
   kx_grid = bz[:,1]
   ky_grid = bz[:,2]
        
   M_1 = np.zeros(k_pair)
   M_2 = np.zeros(k_pair)
   Mean_ac = np.zeros(k_pair)
   Variance_ac = np.zeros(k_pair)

   Norm_self_energy = np.zeros(k_pair)
   Sig0 = 0.5*(n_den-1)*float(U)

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt
       
   Text_dir_spline_fit = "%s/Cubic_spline_fit"%Text_dir_gf_k
   if not os.path.exists(Text_dir_spline_fit):
       os.makedirs(Text_dir_spline_fit)


   for k in range(len(k_lab)):

       k_dir_mf = "%s/k_point_%s"%(Text_dir_gf_k_mf,str(k))
       if not os.path.exists(k_dir_mf):
          os.makedirs(k_dir_mf)

       k_dir_fit = "%s/k_point_%s"%(Text_dir_spline_fit,str(k))
       if not os.path.exists(k_dir_fit):
          os.makedirs(k_dir_fit)

       k_pt = k_lab[k]
       kx = (2*np.pi/int(N))*kx_grid[k]
       ky = (2*np.pi/int(N))*ky_grid[k]

       ep_k = -2*(np.cos(kx)+np.cos(ky))-float(Mu)
       m1 = ep_k-float(Mu)-0.5*float(U)+0.5*n_den*float(U)
       m2 = (ep_k-float(Mu)-0.5*float(U))**2+float(U)*(ep_k-float(Mu)-0.5*float(U))*n_den+0.5*float(U)*float(U)*n_den
       M_1[k] = m1
       M_2[k] = m2
       Mean_ac[k] = m1
       Variance_ac[k] = 2*np.sqrt(m2-m1*m1)
       Norm_self_energy[k] = float(U)*float(U)*(n_den/2)*(1-n_den/2)


       filename_gf_k_0 = '%s/Retarded_green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
       
       if os.path.exists(filename_gf_k_0):

           with open(filename_gf_k_0, 'rb') as infile:
             green_function_k_0 = pickle.load(infile)

           filename_gf_k_variation_0 = '%s/Retarded_green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
           with open(filename_gf_k_variation_0, 'rb') as infile:
             green_function_k_std_0 = pickle.load(infile)

       
           filename_fit = "%s/Spline_fit_G_tau_k_point_%s_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat"%(k_dir_fit,str(k_lab[k]),N,U,Mu,str(dtau),str(L))
           g_k_data = np.copy(green_function_k_0[:,k])
           g_k_std_data = np.copy(green_function_k_std_0[:,k])

           g_k_data[-1] = -1*(1+g_k_data[0])
           g_k_std_data[-1] = g_k_std_data[0]

           omega_n, g_mf_re, g_mf_im, g_mf_re_std, g_mf_im_std = scipy_spline_fourier_transform(float(dtau),int(L),Tau,g_k_data,g_k_std_data,n_omeg,filename_fit)

           Green_function_k_mf_re = g_mf_re.copy()
           Green_function_k_mf_im = g_mf_im.copy()

           g_mf_data = np.stack((omega_n,g_mf_re,g_mf_im),axis=1)
           filename_g_mf = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_point_%s_n_max_%s_r_0.dat'%(k_dir_mf,N,U,Mu,dtau,L,str(k),str(n_omeg))
           np.savetxt(filename_g_mf,g_mf_data)

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

 
                 filename_fit = "%s/Spline_fit_G_tau_k_point_%s_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(k_dir_fit,str(k_lab[k]),N,U,Mu,str(dtau),str(L),realization)
                 g_k_data = np.copy(green_function_k[:,k])
                 g_k_std_data = np.copy(green_function_k_std[:,k])
   
                 g_k_data[-1] = -1*(1+g_k_data[0])
                 g_k_std_data[-1] = g_k_std_data[0]
 
                 omega_n, g_mf_re, g_mf_im, g_mf_re_std, g_mf_im_std = scipy_spline_fourier_transform(float(dtau),int(L),Tau,g_k_data,g_k_std_data,n_omeg,filename_fit)
                 
                 
                 X_Gf_k_mf_re = mat_append(Green_function_k_mf_re,g_mf_re)
                 X_Gf_k_mf_im = mat_append(Green_function_k_mf_im,g_mf_im)
                 
                 Green_function_k_mf_re = X_Gf_k_mf_re.copy()
                 Green_function_k_mf_im = X_Gf_k_mf_im.copy()

                 #======= Saving spline interpolated green function in matsuabara space for each realization =====================

                 g_mf_data = np.stack((omega_n,g_mf_re,g_mf_im),axis=1)
                 filename_g_mf = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_point_%s_n_max_%s_r_%s.dat'%(k_dir_mf,N,U,Mu,dtau,L,str(k),str(n_omeg),realization)
                 np.savetxt(filename_g_mf,g_mf_data)


       #========================================== Saving data using standard sample statistics ===========================================

       print("G_mf_shape",Green_function_k_mf_re.shape)   
       Green_function_mf_re_k_avg = np.mean(Green_function_k_mf_re,axis = 1)
       Green_function_mf_im_k_avg = np.mean(Green_function_k_mf_im,axis = 1)

       Green_function_mf_re_k_std_avg = (1/np.sqrt(len(Green_function_k_mf_re[0,:])))*np.std(Green_function_k_mf_re,axis = 1,ddof=1)
       Green_function_mf_im_k_std_avg = (1/np.sqrt(len(Green_function_k_mf_re[0,:])))*np.std(Green_function_k_mf_im,axis = 1,ddof=1)

       Gf_data = np.stack((omega_n,Green_function_mf_re_k_avg,Green_function_mf_re_k_std_avg,Green_function_mf_im_k_avg,Green_function_mf_im_k_std_avg),axis = 1)
       filename_g_k_point = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_g_k_point,Gf_data)

       Self_en_real_k = np.zeros(len(omega_n))
       Self_en_imag_k = np.zeros(len(omega_n))
       Self_en_real_k_std = np.zeros(len(omega_n))
       Self_en_imag_k_std = np.zeros(len(omega_n))

       for n in range(len(omega_n)):
           g_mf = Green_function_mf_re_k_avg[n]+1j*Green_function_mf_im_k_avg[n]
           sig = 1j*omega_n[n]-ep_k-1/g_mf
           Self_en_real_k[n] = np.real(sig)
           Self_en_imag_k[n] = np.imag(sig)
           del_g = Green_function_mf_re_k_std_avg[n]+1j*Green_function_mf_im_k_std_avg[n]

           ### ==================== Using error propagation through dyson equation by addition/multipication formulas==============================================================================

           ## labeling y = Re G^2+Im G^2 = a^2+b^2
           ## Re Sig = a/y
           ## Im Sig = b/y
           ## Var y = 4* Re G^2*Var Re G + 4* Im G^2 * Var Im G
           ## Var Re Sig = Re Sig^2(Var a/a^2+Var y/y^2)
           ## Var Im Sig = Im Sig^2(Var b/b^2+Var y/y^2)

           y = Green_function_mf_re_k_avg[n]*Green_function_mf_re_k_avg[n]+Green_function_mf_im_k_avg[n]*Green_function_mf_im_k_avg[n]
           y_var = (2*Green_function_mf_re_k_avg[n]*Green_function_mf_re_k_std_avg[n])**2+(2*Green_function_mf_im_k_avg[n]*Green_function_mf_im_k_std_avg[n])**2

           #Self_en_real_k_std[n] = np.sqrt((Green_function_mf_re_k_std_avg[n]*Green_function_mf_re_k_std_avg[n]*y*y+Green_function_mf_re_k_avg[n]*Green_function_mf_re_k_avg[n]*y_var)/(y**4))
           #Self_en_imag_k_std[n] = np.sqrt((Green_function_mf_im_k_std_avg[n]*Green_function_mf_im_k_std_avg[n]*y*y+Green_function_mf_im_k_avg[n]*Green_function_mf_im_k_avg[n]*y_var)/(y**4))

           Self_en_real_k_std[n] = np.sqrt(Self_en_real_k[n]*Self_en_real_k[n]*((Green_function_mf_re_k_std_avg[n]/Green_function_mf_re_k_avg[n])**2+y_var/(y*y)))
           Self_en_imag_k_std[n] = np.sqrt(Self_en_imag_k[n]*Self_en_imag_k[n]*((Green_function_mf_im_k_std_avg[n]/Green_function_mf_im_k_avg[n])**2+y_var/(y*y)))

           ### ==================== Using error propagation through del(Sigma) = (1/G^2)*del(G)

           #del_sig = del_g/(g_mf*g_mf)
           #Self_en_real_k_std[n] = np.real(del_sig)
           #Self_en_imag_k_std[n] = np.imag(del_sig)


       Self_en_real_k_hartree = np.zeros(len(omega_n))
       Self_en_imag_k_hartree = np.zeros(len(omega_n))

       for n in range(len(omega_n)):
           Self_en_real_k_hartree[n] = Self_en_real_k[n]-Self_en_real_k[-1]
           Self_en_imag_k_hartree[n] = Self_en_imag_k[n]
          
       Self_en_real_k_hartree[-1] = 1e-6


       Self_en_data = np.stack((omega_n,Self_en_real_k,Self_en_real_k_std,Self_en_imag_k,Self_en_imag_k_std),axis = 1)
       filename_self_en_k_point = '%s/Spline_interpolated_self_energy_momentum_space_avg_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_k_point,Self_en_data)

       Self_en_hartree_data = np.stack((omega_n,Self_en_real_k_hartree,Self_en_real_k_std,Self_en_imag_k_hartree,Self_en_imag_k_std),axis = 1)
       filename_self_en_hartree_k_point = '%s/Spline_interpolated_hartree_corrected_self_energy_momentum_space_avg_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_hartree_k_point,Self_en_hartree_data)


       kx = (2*np.pi/int(N))*kx_grid[k]
       ky = (2*np.pi/int(N))*ky_grid[k]

       Gf_data = np.stack((omega_n,Green_function_mf_re_k_avg,Green_function_mf_re_k_std_avg,Green_function_mf_im_k_avg,Green_function_mf_im_k_std_avg),axis = 1)
       filename_g_k_point = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_g_k_point,Gf_data)


       #######================================ Perform jackknife statistics of sample, save data =======================

       G_mf_re_jk = np.zeros(len(omega_n))
       G_mf_re_jk_std = np.zeros(len(omega_n))
       G_mf_im_jk = np.zeros(len(omega_n))
       G_mf_im_jk_std = np.zeros(len(omega_n))

       for n in range(len(omega_n)):
           
           g_mf_re_data = Green_function_k_mf_re[n,:]
           g_mf_im_data = Green_function_k_mf_im[n,:]
           G_mf_re_jk[n],G_mf_re_jk_std[n] = jackknife_sampling(g_mf_re_data)
           G_mf_im_jk[n],G_mf_im_jk_std[n] = jackknife_sampling(g_mf_im_data)

       Gf_jk_data = np.stack((omega_n,G_mf_re_jk,G_mf_re_jk_std,G_mf_im_jk,G_mf_im_jk_std),axis = 1)
       filename_g_k_point_jk = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_jackknife_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_g_k_point_jk,Gf_jk_data)

       Sig_mf_re_jk = np.zeros(len(omega_n))
       Sig_mf_re_jk_hc = np.zeros(len(omega_n))
       Sig_mf_re_jk_std = np.zeros(len(omega_n))
       Sig_mf_im_jk = np.zeros(len(omega_n))
       Sig_mf_im_jk_std = np.zeros(len(omega_n))

       for n in range(len(omega_n)):

           g_mf = G_mf_re_jk[n]+1j*G_mf_im_jk[n]
           sig = 1j*omega_n[n]-ep_k-1/g_mf
           Sig_mf_re_jk[n] = np.real(sig)
           Sig_mf_im_jk[n] = np.imag(sig)

           y_jk = G_mf_re_jk[n]*G_mf_re_jk[n]+G_mf_im_jk[n]*G_mf_im_jk[n]
           y_var_jk = (2*G_mf_re_jk[n]*G_mf_re_jk_std[n])**2+(2*G_mf_im_jk[n]*G_mf_im_jk_std[n])**2

           #Self_en_real_k_std[n] = np.sqrt((Green_function_mf_re_k_std_avg[n]*Green_function_mf_re_k_std_avg[n]*y*y+Green_function_mf_re_k_avg[n]*Green_function_mf_re_k_avg[n]*y_var)/(y**4))
           #Self_en_imag_k_std[n] = np.sqrt((Green_function_mf_im_k_std_avg[n]*Green_function_mf_im_k_std_avg[n]*y*y+Green_function_mf_im_k_avg[n]*Green_function_mf_im_k_avg[n]*y_var)/(y**4))

           Sig_mf_re_jk_std[n] = np.sqrt(Sig_mf_re_jk[n]*Sig_mf_re_jk[n]*((G_mf_re_jk_std[n]/G_mf_re_jk[n])**2+y_var/(y*y)))
           Sig_mf_im_jk_std[n] = np.sqrt(Sig_mf_im_jk[n]*Sig_mf_im_jk[n]*((G_mf_im_jk_std[n]/G_mf_im_jk[n])**2+y_var/(y*y)))

       for n in range(len(omega_n)):
           Sig_mf_re_jk_hc[n] = Sig_mf_re_jk[n]-Sig_mf_re_jk[-1]    #Sig0

       Sig_mf_re_jk_hc[-1] = 1e-6
       Self_en_data_jk = np.stack((omega_n,Sig_mf_re_jk,Sig_mf_re_jk_std,Sig_mf_im_jk,Sig_mf_im_jk_std),axis = 1)
       filename_self_en_k_point_jk = '%s/Spline_interpolated_self_energy_momentum_space_jackknife_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_k_point_jk,Self_en_data_jk)

       Self_en_hartree_data_jk = np.stack((omega_n,Sig_mf_re_jk_hc,Sig_mf_re_jk_std,Sig_mf_im_jk,Sig_mf_im_jk_std),axis = 1)
       filename_self_en_hartree_k_point_jk = '%s/Spline_interpolated_hartree_corrected_self_energy_momentum_space_jackknife_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_hartree_k_point_jk,Self_en_hartree_data_jk)



######=================================================================================================

   np.savetxt("%s/Default_model_mean_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L),Mean_ac)
   np.savetxt("%s/Default_model_variance_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L),Variance_ac)
   np.savetxt("%s/Self_energy_norm_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat"%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L),Norm_self_energy)

#=====================================================================================================================================================================

#============================== The following function computes the spline interpolation and fourier transform of G(tau) to get G(iomega_n), computes Sigma(i\omega_n) from Dyson equation, then compute average over all realizations.===========================

#===================================================================================================================================================


def green_function_spline_FT_matsubara_averaging_momentum_space_v2(Text_dir_gf_k,Text_dir_gf_k_mf,Text_dir_gf_k_avg_ac,Text_dir_gf_k_mf_avg_ac,Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,n_den,run_no,n_omeg):

   bz, k_pair = k_point_grid_upper_half_bz(N)

   k_lab = bz[:,0]
   kx_grid = bz[:,1]
   ky_grid = bz[:,2]

   M_1 = np.zeros(k_pair)
   M_2 = np.zeros(k_pair)
   Mean_ac = np.zeros(k_pair)
   Variance_ac = np.zeros(k_pair)

   Norm_self_energy = np.zeros(k_pair)

   num_den = float(n_den)
   Sig0 = 0.5*float(U)*(num_den-1)

   timeslices = int(L)+1
   Tau = np.zeros(timeslices)
   for tt in range(timeslices):
       Tau[tt] = float(dtau)*tt

   Text_dir_spline_fit = "%s/Cubic_spline_fit"%Text_dir_gf_k
   if not os.path.exists(Text_dir_spline_fit):
       os.makedirs(Text_dir_spline_fit)


   for k in range(len(k_lab)):

       k_dir_mf = "%s/k_point_%s"%(Text_dir_gf_k_mf,str(k))
       if not os.path.exists(k_dir_mf):
          os.makedirs(k_dir_mf)

       k_dir_fit = "%s/k_point_%s"%(Text_dir_spline_fit,str(k))
       if not os.path.exists(k_dir_fit):
          os.makedirs(k_dir_fit)


       k_pt = k_lab[k]
       kx = (2*np.pi/int(N))*kx_grid[k]
       ky = (2*np.pi/int(N))*ky_grid[k]
       ep_k = -2*(np.cos(kx)+np.cos(ky))-float(Mu)

       ep_k = -2*(np.cos(kx)+np.cos(ky))-float(Mu)
       m1 = ep_k-float(Mu)-0.5*float(U)+0.5*n_den*float(U)
       m2 = (ep_k-float(Mu)-0.5*float(U))**2+float(U)*(ep_k-float(Mu)-0.5*float(U))*n_den+0.5*float(U)*float(U)*n_den
       M_1[k] = m1
       M_2[k] = m2
       Mean_ac[k] = m1
       Mean_ac[k] = m1
       Variance_ac[k] = 2*np.sqrt(m2-m1*m1)
       Norm_self_energy[k] = float(U)*float(U)*(n_den/2)*(1-n_den/2)


       filename_gf_k_0 = '%s/Retarded_green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)

       if os.path.exists(filename_gf_k_0):

           with open(filename_gf_k_0, 'rb') as infile:
             green_function_k_0 = pickle.load(infile)

           filename_gf_k_variation_0 = '%s/Retarded_green_function_momentum_space_standard_deviation_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.pkl' %(Text_dir_gf_k,N,U,Mu,dtau,L)
           with open(filename_gf_k_variation_0, 'rb') as infile:
             green_function_k_std_0 = pickle.load(infile)


           filename_fit = "%s/Spline_fit_G_tau_k_point_%s_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_0.dat"%(k_dir_fit,str(k_lab[k]),N,U,Mu,str(dtau),str(L))
           g_k_data = np.copy(green_function_k_0[:,k])
           g_k_std_data = np.copy(green_function_k_std_0[:,k])

           g_k_data[-1] = -1*(1+g_k_data[0])
           g_k_std_data[-1] = g_k_std_data[0]

           omega_n, g_mf_re, g_mf_im, g_mf_re_std, g_mf_im_std = scipy_spline_fourier_transform(float(dtau),int(L),Tau,g_k_data,g_k_std_data,n_omeg,filename_fit)
           Green_function_k_mf_re = g_mf_re.copy()
           Green_function_k_mf_im = g_mf_im.copy()

           # ===================================== Dyson equation to compute self energy ========================

           self_en_real = []
           self_en_imag = []
           for n in range(len(omega_n)):
               g_mf = g_mf_re[n]+1j*g_mf_im[n]
               sig = 1j*omega_n[n]-ep_k-1/g_mf
               self_en_real.append(np.real(sig))
               self_en_imag.append(np.imag(sig))

 
           #self_en_real_hartree = []       #============= This block calculates the hartree shifted self energy ===================

           #for n in range(len(self_en_real)):
           #    hc = self_en_real[n]-self_en_real[-1]
           #    self_en_real_hartree.append(hc)

           Self_en_k_real = np.copy(np.asarray(self_en_real))
           Self_en_k_imag = np.copy(np.asarray(self_en_imag))
           #Self_en_hc_k_real = np.copy(np.asarray(self_en_real_hartree))

           # ====================================================================================================

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


                 filename_fit = "%s/Spline_fit_G_tau_k_point_%s_N_%s_U_%s_mu_%s_dtau_%s_L_%s_r_%s.dat"%(k_dir_fit,str(k_lab[k]),N,U,Mu,str(dtau),str(L),realization)
                 g_k_data = np.copy(green_function_k[:,k])
                 g_k_std_data = np.copy(green_function_k_std[:,k])

                 g_k_data[-1] = -1*(1+g_k_data[0])
                 g_k_std_data[-1] = g_k_std_data[0]

                 omega_n, g_mf_re, g_mf_im, g_mf_re_std, g_mf_im_std = scipy_spline_fourier_transform(float(dtau),int(L),Tau,g_k_data,g_k_std_data,n_omeg,filename_fit)


                 X_Gf_k_mf_re = mat_append(Green_function_k_mf_re,g_mf_re)
                 X_Gf_k_mf_im = mat_append(Green_function_k_mf_im,g_mf_im)

                 Green_function_k_mf_re = X_Gf_k_mf_re.copy()
                 Green_function_k_mf_im = X_Gf_k_mf_im.copy()

                 #======= Saving spline interpolated green function in matsuabara space for each realization =====================


                 #g_mf_data = np.stack((omega_n,g_mf_re,g_mf_im),axis=1)
                 #filename_g_mf = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_point_%s_n_max_%s_r_%s.dat'%(k_dir_mf,N,U,Mu,dtau,L,str(k),str(n_omeg),realization)
                 #np.savetxt(filename_g_mf,g_mf_data)

                 # ===================================== Dyson equation to compute self energy ========================
 
                 self_en_real = []
                 self_en_imag = []
                 
                 for n in range(len(omega_n)):
                     g_mf = g_mf_re[n]+1j*g_mf_im[n]
                     sig = 1j*omega_n[n]-ep_k-1/g_mf
                     self_en_real.append(np.real(sig))
                     self_en_imag.append(np.imag(sig))

                 # ==================================== Calculates the Hartree shifted form ============================
                 self_en_real_hartree = []

                 #for n in range(len(self_en_real)):
                 #    hc = self_en_real[n]-self_en_real[-1]
                 #    self_en_real_hartree.append(hc)

                 x_self_en_k_real = np.copy(np.asarray(self_en_real))
                 x_self_en_k_imag = np.copy(np.asarray(self_en_imag))
                 #x_self_en_hc_k_real = np.copy(np.asarray(self_en_real_hartree))

                 print("shape_ re sig",x_self_en_k_real.shape)
                 #print("shape_ re sig hartree",x_self_en_hc_k_real.shape)
                 X_Self_en_k_real = mat_append(Self_en_k_real,x_self_en_k_real) 
                 X_Self_en_k_imag = mat_append(Self_en_k_imag,x_self_en_k_imag)
                 #X_Self_en_hc_k_real = mat_append(Self_en_hc_k_real,x_self_en_hc_k_real)

                 Self_en_k_real = np.copy(X_Self_en_k_real)
                 Self_en_k_imag = np.copy(X_Self_en_k_imag)
                 #Self_en_hc_k_real = np.copy(X_Self_en_hc_k_real)

               # ====================================================================================================


       ################ Performing mean estimation by standard mean, variance ======================================

       #===== Stat of G(k) is skipped, since _v1 func already does it

       #Green_function_mf_re_k_avg = np.mean(Green_function_k_mf_re,axis = 1)
       #Green_function_mf_im_k_avg = np.mean(Green_function_k_mf_im,axis = 1)

       #Green_function_mf_re_k_std_avg = np.std(Green_function_k_mf_re,axis = 1)
       #Green_function_mf_im_k_std_avg = np.std(Green_function_k_mf_im,axis = 1)

       Self_en_k_real_avg = np.mean(Self_en_k_real,axis = 1)
       Self_en_k_imag_avg = np.mean(Self_en_k_imag,axis = 1)

       Self_en_k_real_std = (1/np.sqrt(len(Self_en_k_real[0,:])))*np.std(Self_en_k_real,axis = 1,ddof=1)
       Self_en_k_imag_std = (1/np.sqrt(len(Self_en_k_real[0,:])))*np.std(Self_en_k_imag,axis = 1,ddof=1)

       Self_en_hc_k_real_avg = np.zeros(len(omega_n))
       Self_en_hc_k_real_std = np.zeros(len(omega_n))

       for n in range(len(omega_n)):
           Self_en_hc_k_real_avg[n] = Self_en_k_real_avg[n]-Self_en_k_real_avg[-1]   #np.mean(Self_en_hc_k_real,axis = 1)
           Self_en_hc_k_real_std[n] = Self_en_k_real_std[n]                          #(1/np.sqrt(len(Self_en_k_real[0,:])))*np.std(Self_en_hc_k_real,axis = 1,ddof=1)

       #Self_en_hc_k_real
       #=======================================Saving data for analytic continuation ==========================================================


       #Gf_data = np.stack((omega_n,Green_function_mf_re_k_avg,Green_function_mf_re_k_std_avg,Green_function_mf_im_k_avg,Green_function_mf_im_k_std_avg),axis = 1)
       #filename_g_k_point = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       #np.savetxt(filename_g_k_point,Gf_data)


       #Self_en_k_real_hartree = np.zeros(len(omega_n))
       #Self_en_k_imag_hartree = np.zeros(len(omega_n))

       #for n in range(len(omega_n)):
       #    Self_en_k_real_hartree[n] = Self_en_k_real_avg[n]-0.5*float(U)*(n_den-1)
       #    Self_en_k_imag_hartree[n] = Self_en_k_imag_avg[n]



       # Checks if the standard deviation is too small for the self energies, this will cause problems in analytic continuation

       for n in range(len(Self_en_k_real_avg)):

           if Self_en_k_real_std[n] < 1e-6:
              Self_en_k_real_std[n] = 1e-6
           if Self_en_hc_k_real_std[n] < 1e-6:
              Self_en_hc_k_real_std[n] = 1e-6
           if Self_en_k_imag_std[n] < 1e-6:
              Self_en_k_imag_std[n] = 1e-6 

       Self_en_data = np.stack((omega_n,Self_en_k_real_avg,Self_en_k_real_std,Self_en_k_imag_avg,Self_en_k_imag_std),axis = 1)
       filename_self_en_k_point = '%s/Spline_interpolated_self_energy_momentum_space_avg_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat' %(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_k_point,Self_en_data)

       Self_en_hartree_data = np.stack((omega_n,Self_en_hc_k_real_avg,Self_en_hc_k_real_std,Self_en_k_imag_avg,Self_en_k_imag_std),axis = 1)
       filename_self_en_hartree_k_point = '%s/Spline_interpolated_hartree_corrected_self_energy_momentum_space_avg_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat' %(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_hartree_k_point,Self_en_hartree_data)


       ############################# Performing sample statistics estimation by jackknife sampling ============================================

       #######================================ Perform jackknife statistics of sample, save data =======================

       #===== Jacknife stat of G(k) is skipped, since _v1 func already does it
       #G_mf_re_jk = np.zeros(len(omega_n))
       #G_mf_re_jk_std = np.zeros(len(omega_n))
       #G_mf_im_jk = np.zeros(len(omega_n))
       #G_mf_im_jk_std = np.zeros(len(omega_n))

       #for n in range(len(omega_n)):

       #    g_mf_re_data = Green_function_k_mf_re[n,:]
       #    g_mf_im_data = Green_function_k_mf_im[n,:]
       #    G_mf_re_jk[n],G_mf_re_jk_std[n] = jackknife_sampling(g_mf_re_data)
       #    G_mf_im_jk[n],G_mf_im_jk_std[n] = jackknife_sampling(g_mf_im_data)

       #Gf_jk_data = np.stack((omega_n,G_mf_re_jk,G_mf_re_jk_std,G_mf_im_jk,G_mf_im_jk_std),axis = 1)
       #filename_g_k_point_jk = '%s/Spline_matsubara_frequency_retarded_green_function_momentum_space_jackknife_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat' %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       #np.savetxt(filename_g_k_point_jk,Gf_jk_data)


       #######================================ Perform jackknife statistics of sample, save data =======================

       Sig_mf_re_jk = np.zeros(len(omega_n))
       Sig_mf_re_jk_hc = np.zeros(len(omega_n))
       Sig_mf_re_jk_std = np.zeros(len(omega_n))
       Sig_mf_im_jk = np.zeros(len(omega_n))
       Sig_mf_im_jk_std = np.zeros(len(omega_n))

       for n in range(len(omega_n)):

           sig_mf_re_data = Self_en_k_real[n,:]
           sig_mf_im_data = Self_en_k_imag[n,:]
           Sig_mf_re_jk[n],Sig_mf_re_jk_std[n] = jackknife_sampling(sig_mf_re_data)
           Sig_mf_im_jk[n],Sig_mf_im_jk_std[n] = jackknife_sampling(sig_mf_im_data)

       for n in range(len(omega_n)):
           Sig_mf_re_jk_hc[n] = Sig_mf_re_jk[n]-Sig_mf_re_jk[-1]

       Sig_mf_re_jk_hc[-1] = 1e-6


       Self_en_data_jk = np.stack((omega_n,Sig_mf_re_jk,Sig_mf_re_jk_std,Sig_mf_im_jk,Sig_mf_im_jk_std),axis = 1)
       filename_self_en_k_point_jk = '%s/Spline_interpolated_self_energy_momentum_space_jackknife_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_k_point_jk,Self_en_data_jk)

       Self_en_hartree_data_jk = np.stack((omega_n,Sig_mf_re_jk_hc,Sig_mf_re_jk_std,Sig_mf_im_jk,Sig_mf_im_jk_std),axis = 1)
       filename_self_en_hartree_k_point_jk = '%s/Spline_interpolated_hartree_corrected_self_energy_momentum_space_jackknife_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(k),str(n_omeg))
       np.savetxt(filename_self_en_hartree_k_point_jk,Self_en_hartree_data_jk)


   np.savetxt("%s/Default_model_mean_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L),Mean_ac)
   np.savetxt("%s/Default_model_variance_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat" %(Text_dir_gf_k_mf_avg_ac,N,U,Mu,dtau,L),Variance_ac)
   np.savetxt("%s/Self_energy_norm_N_%s_U_%s_mu_%s_dtau_%s_L_%s.dat"%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L),Norm_self_energy)






#============================================ This function evaluates the spline interpolation of the already averaged data across several realizations========================================
#
#============================================ Also computes the self energy. Compared to above function, there is no variance in G(i\omega_n) and Sigma(i\omega_n), assigned to be a default number = 1/sqrt(measurement_sweep = 15000)====================


def green_function_self_energy_spline_FT_matsubara_calc_v2(Text_dir_gf_k_avg_ac,Text_dir_gf_k_mf_avg_ac,Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,num_den,n_omeg):


       beta = float(dtau)*float(L) 
       bz,k_points = k_point_grid_upper_half_bz(N)
       k_label = bz[:,0]
       Kx_grid = bz[:,1]
       Ky_grid = bz[:,2]


       Text_dir_spline_fit = "%s/Cubic_spline_fit"%Text_dir_gf_k_avg_ac
       if not os.path.exists(Text_dir_spline_fit):
          os.makedirs(Text_dir_spline_fit)

       for kk in range(k_points):

           kx = (2*np.pi/int(N))*Kx_grid[kk]
           ky = (2*np.pi/int(N))*Ky_grid[kk]
       
       
           #k_lab[k] = int(k_label[kk])
           #k_point = str(int(k_label[kk]))

           filename_g_k_point = '%s/Retarded_green_function_momentum_space_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat'%(Text_dir_gf_k_avg_ac,N,U,Mu,dtau,L,str(kk))
           Tau, G_k_data, G_k_std_data = np.loadtxt(filename_g_k_point,unpack = 'True', usecols = [0,1,2])
           G_k_data[-1] = -1*(1+G_k_data[0])
           G_k_std_data[-1] = G_k_std_data[0]

           #============================== Calculating G(iw_n) from Fourier transforming G(tau)===================================

           filename_fit = "%s/Spline_fit_G_tau_avg_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_point_%s.dat"%(Text_dir_spline_fit,N,U,Mu,str(dtau),str(L),str(kk))
           Omega_n, Re_G_mf, Im_G_mf, Re_G_mf_std, Im_G_mf_std = scipy_spline_fourier_transform(float(dtau),int(L),Tau,G_k_data,G_k_std_data,n_omeg,filename_fit)


           filename_G_mf_spl = "%s/Retarded_green_function_momentum_space_avg_spline_fourier_transform_to_matsubara_frequency_n_max_%s_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s.dat"%(Text_dir_gf_k_mf_avg_ac,n_omeg,N,U,Mu,dtau,L,str(kk)) 
           data_gf_mf_spl = np.stack((Omega_n,Re_G_mf,np.abs(Re_G_mf_std),Im_G_mf,np.abs(Im_G_mf_std)),axis = 1)
           np.savetxt(filename_G_mf_spl,data_gf_mf_spl)


           #==================================== Calculating the self energy from Fourier transformed G(i\omega_n) data ============================================


           ep_k = -2*(np.cos(kx)+np.cos(ky))-float(Mu)

           Sigma_k_real = np.zeros(len(Omega_n))
           Sigma_k_imag = np.zeros(len(Omega_n))

           Sigma_k_real_hartree = np.zeros(len(Omega_n))
           Sigma_k_imag_hartree = np.zeros(len(Omega_n))

           #Sigma_k_real_hartree_nrm = np.zeros(len(Omega_n))
           #Sigma_k_imag_hartree_nrm = np.zeros(len(Omega_n))

           Sigma_k_real_std = np.zeros(len(Omega_n))
           Sigma_k_imag_std = np.zeros(len(Omega_n))

           Sigma_k_real_hartree_std = np.zeros(len(Omega_n))
           Sigma_k_imag_hartree_std = np.zeros(len(Omega_n))


           for n in range(len(Omega_n)):

               G_k_mf = Re_G_mf[n]+1j*Im_G_mf[n]
               sig = 1j*Omega_n[n]-ep_k-1/G_k_mf

               Sigma_k_real[n] = np.real(sig) 
               Sigma_k_imag[n] = np.imag(sig)


               Sigma_k_real_std[n] = np.sqrt(1/15000) 
               Sigma_k_imag_std[n] = np.sqrt(1/15000) 


           for n in range(len(Omega_n)):

               Sigma_k_real_hartree[n] = Sigma_k_real[n]-0.5*float(U)*(num_den-1)
               Sigma_k_imag_hartree[n] = Sigma_k_imag[n]
               Sigma_k_real_hartree_std[n] = np.sqrt(1/15000)                # this is just a proxy error estimate, 15000 = number of samples
               Sigma_k_imag_hartree_std[n] = np.sqrt(1/15000)

               #Sigma_k_real_hartree_nrm[n] = (Sigma_k_real[n]-Sig0)/Sig1
               #Sigma_k_imag_hartree_nrm[n] = (Sigma_k_imag[n])/Sig1
               #Sigma_k_real_nrm_std[n] = Sigma_k_real_std[n]/Sig1
               #Sigma_k_imag_nrm_std[n] = Sigma_k_imag_std[n]/Sig1

           Sigma_data = np.stack((Omega_n,Sigma_k_real,Sigma_k_real_std,Sigma_k_imag,Sigma_k_imag_std),axis = 1)
           filename_sigma_k_point = '%s/Spline_interpolation_avg_data_Self_energy_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk),str(n_omeg))
           np.savetxt(filename_sigma_k_point,Sigma_data)

           Sigma_data_hartree = np.stack((Omega_n,Sigma_k_real_hartree,Sigma_k_real_std,Sigma_k_imag_hartree,Sigma_k_imag_std),axis =1)
           filename_sigma_k_point_hartree = '%s/Spline_interpolation_avg_data_Self_energy_hartree_corrected_matsubara_frequency_N_%s_U_%s_mu_%s_dtau_%s_L_%s_k_label_%s_n_max_%s.dat'%(Text_dir_self_en_k_mf_avg_ac,N,U,Mu,dtau,L,str(kk),str(n_omeg))
           np.savetxt(filename_sigma_k_point_hartree,Sigma_data_hartree)



def main(total,cmdargs):
    if(total!=8):
        raise ValueError('missing args')

    N = cmdargs[1]
    U = cmdargs[2]
    mu = cmdargs[3]
    L = cmdargs[4]
    Dtau = cmdargs[5]
    run_no = int(cmdargs[6])
    N_omeg = int(cmdargs[7])

    print(U,"U")
    print(L,"L")
    Beta = str((float(L))*(float(Dtau)))
    print(Beta,"Beta")

    #==========================Diretories for real space correlation averages =============================================

    Text_dir_gf_r = '../../../Text_files_longer_sweep/Text_files_N_%s_real_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_real_space'%(N,N,U,Dtau,mu,Dtau,L)


    #==================== Directories for momentum space correlation averages ===================================================


    Text_dir_gf_k = '../../../Text_files_longer_sweep/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_gf_k_mf = '../../../Text_files_longer_sweep/Text_files_N_%s_momentum_space_correlations/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space_matsubara_frequency'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_k_mf):
       os.makedirs(Text_dir_gf_k_mf) 


     #===========================Directories for correlations for analytic continuation averages, matsubara frequency space =========================================

    Text_dir_gf_k_avg_ac = '../../../Text_files_longer_sweep/Text_files_N_%s_analytic_continuation/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space_averaged'%(N,N,U,Dtau,mu,Dtau,L)

    Text_dir_gf_k_mf_avg_ac = '../../../Text_files_longer_sweep/Text_files_N_%s_analytic_continuation/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Retarded_green_functions_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_gf_k_mf_avg_ac):
        os.makedirs(Text_dir_gf_k_mf_avg_ac)


    Text_dir_self_en_k_mf_avg_ac_v1 = '../../../Text_files_longer_sweep/Text_files_N_%s_analytic_continuation_self_energy_v1/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Self_energy_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_self_en_k_mf_avg_ac_v1):
        os.makedirs(Text_dir_self_en_k_mf_avg_ac_v1)

    Text_dir_self_en_k_mf_avg_ac_v2 = '../../../Text_files_longer_sweep/Text_files_N_%s_analytic_continuation_self_energy_v2/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Self_energy_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_self_en_k_mf_avg_ac_v2):
        os.makedirs(Text_dir_self_en_k_mf_avg_ac_v2)
        
    Text_dir_self_en_k_mf_avg_ac = '../../../Text_files_longer_sweep/Text_files_N_%s_analytic_continuation_self_energy_avg_data/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Self_energy_momentum_space_matsubara_frequency_averaged'%(N,N,U,Dtau,mu,Dtau,L)
    if not os.path.exists(Text_dir_self_en_k_mf_avg_ac):
        os.makedirs(Text_dir_self_en_k_mf_avg_ac)

    Text_dir_eqm = "../../../Text_files_longer_sweep/Text_files_N_%s/Text_files_N_%s_U_%s_dtau_%s/Mu_%s/dtau_%s_L_%s/Thermodynamic_measurements"%(N,N,U,Dtau,mu,Dtau,L)
    
    filename_eqm_avg = '%s/Thermodynamic_measurements_normal_averaged_dictionary_N_%s_U_%s_mu_%s_dtau_%s_L_%s.pkl' %(Text_dir_eqm,N,U,mu,Dtau,L)
    with open(filename_eqm_avg, 'rb') as infile:
         sys_measure_avg = pickle.load(infile)

    Nden = sys_measure_avg['Density averaged']
    print("avg sign, density, mu", Nden, mu)
    
    
    # ================================================= Spline interpolation, Matsubara frequency space data averaging   ================================================================
    
    green_function_spline_FT_matsubara_averaging_momentum_space_v1(Text_dir_gf_k,Text_dir_gf_k_mf,Text_dir_gf_k_avg_ac,Text_dir_gf_k_mf_avg_ac,Text_dir_self_en_k_mf_avg_ac_v1,N,U,mu,Dtau,L,Nden,run_no,N_omeg)
    green_function_spline_FT_matsubara_averaging_momentum_space_v2(Text_dir_gf_k,Text_dir_gf_k_mf,Text_dir_gf_k_avg_ac,Text_dir_gf_k_mf_avg_ac,Text_dir_self_en_k_mf_avg_ac_v2,N,U,mu,Dtau,L,Nden,run_no,N_omeg)

    #================================================= Spline interpolation of averaged data, Matsurbara_frequency space data calc ==================================================

    green_function_self_energy_spline_FT_matsubara_calc_v2(Text_dir_gf_k_avg_ac,Text_dir_gf_k_mf_avg_ac,Text_dir_self_en_k_mf_avg_ac,N,U,mu,Dtau,L,Nden,N_omeg)



if __name__ == '__main__':
    sys.argv
    total = len(sys.argv)
    print("No of sys arguments",total)
    cmdargs = sys.argv
    main(total,cmdargs)
