# -*- coding: utf-8 -*-

#use scipy DLS to converge Psi Ellipsometric data to model
#arl92@case.edu 20190926

#needs: data for model materials (n,k 300-5-1000nm)
#       experimental data(psi: 300-5-1000nm)
#       TMM module in working folder


import TMM as tmm
import numpy as np
import matplotlib.pyplot as plt
#from scipy.optimize import leastsq #legacy version of LMA
from scipy.optimize import least_squares

###################################################################

#import materials data from csv files, assuming [n k]
a = np.loadtxt('mats/7059g.csv',delimiter=',',skiprows=1)
n_7059 = a[:,1]
k_7059 = a[:,2]
a = np.loadtxt('mats/al2o3.csv',delimiter=',',skiprows=1)
n_al2o3 = a[:,1]
k_al2o3 = a[:,2]
a = np.loadtxt('mats/ag-g.csv',delimiter=',',skiprows=1)
n_ag = a[:,1]
k_ag = a[:,2]
a = np.loadtxt('mats/a-ge.csv',delimiter=',',skiprows=1)
wavelength = a[:,0]
n_ge = a[:,1]
k_ge = a[:,2]

#import experimental data
a = np.loadtxt('data/experimental_1lay_enz_psi.txt')
psi_e = a[:,1]

###################################################################

#model

#  --> AIR [ al2o3 ][ ag ][ ge ] GLASS

#index of refraction for layers, cover(above), and substrate(below)
n = np.array([n_al2o3+1j*k_al2o3, n_ag+1j*k_ag, n_ge+1j*k_ge])
n_cover = 1 #air
n_subst = n_7059+1j*k_7059 #corning 7059 glass slide

#user defined thickness of each layer, m 
l = np.array([33E-9,10E-9,1E-10]) #initial guess

#angle of incidence, degrees
ang_of_inc = 45

#wavelength of light incident in m
wavelength *= 1E-9 #use wavelength from data files (300-5-1000nm)

#Setup for the MC errorbars
#if you dont want this set num_mc=1,sig=0
num_mc = 1 #number of MC iterations
sig = 0 #nm,stdev of dist to search
np.random.seed(5497) #Seed the RNG
#initialization
thick_mc = np.zeros((num_mc,l.size))
mse = np.zeros(num_mc)
l_er = np.ones(l.size)

###################################################################

#DLS thickness fitting

#residuals func def for psi experimental
def residuals_fcn(l_fit,psi_e,wavelength):
    psi_t = np.zeros(wavelength.size)
    delta_t = np.zeros(wavelength.size)
    for i in range(0,wavelength.size):
        (psi_t[i],delta_t[i]) = tmm.ellips(ang_of_inc, wavelength[i], n[:,i], l_fit, n_cover, n_subst[i])
    #calc the residuals
    return psi_e-psi_t

#least squares optimization in scipy for fitting layer thicknesses
#loop for the MC errorbars
for ii in range(0,num_mc):
    #add a random number to the thickness guesses
    for el in range(0,l.size):
        l_er[el] = l[el]+np.random.normal(loc=0.0,scale=sig)
        #no negative thicknesses
        if (l_er[el]) <= 0.0:
            l_er[el] = 1E-10
    #plsq = leastsq(residuals_fcn,l,args = (psi_e,wavelength),maxfev = 2000) #legacy LMA implementation without boundaries
    #fit thicknesses with SciPy TRF
    plsq = least_squares(residuals_fcn,l_er,args = (psi_e,wavelength),bounds = (0,np.inf),method = 'trf',gtol = 1E-10,loss = 'soft_l1',max_nfev = 1000) #TRF ls optimization routine, not sure how robust
    #save thickness data and MSE
    for m in range(0,l.size):
        thick_mc[ii,m] = plsq.x[m]
    mse[ii] = sum(np.abs(residuals_fcn(thick_mc[ii,:],psi_e,wavelength)))/psi_e.size
    print('Layer %1d / %1d ...'%(ii+1,num_mc))

#find the thickness guess that minimizes the MSE, use the error on mean as uncertainty
it = np.argmin(mse)
mse_m = mse[it]
thick_m = thick_mc[it,:]
thick_err = np.std(thick_mc,axis=0)

#calculate the thickness of the whole stack (could be measured expeirmentally)
stack_thick = sum(thick_m)*1E9
stack_stderr = np.sqrt(sum(thick_err**2))*1E9/np.sqrt(num_mc)

###################################################################

#output to terminal

#calculate the best fit (by MSE) psi and delta
psi_f = np.zeros(wavelength.size)
delta_f = np.zeros(wavelength.size)
for i in range(0,wavelength.size):
    (psi_f[i],delta_f[i]) = tmm.ellips(ang_of_inc, wavelength[i], n[:,i],thick_m, n_cover, n_subst[i])
 
#Plot experimental and fit Psi for comparison
mse= sum(np.sqrt((psi_e-psi_f)**2))/psi_e.size
plt.plot(wavelength*1E9,psi_e,'k.',wavelength*1E9,psi_f,'b--')
plt.title('Psi TRF fit thickness')
plt.xlabel('Wavelength (nm)')
plt.ylabel('Angle (deg)')
plt.legend(['Experimental Data','Psi fit'])
plt.show()

#Print relevant data to terminal
print('Thicknesses (nm):',thick_m*1E9)
print('Errorbars   (nm):',thick_err*1E9)
print('Stack Thick (nm): %1.2f +/- %1.1f'%(stack_thick,stack_stderr))
print('MSE:%3.4f'%(mse_m))
