# -*- coding: utf-8 -*-
#import sys
#sys.path.append('D:/Desktop/tmm_program/mods')
import TMM as tmm
import numpy as np
import tables
import matplotlib.pyplot as plt


a = np.loadtxt('glass_n_k_M2k.txt',skiprows=1)
n_7059 = a[:,1]
k_7059 = a[:,2]
a = np.loadtxt('ITO_n_k.txt',skiprows=1)
n_ito = a[:,1]
k_ito = a[:,2]
a = np.loadtxt('Ag_n_k.txt',skiprows=1)
n_ag = a[:,1]
k_ag = a[:,2]
wavelength = a[:,0]



#################################################
#user defined index of refraction for layers, cover(above), and substrate(below)
n = np.array([n_ito+1j*k_ito, n_ag+1j*k_ag, n_ito+1j*k_ito, n_ag+1j*k_ag])
#n = np.array([n_al2o3+1j*k_al2o3, n_ag+1j*k_ag, n_ge+1j*k_ge,n_al2o3+1j*k_al2o3, n_ag+1j*k_ag, n_ge+1j*k_ge,n_al2o3+1j*k_al2o3, n_ag+1j*k_ag, n_ge+1j*k_ge,n_al2o3+1j*k_al2o3, n_ag+1j*k_ag, n_ge+1j*k_ge,n_al2o3+1j*k_al2o3, n_ag+1j*k_ag, n_ge+1j*k_ge]) 
n_cover = 1 #ambient
n_subst = n_7059+1j*k_7059 #corning 7059 glass slide
#user defined thickness of each layer, m
l = np.array([1E-8,1E-8,30E-8,1E-8])
#l = np.array([39E-9,7E-9,1E-9,43E-9,7E-9,1E-9,38E-9,7E-9,1E-9,34E-9,7E-9,1E-9,43E-9,7E-9,1E-9])
#angle of incidence, degrees
ang_of_inc = 50
# wavelength of light incident in m
wavelength *= 1E-9

#parameters for generating random data set
sig1 = 2E-8
sig2 = 20E-8
np.random.seed(839137)
set_length = 60000

#filesave for dataset
filename = 'traindata_itocav_60000n_DBW.h5'
f = tables.open_file(filename, mode='w')
atom = tables.Float64Atom()
array_c = f.create_earray(f.root, 'data',atom, (0,l.size + 2*wavelength.size))
#################################################
psi = np.zeros(wavelength.size)
delta = np.zeros(wavelength.size)


l_no = np.zeros(l.size)
for g in range(0,set_length):
    for el in range(0,l.size):
        l_no[el] = l[el]+np.random.normal(loc=0.0,scale=sig1)
        if el == 2:
            l_no[el] = l[el]+np.random.normal(loc=0.0,scale=sig2)
        if (l_no[el]) < 1E-10:
            l_no[el] = 1E-10
    for i in range(0,wavelength.size):
        (psi[i],delta[i]) = tmm.ellips(ang_of_inc, wavelength[i], n[:,i], l_no, n_cover, n_subst[i])
   
    array_c.append(np.reshape(np.concatenate((l_no,psi, delta)),(1,l.size + 2*wavelength.size)))
    if (g%1000==1):
       print('Processing System %1d/%1d ...'%(g,set_length))
      
f.close()
print('Completed!')

