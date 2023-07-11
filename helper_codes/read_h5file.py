# -*- coding: utf-8 -*-
"""
Created on Thu Sep 26 21:24:28 2019

@author: Andy
"""
import tables
import numpy as np
import matplotlib.pyplot as plt

f = tables.open_file('../data/1layenz_CNN_data/traindata_1layenz_60000n_3923s_20190926.h5',mode='r')
a = (f.root.data[50000:50500,:])
f.close()

plt.plot(np.arange(295,1000,5),a[1,3:144])
plt.show()