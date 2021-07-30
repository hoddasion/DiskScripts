# -*- coding: utf-8 -*-
"""
Created on Fri Nov  1 17:39:10 2019

@author: kse18nru
"""



import numpy as np
import matplotlib.pyplot as plt

#%%
ks = np.linspace(10**(-2),1,1000)
print(ks)
ks_log = np.log(ks)
k3 = ks**(-3)
k3_log = np.log(k3)
print('y-intercept',k3[0])

#fig, ax = plt.subplots(1,1, figsize = (10,15))

#q1 = ax.plot(ks_log, k3)


#%%
x = 15
if x % 3 == 0:
    print('\n\n\n\nDivisor 3')
elif x > 0:
    print('\n\n\n\nis positive')
else:
    print('\n\n\n\nNone of the above')