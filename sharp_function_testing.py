# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 12:59:10 2022

@author: kse18nru
"""
#%%
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams.update({'font.size': 24})
#%% define sharp functions
def calc_A(Rit, g0 = 10):
    return (1 - g0*Rit) / (1 - g0*Rit/2)**2

def calc_B(Rit, g0 = 10):
    return (g0/2) / (1 - g0*Rit/2)**2

def sharp_stability_function(Ri, Rit, g0 = 10, return_concatenated = True):
    """
    The family of sharp stabiliy functions, with SHARPEST transition richardson number.
    Only applies to Ri values above 0

    Parameters
    ----------
    Ri : must be 1D nparray of rischardson numbers
    Rit : integer, transition richardson number
    g0 :  The default is 10.

    Returns
    -------
    Two nparray f stability function values: first array is values below or equal to Rit, second is for above

    """
    f_below_threshold = (1 - 5*Ri[Ri <= Rit])**2
    f_above_threshold = (calc_A(Rit, g0=g0) + calc_B(Rit, g0=g0)*Ri[Ri > Rit])**(-2)
    if return_concatenated:
        f_full = np.concatenate((f_below_threshold, f_above_threshold))
        return f_full
    else:    
        return f_below_threshold, f_above_threshold
    
def long_tail(Ri):
    return (1+10*Ri)**(-1)

#%%


Ri_values = np.arange(0,41)/100


f_rit0p1_g10 = sharp_stability_function(Ri_values, 0.1)
f_rit0p15_g10 = sharp_stability_function(Ri_values, 0.15)
f_rit0p05_g10 = sharp_stability_function(Ri_values, 0.05)
f_rit0p01_g10 = sharp_stability_function(Ri_values, 0.01)

f_rit0p1_g5 = sharp_stability_function(Ri_values, 0.1, g0 = 5)
f_rit0p15_g5 = sharp_stability_function(Ri_values, 0.15, g0 = 5)
f_rit0p05_g5 = sharp_stability_function(Ri_values, 0.05, g0 = 5)
f_rit0p01_g5 = sharp_stability_function(Ri_values, 0.01, g0 = 5)

long_tail_values = long_tail(Ri_values)

#%% plot

fig, ax = plt.subplots(1,1, figsize = (10,16))

ax.plot(Ri_values, f_rit0p15_g10, label = r'$Ri_t = 0.15$')
ax.plot(Ri_values, f_rit0p1_g10, label = r'$Ri_t = 0.1$')
ax.plot(Ri_values, f_rit0p05_g10, label = r'$Ri_t = 0.05$')
ax.plot(Ri_values, f_rit0p01_g10, label = r'$Ri_t = 0.01$')
ax.plot(Ri_values, long_tail_values, label = 'Long tail ')
ax.grid()
ax.set_xlabel(r'$Ri$')
ax.set_ylabel(r'$f_{stable}$')

ax.set_xticks([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
ax.set_xticklabels([0.0,'',0.1,'',0.2,'',0.3,'',0.4])

ax.set_yticks([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])
ax.set_yticklabels([0.0,'',0.1,'',0.2,'',0.3,'',0.4,'',0.5,'',0.6,'',0.7,'',0.8,'',0.9,'',1.0])
ax.set_title(r'SHARPEST functions $g_0$ = 10')
plt.legend()

# =============================================================================
# ax1.plot(Ri_values, f_rit0p15_g5, label = r'$Ri_t = 0.15$')
# ax1.plot(Ri_values, f_rit0p1_g5, label = r'$Ri_t = 0.1$')
# ax1.plot(Ri_values, f_rit0p05_g5, label = r'$Ri_t = 0.05$')
# ax1.plot(Ri_values, f_rit0p01_g5, label = r'$Ri_t = 0.01$')
# ax1.grid()
# ax1.set_xlabel(r'$Ri$')
# #ax1.set_ylabel(r'$f_{stable}$')
# 
# ax1.set_xticks([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
# ax1.set_xticklabels([0.0,'',0.1,'',0.2,'',0.3,'',0.4])
# 
# ax1.set_yticks([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])
# ax1.set_yticklabels([0.0,'',0.1,'',0.2,'',0.3,'',0.4,'',0.5,'',0.6,'',0.7,'',0.8,'',0.9,'',1.0])
# ax1.set_title(r'$g_0$ = 5')
# =============================================================================

#fig.suptitle('SHARPEST functions')

#plt.savefig('D:/Project/Figures/PNG/Theory/SHARPEST_functions_fvsRi.png')

#%% simplified for poster


  #%% plot

fig, ax = plt.subplots(1,1, figsize = (10,16))

ax.plot(Ri_values, f_rit0p15_g10, label = 'SHAPREST', linewidth = 3)
#ax.plot(Ri_values, f_rit0p1_g10, label = r'$Ri_t = 0.1$')
#ax.plot(Ri_values, f_rit0p05_g10, label = r'$Ri_t = 0.05$')
#ax.plot(Ri_values, f_rit0p01_g10, label = r'$Ri_t = 0.01$')
ax.plot(Ri_values, long_tail_values, label = 'Long-tail', linewidth = 3)
#ax.grid()
ax.set_xlabel(r'$Ri$')
ax.set_ylabel(r'$f_{stable}$')

ax.set_xticks([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4])
ax.set_xticklabels([0.0,'',0.1,'',0.2,'',0.3,'',0.4])

#ax.set_yticks([0.0,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7,0.75,0.8,0.85,0.9,0.95,1])
#ax.set_yticklabels([0.0,'',0.1,'',0.2,'',0.3,'',0.4,'',0.5,'',0.6,'',0.7,'',0.8,'',0.9,'',1.0])

plt.legend()  