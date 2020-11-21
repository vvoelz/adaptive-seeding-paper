import os, sys
import  numpy as np

# See also: http://matplotlib.org/users/customizing.html

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

fontfamily={'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)


fdata = np.loadtxt('long_traj_DTRAM.free_energies')   # in kcal/mol

"""
#k  bootstrap f(kcal/mol) ... 
0   0    2.35503014  2.11610359  2.04031726  2.12146301  2.35695322  2.74649161 3.28983793  3.97657149  4.69702787  4.62023371  3.80035372  2.81729098 1.95185535  1.25115304  0.70521978  0.31346704  0.07856213  0. 0.07845942  0.3136316
1   0    2.81594434  2.53132883  2.44031723  2.53668797  2.8178647   3.2835708 3.93344241  4.75532526  5.62041096  5.55079617  4.54440571  3.36626214 2.33307481  1.49501036  0.84232504  0.37438119  0.09378737  0. 0.09368466  0.3745458 
2   0    3.27631529  2.94641618  2.8403172   2.95177505  3.27823304  3.81945968 4.57502257  5.53128847  6.54211422  6.4796022   5.28368919  3.91109925  2.71123532  1.73682056  0.97823825  0.4347521   0.10
...
"""

from BinManager import *
nbins = 20
bins = BinManager(nbins=nbins)

lam_values = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
colors = ['m','y','c','r','g','b','k']
K = len(lam_values)

plt.figure( figsize=(3.3, 3) )
xmin, xmax = 1.5, 5.5
for k in range(K):

    kslice = fdata[ (fdata[:,0]==k), 2:]
    #print 'kslice', kslice
    
    f = kslice.mean(axis=0)
    df = kslice.std(axis=0)
    
    print 'k =', k, df

    plt.plot(bins.centers,f, label='$\lambda$ = %2.1f'%lam_values[k], color=colors[k], linewidth=1)
    plt.fill_between(bins.centers, f,  f-df, f+df, facecolor=colors[k], alpha=0.25)

    plt.xlabel('$x$')
    plt.ylabel('$U_{\lambda}(x)$ (kcal/mol)')
    plt.xlim(xmin, xmax)
    plt.legend(loc='best', fontsize=8)
    plt.tight_layout()

plt.savefig('dtram_free_energies.pdf')

##### timescales ####
tau_data  = np.loadtxt('long_traj_DTRAM.timescales')[:,2]   # in kcal/mol
tau_mean = tau_data.mean(axis=0)
tau_std = tau_data.std(axis=0)

print 'tau =', tau_mean, '+/-', tau_std, 'steps.'

plt.figure( figsize=(3.3, 3) )
counts, bin_edges = np.histogram(tau_data)
bin_centers = (bin_edges[0:-1] + bin_edges[1:])/2.0
plt.bar(bin_centers, counts)
plt.title( 'tau = %3.2e +/- %3.2e steps'%(tau_mean,tau_std), fontsize=10)
plt.tight_layout()
plt.savefig('dtram_timescales.pdf')



