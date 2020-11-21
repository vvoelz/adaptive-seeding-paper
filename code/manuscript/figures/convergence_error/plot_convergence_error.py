# plot_convergence.py            

import os, sys
sys.path.append('../../../')

import numpy as np

from BinManager import *
from twostate import *

ntrials = 20
nseeds = [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000] 
trajlengths = [10, 100, 1000]  #, 1000]

nbins = 20
bins = BinManager(nbins=nbins)

lagtime = 100

# arrays to store the implied timescale results
## NOTE: inner dimension is tau_MLE, tau_MLE_known_pi, tau_weighted_MLE, tau_rownorm
tau = np.ones( (4, ntrials, len(trajlengths), len(nseeds)) )

# arrays to store the delta_f results (i.e. the two state $\Delta F of the left basin)
## NOTE: inner dimension is delF_MLE, delF_MLE_known_pi, delF_weighted_MLE, delF_rownorm (all in kcal/mol)
delta_f = np.ones( (4, ntrials, len(trajlengths), len(nseeds)) )

# arrays to store the PMFs
## NOTE inner dimension is pmf_MLE, pmf_MLE_known_pi, pmf_weighted_MLE, pmf_rownorm
pmf = np.ones( (4, ntrials, len(trajlengths), len(nseeds), nbins) )

actual_timescale = 9.66e6
std_actual_timescale = 1.37e6

actual_delta_f = 4.0

#### Read in all the results #####
for trial in range(ntrials):

    for seed in range(len(nseeds)):
        s = nseeds[seed]

        for itrajlength in range(len(trajlengths)):
            n = trajlengths[itrajlength]

            indir = '../../../results/seeding_tests/trajlength%d/nseeds%d'%(n, s)

            # read in timescales
            infile = os.path.join(indir,'trial%d_timescales.txt'%trial)
            tau[:, trial, itrajlength, seed] = np.loadtxt(infile)

            # read in delta_f 
            infile = os.path.join(indir,'trial%d_deltaF.txt'%trial)
            delta_f[:, trial, itrajlength, seed] = np.loadtxt(infile)

            # read in pmf 
            infile = os.path.join(indir,'trial%d_pmf.txt'%trial)
            data = np.loadtxt(infile)
            pmf[:, trial, itrajlength, seed, :] = data.transpose() 


print 'tau', tau
print 'delta_f', delta_f
print 'pmf', pmf

#### Calcluate the mean and std across all the trials

# for the timescales, average the log
unsigned_err = np.abs(tau - actual_timescale)
mean_tau_error = unsigned_err[(tau != 1.0)].mean(axis=1)
std_tau_error  = unsigned_err[(tau!=1.0)].std(axis=1)
### mean log(tau) needed for error bars
mean_log_tau_error = np.log(unsigned_err).mean(axis=1)
std_log_tau_error = np.log(unsigned_err).std(axis=1) 


mean_delta_f_error = (delta_f - actual_delta_f).mean(axis=1)
std_delta_f_error = (delta_f - actual_delta_f).std(axis=1)

mean_pmf = pmf.mean(axis=1)
std_pmf = pmf.std(axis=1)

print 'mean_tau_error', mean_tau_error, 'mean_tau_error.shape', mean_tau_error.shape
print 'std_tau_error', std_tau_error, 'std_tau_error.shape', std_tau_error.shape


# See also: http://matplotlib.org/users/customizing.html

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

fontfamily={'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)


labels = ['MLE', 'MLE known_pi', 'MLE pop-weighted', 'row-norm']
colors = ['r','b','m','y']
linewidth = 1.5
markersize = 4

plt.figure(figsize=(8,8))

for ntraj_index in range(len(trajlengths)):

    plt.subplot(len(trajlengths),2,2*ntraj_index+1)

    for i in range(4): 
        plt.plot(nseeds, mean_tau_error[i,ntraj_index,:], colors[i]+'o-', ms=markersize, lw=linewidth, label=labels[i]) 
        lower_bound = mean_tau_error[i,ntraj_index,:]/np.exp(std_log_tau_error[i,ntraj_index,:])
        upper_bound = mean_tau_error[i,ntraj_index,:]*np.exp(std_log_tau_error[i,ntraj_index,:])
        #plt.fill_between(nseeds, mean_tau[i,ntraj_index,:]-std_tau[i,ntraj_index,:], mean_tau[i,ntraj_index,:]+std_tau[i,ntraj_index,:], facecolor=colors[i], alpha=0.25)
        plt.fill_between(nseeds, lower_bound, upper_bound, facecolor=colors[i], alpha=0.25)

    #plt.plot(nseeds, [actual_timescale]*len(nseeds), 'k-', lw=linewidth, label='actual')
    #plt.fill_between(nseeds, [actual_timescale - std_actual_timescale]*len(nseeds), [actual_timescale + std_actual_timescale]*len(nseeds), facecolor='k', alpha=0.25)
    
    plt.xscale('log')
    #plt.xticks(trajlengths)
    plt.xlabel('numbers of trajectories')
    plt.xlim(10,1000)

    plt.yscale('log')
    #plt.ylim(-5,10)
    plt.ylabel('mean unsigned error in $\tau$')
    plt.ylim(1e-5,1e11)

    plt.title('$10^{%d}$ steps'%np.log10(trajlengths[ntraj_index]*lagtime) )
    plt.legend(loc='best', fontsize=9)
    
for ntraj_index in range(len(trajlengths)):

    plt.subplot(len(trajlengths),2,2*ntraj_index+2)

    for i in range(4):
        plt.plot(nseeds, mean_delta_f_error[i,ntraj_index,:], colors[i]+'o-', ms=markersize, lw=linewidth, label=labels[i])
        plt.fill_between(nseeds, mean_delta_f_error[i,ntraj_index,:]-std_delta_f_error[i,ntraj_index,:], mean_delta_f_error[i,ntraj_index,:]+std_delta_f_error[i,ntraj_index,:], facecolor=colors[i], alpha=0.25)

    plt.plot(nseeds, [0.0]*len(nseeds), 'k-', lw=linewidth)  #, label='actual')

    plt.xscale('log')
    plt.xlim(10,1000)
    #plt.xticks(trajlengths)
    plt.xlabel('numbers of trajectories')

    #plt.yscale('log')
    plt.ylabel('error in $\Delta F$ estimate (kcal/mol)')
    plt.ylim(-5,5)

    plt.title('$10^{%d}$ steps'%np.log10(trajlengths[ntraj_index]*lagtime) )
    plt.legend(loc='best', fontsize=9)

plt.tight_layout()    
plt.savefig('convergence_error.pdf')


