# plot_convergence.py            

import os, sys
sys.path.append('../../../')

import numpy as np

from BinManager import *
from twostate import *

ntrials = 20
nseeds = [5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200, 300, 500, 750, 1000] 
#nseeds = [5,10,20,50]
trajlengths = [10, 100, 1000]  #, 1000]

nbins = 20
bins = BinManager(nbins=nbins)

lagtime = 100

# arrays to store the implied timescale results
## NOTE: inner dimension is tau_MLE, tau_MLE_known_pi, tau_weighted_MLE, tau_rownorm, tau_corrected, tau_OOM, rankM
tau = np.ones( (7, ntrials, len(trajlengths), len(nseeds)) )

# arrays to store the delta_f results (i.e. the two state $\Delta F of the left basin)
## NOTE: inner dimension is delF_MLE, delF_MLE_known_pi, delF_weighted_MLE, delF_rownorm, delF_corrected (all in kcal/mol)
delta_f = np.ones( (5, ntrials, len(trajlengths), len(nseeds)) )

# arrays to store the PMFs
## NOTE inner dimension is pmf_MLE, pmf_MLE_known_pi, pmf_weighted_MLE, pmf_rownorm, pmf
pmf = np.ones( (5, ntrials, len(trajlengths), len(nseeds), nbins) )

#actual_timescale = 9.66e6
actual_timescale = 5.0e6   # Jun 6 estimate
std_actual_timescale = 1.37e6


#### Read in all the results #####
for trial in range(ntrials):

    for seed in range(len(nseeds)):
        s = nseeds[seed]

        for itrajlength in range(len(trajlengths)):
            n = trajlengths[itrajlength]

            indir = '../../../results-oct2018/seeding_tests/trajlength%d/nseeds%d'%(n, s)

            # read in timescales
            infile = os.path.join(indir,'trial%d_timescales.txt'%trial)
            tau[:, trial, itrajlength, seed] = np.loadtxt(infile)
            #print trial, itrajlength, seed, 'rankM', tau[6, trial, itrajlength, seed]

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
mean_tau = tau.mean(axis=1)
std_tau = tau.std(axis=1)
### mean log(tau) needed for error bars
mean_log_tau = np.log(tau).mean(axis=1)
std_log_tau = np.log(tau).std(axis=1) 


mean_delta_f = delta_f.mean(axis=1)
std_delta_f = delta_f.std(axis=1)

mean_pmf = pmf.mean(axis=1)
std_pmf = pmf.std(axis=1)

print 'mean_tau', mean_tau, 'mean_tau.shape', mean_tau.shape
print 'std_tau', std_tau, 'std_tau.shape', std_tau.shape


# See also: http://matplotlib.org/users/customizing.html

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

fontfamily={'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)


labels = ['MLE', 'MLE known_pi', 'MLE pop-weighted', 'row-norm', 'corrected', 'OOM', 'rank M']
colors = ['r','b','m','y','c','g','g']  # last one for rankM 
linewidth = 1.5
markersize = 4

UseLegend = False

plt.figure(figsize=(8,8))

for ntraj_index in range(len(trajlengths)):

    plt.subplot(len(trajlengths),2,2*ntraj_index+1)

    #for i in range(len(labels)): 
    for i in [6]:
        plt.plot(nseeds, mean_tau[i,ntraj_index,:], colors[i]+'o-', ms=markersize, lw=linewidth, label=labels[i]) 
        lower_bound = mean_tau[i,ntraj_index,:]/np.exp(std_log_tau[i,ntraj_index,:])
        upper_bound = mean_tau[i,ntraj_index,:]*np.exp(std_log_tau[i,ntraj_index,:])
        #plt.fill_between(nseeds, mean_tau[i,ntraj_index,:]-std_tau[i,ntraj_index,:], mean_tau[i,ntraj_index,:]+std_tau[i,ntraj_index,:], facecolor=colors[i], alpha=0.25)
        plt.fill_between(nseeds, lower_bound, upper_bound, facecolor=colors[i], alpha=0.25)

    #plt.plot(nseeds, [actual_timescale]*len(nseeds), 'k-', lw=linewidth, label='actual')
    #plt.fill_between(nseeds, [actual_timescale - std_actual_timescale]*len(nseeds), [actual_timescale + std_actual_timescale]*len(nseeds), facecolor='k', alpha=0.25)

    plt.xscale('log')
    #plt.xticks(trajlengths)
    plt.xlabel('numbers of trajectories')
    plt.xlim(5,1000)

    #plt.yscale('log')
    plt.ylim(0,12)
    plt.ylabel('implied timescale estimate')
    #plt.ylim(1e3,1e11)

    plt.title('$10^{%d}$ steps'%np.log10(trajlengths[ntraj_index]*lagtime) )
    #if (ntraj_index==2):
    #    plt.legend(loc='best', fontsize=7)
    

plt.tight_layout()    
plt.savefig('rankM_oct2018.pdf')


