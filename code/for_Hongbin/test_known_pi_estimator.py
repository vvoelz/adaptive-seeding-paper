import os, sys
sys.path.append('../')

import numpy as np

from twostate import *
from estimators import *
from BinManager import *


# Here's a test of the MLE_known_pi estimator

# First, let's read in one of the count matrices that I built
# from your long MC runs on cb2rr
C = np.load('counts_20bins_229.npy')
lagtime = 100   # steps

# The count matrix should be a 20 x 20 sqaure matrix
assert C.shape[0] == 20
assert C.shape[0] == C.shape[1] 

# I have a helper object called a BinManager that just keeps track of the bin sizes, edges, and bin counts
nbins = 20
bins = BinManager(nbins=nbins)

### Calculate the known pi of each bin

# Here we need to be a bit careful, because the true populations arent just the the
# Boltzmann weight at the bin centers, but rather the integrated area between the bin edges.
# I will numerically approxximate this use a very fine-graining binning and summing up the
# Boltzmann weights

fine_grain = 1000
nbins_fine = fine_grain*nbins
bins_fine = BinManager(nbins=nbins_fine)
pi_fine = np.exp(-U(bins_fine.centers)/0.596)

# sum the fine-grain populations to get the coarse-grained pops (nbins=20) 
pi = []
for i in range(nbins):
    pi.append(pi_fine[i*fine_grain:(i+1)*fine_grain].sum())
pi = np.array(pi)
pi = pi/pi.sum()

# Use the MLE_tProb_known_pi() estimator to get tProb
tProb_known_pi = MLE_tProb_known_pi(C, pi, verbose=True)

# calculate equil populations from tProb -- this should match up with the known result!
pi_known_pi = get_stationary_pops(tProb_known_pi)
f_known_pi = -0.596*np.log(pi)  # free energies in kcal/mol
f_known_pi = f_known_pi - f_known_pi.min() # make lowest free energy value 0
print '----'
print 'pi (true):                        ', pi
print 'pi_known_pi (from tProb_known_pi):', pi_known_pi
print '----'
# calcule the free energy difference of the left basin versus the right basin
p_left = pi_known_pi[0:nbins/2].sum()
delta_f_known_pi = -0.596*np.log(p_left/(1.-p_left))  # \Delta F in kcal/mol
print 'Delta F (true)      : 4.0 kcal/mol'
print 'Delta F (estimate)  :', delta_f_known_pi, 'kcal/mol'
print '----'

# get the timescales from tProb_known_pi
evals, evecs = get_evals_evecs(tProb_known_pi)
implied_timescale_known_pi = max(1.,-lagtime/np.log(evals[1]))
print 'implied timescale (true)        : 9.66e7 +/- 1.5e7 steps'
print 'implied_timescale (estimate)    :', implied_timescale_known_pi, 'steps.',
