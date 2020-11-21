from pyemma.msm.estimators import OOMReweightedMSM
from pyemma.msm import estimate_markov_model

import numpy as np


"""
method 1
"""
#msmrev=OOMReweightedMSM(lag=150,sparse=True,reversible=False,rank_Ct='bootstrap_trajs')
#msmrev=OOMReweightedMSM(lag=150,sparse=True,reversible=False)

#tol_rank=10.0  or smaller? 
#sparse=True/False
#reversible=True


"""
method2
"""
sequence=np.load('all_faked_trajs_0.npy')

dtrajs=[sequence[i] for i in range(len(sequence))]

#msmrev_fit=msmrev.fit(dtrajs)


msm = estimate_markov_model(dtrajs, lag=200, weights='oom')
np.save('msm_timescales.npy',msm.timescales())

#msm.stationary_distribution
