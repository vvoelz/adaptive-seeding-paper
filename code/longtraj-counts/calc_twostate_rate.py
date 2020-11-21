import os, sys, glob
import numpy as np


### NOTE ###
"""
Each of these datafiles are 10k x 10k arrays that contain one continuous trajectory.
It can be accessed by the np.concatenate() transformation

The time step of the sampling was 0.05, and snapshots were saved every 100 steps (lagtime = 100 steps)

"""
infiles = glob.glob('/home/tuf74538/for_vince/DHAM-detailed_balance/data/*.npy')

counts = np.zeros( (20,20) )

tally = 0
for infile in infiles:
    print 'Loading', infile, '...',
    f = np.load(infile) 
    print 'concatenating...',
    g = np.concatenate(f)  # a 100 million-frame traectory!!!

    print 'converting...', 
    # convert the trajectory to 20 bins
    h = ((g -1.5)*20.0/4.0).astype(int)

    print '...compiling counts...'
    for i in range(h.shape[0]-1):
        counts[h[i],h[i+1]] += 1.0

    outfile = 'counts_20bins_%d.npy'%tally
    print '...saving', outfile, '...',
    np.save('counts_20bins_%d.npy'%tally,counts)  # save cumulative tally
    print '...Done.'

    tally += 1
