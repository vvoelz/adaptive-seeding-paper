import os, sys
import  numpy as np

# See also: http://matplotlib.org/users/customizing.html

import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt

fontfamily={'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)

######## Import the two-state potential ########
from twostate import *

from BinManager import *
nbins = 20
bins = BinManager(nbins=nbins)

dx = 0.01
xmin, xmax = 1.5, 5.5
x = np.arange(xmin, xmax+dx, dx)

plt.figure(figsize=(3.3,3))

# plot the bin borders
for e in bins.edges:
    plt.plot([e,e],[-1.,10], '-', color='#808080', lw=0.6)

plt.plot(x,U(x))
plt.xticks([1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5], fontsize = 10)
plt.xlabel('$x$')
plt.ylabel('$U(x)$ (kcal/mol)')
plt.xlim(xmin, xmax)
plt.ylim(-1,10)

plt.tight_layout()
plt.savefig('1D-potential.pdf')


