from matplotlib import pyplot as plt
import numpy as np

import os, sys

if len(sys.argv) < 3:
   print "Usage: python plot_results.py <datafile> stride\n"
   sys.exit(1)

infile = sys.argv[1]
stride = int(sys.argv[2])


data =  np.loadtxt(infile)
plt.figure()
plt.plot(data[::stride,0], data[::stride,1])
plt.xlabel('number of steps')
plt.ylabel('position')
plt.show()

