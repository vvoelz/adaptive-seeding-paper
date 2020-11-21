import numpy as np

class BinManager(object):
    """A class to manage 1D bin partitions"""
    
    def __init__(self, xmin=1.5, xmax=5.5, nbins=10):
        """Initialize the class with uniform bins."""
        
        self.xmin = xmin
        self.xmax = xmax
        self.nbins = nbins
        self.nedges = nbins + 1
        
        # Define initial bin_edges to be uniformly spaced
        dx = (xmax - xmin)/nbins
        self.edges = np.arange(self.xmin, self.xmax + dx, dx)
        assert len(self.edges) == self.nedges
        
        # bin widths
        self.update_widths()
        
        # bin centers
        self.update_centers()
        
        # neighbor list
        self.update_neighbors()
        

    def update_widths(self):
        """Calculate the bin widths from the bin edges"""
        self.widths = self.edges[1:] - self.edges[0:self.nedges-1]
            
    def update_centers(self):
        """Calculate the bin centers from the bin edges"""
        self.centers = (self.edges[1:] + self.edges[0:self.nedges-1])/2.0
        
    def update_neighbors(self):
        """Calculate the neighbors."""
        self.neighbors = []
        self.neighbors.append([1]) # the left-most bin (index 0) has only 1 neighbor
        for i in range(1, self.nbins-1):
            self.neighbors.append([i-1,i+1]) 
        self.neighbors.append([self.nbins-2])
        
    def bin_index(self, x):
        """returns the bin index of an input value x.
        
        NOTES
        Will return     -1      if outside the left bin edge
                        nbins   if outside the right bin edge.
        """
        
 
        if not isinstance(x, np.ndarray):
            # x is a number
            if x < self.edges[0]:
                return -1
            elif x > self.edges[-1]:
                return self.nbins
            else:
                i = 0
                while i < self.nbins:
                    if (x >= self.edges[i]) and (x < self.edges[i+1]):
                        return i
                    i += 1

        else:
            # x is an array!  ###
            a = ((x - self.xmin)*self.nbins/(self.xmax - self.xmin)).astype(int)
            
            # correct for any bin indices < 0 or  >= self.nbins:
            a2 = np.minimum((self.nbins-1)*np.ones( a.shape, dtype=int), a)
            return np.maximum(np.zeros( a.shape, dtype=int), a2)

              
    def split(self, i):
        """split bin i into two states.
        
        NOTE: states will be renumbered after this!"""
        
        assert i >= 0
        assert i < self.nbins
        
        left_edge, right_edge = self.edges[i], self.edges[i+1]
        
        # the bin center becomes a new edges
        new_edges = np.hstack( (self.edges[0:i+1], self.centers[i:i+1], self.edges[i+1:]) )
        print 'new_edges', new_edges
        
        # replace old center with the midpoints of the split bin
        left_center  = (self.centers[i] + left_edge)/2.0
        right_center = (right_edge + self.centers[i])/2.0
        new_centers = np.hstack( (self.centers[0:i], np.array([left_center, right_center]), self.centers[i+1:]) )
        print 'new_centers', new_centers
        
        self.nbins += 1
        self.nedges += 1
        self.edges = new_edges
        self.centers = new_centers
        self.update_widths()
        self.update_neighbors()

        
