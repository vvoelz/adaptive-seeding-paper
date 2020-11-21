# Functions for the two-state model

import numpy as np

# This is a test of browser-based commits

def U(x, kT=0.596):
    """Return the value of the 1D potential energy surface in kcal/mol:
    
    U(x) = -\frac{2k_BT}{0.596} \ln [ e^{-2(x-2)^2-2} + e^{-2(x-5)^2} ].      (1)
    
    PARAMS
    kT   - the thermal energy in units of kcal/mol 
    """
    return (-2.0*kT/0.596)*np.log( np.exp(-2.0*(x-2)**2 - 2) + np.exp(-2.0*(x-5)**2) )

def sample(U, xinit, nsteps, djump=0.05, xmin=1.5, xmax=5.5,
           kT=0.596, nstride=100, nprint=10000, verbose=False):
    """Perform Monte Carlo sampling of the potential energy surface U
    by 
    
    INPUT
    U        a supplied potential energy function 
    x        the starting position
    nsteps   number of steps of Monte Carlo to perform
    
    PARAMS
    djump    attempt random moves drawn from [-djump, +djump]
    xmin     reject moves x < xmin
    xmax     reject moves x > xmax
    kT       thermal energy in units of kcal/mol (Default: 0.596)
    nstride  frequency of step to subsample the trajectory
    
    Note:  djump=0.005 parameters are from the 2017 Stelzl et al. paper    
    """
    
    x = xinit
    energy = U(x)
    
    step = 0
    accepted_steps = 0
    traj = np.zeros( nsteps/nstride )
    itraj = 0
    
    # pre-calculate random numbers
    r = np.random.random( nsteps )
    s = np.random.random( nsteps )

    while step < nsteps:
        
        xnew = x + djump*(2.0*s[step]-1.0)
        new_energy = U(xnew)
        
        # calculate Metropolis acceptance 
        accept = (r[step] < min(1, np.exp( -1.0*(new_energy-energy)/kT ) ))
        
        # reject moves that bring x outside the range
        accept = accept*(xnew>xmin)*(xnew<xmax)
        
        if accept:
            accepted_steps += 1
            x = xnew
            energy = U(x)
                  
        if step%nstride == 0:
            traj[itraj] = x
            itraj += 1
            
        if verbose:
            if step%nprint == 0:
                print 'step', step, 'of', nsteps, ': x =', x, 'energy =', energy
                                                             
        step += 1
        acc_ratio = float(accepted_steps)/float(step)
        
    return traj
               
