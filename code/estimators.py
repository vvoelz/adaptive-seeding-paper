import numpy as np


def rownormalized_tProb(C, symmetrize=False, verbose=False):
    """Estimates a transition probability matrix tProb
    from a row-nomalized  counts matrix of elements C_ij where i-->j.
    
    RETURNS
    tProb - transition probabiliity matrix with elements T_ij    for i --> j 
    
    """
    m, n = C.shape[0], C.shape[1]
    assert m == n
    
    # based on row-normalized counts
    tProb = np.zeros( (m,n) )
    if symmetrize:
        C_sym = (C + np.transpose(C))/2.0  # symmetrize to make evals real
        rowsums = np.transpose(np.tile( C_sym.sum(axis=1), (n,1)) )
        return C_sym/rowsums
    else:
        rowsums = np.transpose(np.tile( C.sum(axis=1), (n,1)) )
        return C/rowsums

    # print 'tProb', tProb, 'tProb.shape', tProb.shape


def MLE_tProb_reversible(C, verbose=False, maxsteps=1000):
    """Estimates a transition probability matrix tProb
    from a count matrix, using a maximum-likelihood estimator (Eq. 29 and 30 of Wu et al. J Chem. Phys 2014)

    INPUTS
    C     - a count matrix of elements C_ij where i-->j
    
    RETURNS
    tProb - transition probabiliity matrix with elements T_ij    for i --> j 
    
    """

    # The count matrix must be sqaure
    m, n = C.shape[0], C.shape[1]
    assert m == n

    N_i = C.sum(axis=1)

    # set initial guess for eq pops
    pi = np.ones(n)/float(n)

    error = 1.0e99
    tolerance = 1e-6

    Csym = C + C.transpose()   # each element Csym[i,j] is C[i,j]+C[j,i]
    
    step = 0
    maxsteps_reached = False
    while (error > tolerance) and (maxsteps > step):
        
        if (verbose):
            print 'step', step, 'error', error, 'tolerance', tolerance
        
        # estimate new tProbs
        new_pi = np.zeros( pi.shape )
        for i in range(n):
            for j in range(n):
                new_pi[i] += Csym[i,j]/(N_i[i]/pi[i] + N_i[j]/pi[j]) 
     
        error = np.abs(new_pi - pi).sum()
        pi = new_pi
        step += 1

        if maxsteps == step:
            maxsteps_reached = True

    tProb = np.zeros( (n,n) )
    for i in range(n):
        for j in range(n):
            tProb[i,j] = Csym[i,j]*pi[j]/(pi[i]*N_i[j] + pi[j]*N_i[i])

    if maxsteps_reached:
        print 'WARNING: maxsteps =', maxsteps, 'reached without convergence (tol = %f)'%tolerance
    return tProb
 


def MLE_tProb_dTRAM(C_k, gamma, verbose=False):
    """Estimates a transition probability matrix tProb using DTRAM
    from a series (k=1...K) of transition counts C_ij^(k) where i-->j taken from
    thermodynamic ensemble (k), and thermodynamic biases \gamma_i^(k) for each bin k 
    
    Uses Eq. 16 and 17 of Wu et al JCP 2004
    
    INPUTS
    
    C_k   - an array of shape (K, N, N) where C_k[k,i,j] = C_ij^(k)
    gamma - an array of shape (K, N) containing the bias weights \gamma_i^(k) for each state index i
    
    RETURNS
    tProb - a series of transition probabiliity matrix with elements T_ij^(k) for i --> j 
            stored in array of shape (K, N, N).
    
    """
    
    # the count matrices must be square, with n states
    M, N = C_k.shape[1], C_k.shape[2]
    assert M == N
    
    # the number of thermodynamic ensembles must match
    K, L = C_k.shape[0], gamma.shape[0]
    assert K == L

  
    # set initial guess for eq pops
    pi = np.ones( N )/float(N)
  
    # set initial guesses for Lagrange v_i^(k)
    v_i = C_k.sum(axis=2)
    v_i += np.ones(v_i.shape)  # make sure there are no zeros
    
    error = 1.0e99
    tolerance = 1e-12
    maxsteps = 10000
    
    Csym2 = np.zeros( C_k.shape )
    for k in range(K):
        Csym2[k] = C_k[k] + C_k[k].transpose()   # each element Csym2[k,i,j] is C_k[k,i,j]+C_k[k,j,i]

    maxsteps_reached = False    
    step = 0
    while (error > tolerance) and (maxsteps > step):
        
        if (verbose):
            print 'step', step, 'err', error, 'tol', tolerance, 'v_i', v_i, 'pi', pi
        
        # estimate new_v_i
        new_v_i = np.zeros( (K,N) )
        for i in range(N):
            jsum = 0.0
            for j in range(N):
                jsum += Csym2[:,i,j]*gamma[:,j]*pi[j]/(gamma[:,i]*pi[i]*v_i[:,j] +
                                                       gamma[:,j]*pi[j]*v_i[:,i])
            new_v_i[:,i] = v_i[:,i]*jsum                 

        # convert the nans to zero
        new_v_i[np.isnan(new_v_i)] = 0.0
        
        if (verbose):
            print 'new_v_i', new_v_i
        

        # estimate new_pi 
        new_pi = np.zeros( pi.shape )
        for i in range(N):      
            jsum = np.zeros(K)
            for j in range(N):
                #numerator = Csym2[:,i,j]*v_i[:,j]
                #denominator = (gamma[:,i]*pi[i]*v_i[:,j] + gamma[:,j]*pi[j]*v_i[:,i])
                #Ind = np.isnan(denominator)
                jsum += Csym2[:,i,j]*v_i[:,j]/(gamma[:,i]*pi[i]*v_i[:,j] + gamma[:,j]*pi[j]*v_i[:,i])
            print 'jsum', jsum
            new_pi[i] = C_k[:,:,i].sum()/jsum.sum()
 
        print 'new_pi', new_pi
    
        v_i_error = np.abs(new_v_i - v_i).sum()           
        pi_error = np.abs(new_pi - pi).sum()

        error = v_i_error + pi_error

        pi = new_pi
        v_i = new_v_i

        step += 1

        if maxsteps == step:
            maxsteps_reached =  True
 
    tProb = np.zeros( (K, N, N) )
    for i in range(N):  
        for j in range(N):
            tProb[:,i,j] = Csym2[:,i,j]*gamma[:,j]*pi[j]/(gamma[:,i]*pi[i]*v_i[:,j] + 
                                                          gamma[:,j]*pi[j]*v_i[:,i])

    if maxsteps_reached:
        print 'WARNING: maxsteps =', maxsteps, 'reached without convergence (tol = %f)'%tolerance
    return tProb



def MLE_tProb_known_pi(C, pi, verbose=False):
    """Estimates a transition probability matrix tProb
    from a count matrix, using a maximum-likelihood estimator 
    and known equilibirum populations. (Trendelkamp-Schroer and Noe, Phys. Rev. X 2016)

    INPUTS
    C     - a count matrix of elements C_ij where i-->j
    pi    - a known equilibrium distribution
    
    RETURNS
    tProb - transition probabiliity matrix with elements T_ij    for i --> j 
    
    """

    # The count matrix must be sqaure
    m, n = C.shape[0], C.shape[1]
    assert m == n

    # the equil pops vector must be the same size
    assert pi.shape[0] == n

    # make an informed guess for the tProb and initial populations
    # based on row-normalized counts
    tProb = rownormalized_tProb(C, symmetrize=True)

    lam = np.ones(n)   # initial guess for the lambdas

    error = 1.0e99
    tolerance = 1e-12
    maxsteps = 100

    step = 0
    while (error > tolerance) and (maxsteps > step):

        if (verbose):
            print 'step', step, 'error', error, 'tolerance', tolerance

        # estimate new tProbs
        new_tProb = np.zeros( tProb.shape )
        for i in range(n):
            for j in range(n):
                new_tProb[i,j] = (C[i,j] + C[j,i])*pi[j]/(lam[i]*pi[j] + lam[j]*pi[i])

        # new_pi = get_stationary_pops(new_tProb)

        # estimate new lams
        new_lam = np.ones( lam.shape )
        for i in range(n):
            new_lam[i] =  new_tProb[i,:].sum()*lam[i]

        error = np.abs(new_tProb - tProb).sum()

        tProb = new_tProb
        lam = new_lam

        step += 1

    if maxsteps <= step:
        print 'WARNING: maxsteps =', maxsteps, 'reached without convergence (tol = %f)'%tolerance
    else:
        return tProb
 

    
def get_stationary_pops(tProb):
    """Returns the stationary eigenvector of a tProb matrix
    (where T_ij holds transiton probs i-->j) """
    
    evals, evecs = get_evals_evecs(tProb)
    return evecs[:,0]/evecs[:,0].sum()

def get_evals_evecs(tProb):
    """Returns *sorted* eigenvalues and eigenvectors (largest to smallest)
    of a tProb matrix where T_ij holds transiton probs i-->j
    
    RETURNS
    evals, evecs
    """
    
    evals, evecs =  np.linalg.eig(np.transpose(tProb))
    
    # sort the evecs by largest eval
    Ind = np.argsort(-evals) # sort largest-to-smallest
    evecs = evecs[:,Ind]
    evals = evals[Ind]
    
    return evals, evecs

def MLE_tProb_dTRAM_single_ensemble(C, verbose=False):
    """Estimates a transition probability matrix tProb using DTRAM
    from a transitions counts matrix of elements C_ij where i-->j.
    
    Uses Eq. 16 and 17 of Wu et al JCP 2004
    
    NOTE: *** this gives exactly the same result as MLE_tProb_reversible ***
    *** i.e. MLE_tProb_reversible is a special case when the number of ensembles K = 1 ***
    
    RETURNS
    tProb - transition probabiliity matrix with elements T_ij    for i --> j 
    
    """
    
    m, n = C.shape[0], C.shape[1]
    assert m == n
    
    # set initial guess for eq pops
    pi = np.ones(n)/float(n)
  
    # set initial guess for Lagrange v_i
    v_i = C.sum(axis=1)
    
    error = 1.0e99
    tolerance = 1e-12
    maxsteps = 10000
    
    Csym2 = C + C.transpose()   # each element Csym2[i,j] is C[i,j]+C[j,i]

    maxsteps_reached = False    
    step = 0
    while (error > tolerance) and (maxsteps > step):
        
        if (verbose):
            print 'step', step, 'err', error, 'tol', tolerance, 'v_i', v_i, 'pi', pi
        
        # estimate new_v_i
        new_v_i = np.zeros( n )
        for i in range(n):
            jsum = 0.0
            for j in range(n):
                jsum += Csym2[i,j]*pi[j]/(pi[i]*v_i[j] + pi[j]*v_i[i])
            new_v_i[i] = v_i[i]*jsum                 

        # estimate new_pi 
        new_pi = np.zeros( pi.shape )
        for i in range(n):       
            jsum = 0.0
            for j in range(n):
                jsum += Csym2[i,j]*v_i[j]/(pi[i]*v_i[j] + pi[j]*v_i[i])
            new_pi[i] = C[:,i].sum()/jsum
 
        v_i_error = np.abs(new_v_i - v_i).sum()           
        pi_error = np.abs(new_pi - pi).sum()

        error = v_i_error + pi_error

        pi = new_pi
        v_i = new_v_i

        step += 1

        if maxsteps == step:
            maxsteps_reached =  True
 
    tProb = np.zeros( (n,n) )
    for i in range(n):  
        for j in range(n):
            tProb[i,j] = Csym2[i,j]*pi[j]/(pi[i]*v_i[j] + pi[j]*v_i[i])

    if maxsteps_reached:
        print 'WARNING: maxsteps =', maxsteps, 'reached without convergence (tol = %f)'%tolerance
    return tProb





