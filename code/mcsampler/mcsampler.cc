#include <iostream>
#include <fstream>
#include <stdio.h>      /* printf */
#include <stdlib.h>     /* strtol */

#include <cstdlib>
#include <cmath>

#include <ctime>

#include <time.h>


double rnd()
{
      return (double)rand() / ( RAND_MAX );
}

double U(double x)
{
    /* Returns the value of the 1D potential energy surface in kcal/mol:
    U(x) = -\frac{2k_BT}{0.596} \ln [ e^{−2(x−2)^2-2} + e^{−2(x−5)^2} ] 
    kT = 0.596  - the thermal energy in units of kcal/mol */
    
    return -2.0*log( exp(-2.0*pow(x-2.0,2) - 2.0) + exp(-2.0*pow(x-5.0,2) ));
}

float* sample(double xinit, long nsteps, long nstride, char* outfile, double lambda)
{

 /* 
 sample(U, xinit, nsteps, djump=0.05, xmin=1.5, xmax=5.5,
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
    
    Note:  djump=0.005 parameters are from the 2017 Stelzl et al. paper    -- we are doing 10x this step
    """ */

    /* params */ 
    double djump=0.05, xmin=1.5, xmax=5.5, kT=0.596;
    long nprint=1000; 

    double x, xnew, energy, new_energy, acc_ratio;
    long step, accepted_steps;
    bool accept;
 
    x = xinit;
    energy = lambda*U(x);
    
    step = 0;
    accepted_steps = 0;


    float* traj = NULL;   // Pointer to float array, initialize to nothing. Store as floats -- values will range from 1.5 to 5.5
    traj = new float[ nsteps/nstride ];  // Allocate nsteps/nstride ints and save ptr in traj.
    for (int i=0; i<(nsteps/nstride); i++) {
        traj[i] = 0.0;    // Initialize all elements to zero.
    }
    long itraj = 0;
    
    while (step < nsteps)
    {    
        xnew = x + djump*(2.0*rnd()-1.0);
        new_energy = lambda*U(xnew);
        
        /* calculate Metropolis acceptance */
        accept = (rnd() < fmin(1.0, exp( -1.0*(new_energy-energy)/kT ) ));
        
        /* reject moves that bring x outside the range */
        accept = accept && (xnew>xmin) && (xnew<xmax);
        
        if (accept)
        {
            accepted_steps += 1;
            x = xnew;
            energy = lambda*U(x);
        }          
        if ((step % nstride) == 0)
        {
            traj[itraj] = x;
            itraj += 1;
        }    
        if ((step % nprint) == 0)
        {
            printf("step %10ld of %10ld nsteps: x = %7.4lf, energy = %7.4lf\n", step, nsteps, x, energy);
        }                                                         

        step += 1;
        acc_ratio = (double)accepted_steps/(double)step;
    }    

    /* write the trajectory to file */
    printf("Writing to file %s ...", outfile);

    std::ofstream myfile;
    myfile.open(outfile);
    myfile << "#step\tx\n";
    for (int i=0; i<(nsteps/nstride); i++)
    {
        myfile << nstride*i << "\t" << traj[i] << "\n";
    }
    myfile.close();

    printf("...DONE\\n");

    return traj;
}               

 
int main(int argc, char** argv)
{

    std::cout << "You have entered " << argc
         << " arguments:" << "\n";

    if (argc != 7) 
    {
        std::cout << "Usage: mcsampler <lambda> <nsteps> <nstride> <xinit> <outfile> randseed" << "\n";
        return 1;
    }
   
    float lambda = strtof(argv[1], NULL);
    long nsteps = strtof(argv[2], NULL);
    long nstride = strtof(argv[3], NULL);
    float xinit = strtof(argv[4], NULL);

    long myseed = strtof(argv[6], NULL);

    printf("lambda = %7.4lf ; nsteps =  %10ld ; nstride = %10ld ; xinit = %7.4lf; seed = %10ld \n", lambda, nsteps, nstride, xinit, myseed);


    /* reseed the random number generator */
    /* std::srand(std::time(0)); */
    srand(myseed);

    sample((double)xinit, nsteps, nstride, argv[5], (double)lambda);

 
    return 0;
}


