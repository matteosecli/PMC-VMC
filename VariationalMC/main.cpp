/*    CompPhys@SISSA - Variational Monte Carlo Project
 *               Matteo Seclì,  Spring 2017
 *                <secli.matteo@gmail.com>
 *
 *
 * This is part of the second project for the course
 * "Understanding electron correlation with computer
 * simulation" held at SISSA in the AA 2016-2017 by
 * Sandro Sorella and Federico Becca.
 *
 * It's a VMC simulator for a 1D Heisenberg model with PBC.
 * More details in the presentation attached.
 *
 * At every MC step it prints on file the value of the
 * variational parameter, the local energy, the
 * staggered magnetization and the first-neighbor
 * correlation function.
 *
 * TODO: (in descending order of importance)
 *  - Implement the automatic optimization of the
 *    variational parameter (i.e. steepest descent).
 *  - Implement (maybe) some sort of user input. Not
 *    really necessary.
 *
 *
 * The code is hosted at:
 *  <https://github.com/matteosecli/PMC-VMC>
 *
 * Refer to the LICENSE file for license information.
 *
 */

#include <cmath>
#include <armadillo>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <tuple>
#include <ctime>
#include <chrono>
//#include <omp.h>  // For future parallel optimization
#include "lib.h"


using namespace std;
using namespace arma;


unsigned nextElementOf(const unsigned& k, const unsigned& L)
{
    /* === RETURN THE NEXT PBC ELEMENT === */
    return (k+1)%L;
}


double wfRatio(const ivec& x, const unsigned& k, const double& alpha,
               const sp_mat& v, sp_imat& sMatrix, const bool& updateS = false)
{
    /* === CALCULATE THE WF RATIO ===
     *      This function calculates the ratio ψ(x')/ψ(x),
     *      where x' is the configuration obtained from x
     *      by flipping spin k and spin k+1.
     */
    /* Calculate the next element */
        unsigned L = x.n_elem;
        unsigned kNext = nextElementOf(k, L);
    /* Set the ratio of the M signs << "Energy (per site) = " */
        int signM_ratio = -1;
    /* Do the (half) sums in the ratio of the exponentials */
        double expSum = 0;
    /* i = k */
        for ( unsigned j = 0; j < k; j++ )
        {
            if ( j != kNext )
            {
                expSum += -v(k,j)*(double)sMatrix(k,j);
            }
        }
    /* i = kNext */
        for ( unsigned j = 0; j < kNext; j++ )
        {
            if ( j != k )
            {
                expSum += -v(kNext,j)*(double)sMatrix(kNext,j);
            }
        }
    /* j = k */
        for ( unsigned i = k+1; i < L; i++ )
        {
            if ( i != kNext )
            {
                expSum += -v(i,k)*(double)sMatrix(i,k);
            }
        }
    /* j = kNext */
        for ( unsigned i = kNext+1; i < L; i++ )
        {
            if ( i != k )
            {
                expSum += -v(i,kNext)*(double)sMatrix(i,kNext);
            }
        }
    /* Update sMatrix if required */
        if ( updateS == true )
        {
            /* i = k */
                for ( unsigned j = 0; j < k; j++ )
                {
                    if ( j != kNext )
                    {
                        sMatrix(k,j) *= -1;
                    }
                }
            /* i = kNext */
                for ( unsigned j = 0; j < kNext; j++ )
                {
                    if ( j != k )
                    {
                        sMatrix(kNext,j) *= -1;
                    }
                }
            /* j = k */
                for ( unsigned i = k+1; i < L; i++ )
                {
                    if ( i != kNext )
                    {
                        sMatrix(i,k) *= -1;
                    }
                }
            /* j = kNext */
                for ( unsigned i = kNext+1; i < L; i++ )
                {
                    if ( i != k )
                    {
                        sMatrix(i,kNext) *= -1;
                    }
                }
        }

    return (double)signM_ratio*exp(alpha*2.0*expSum);
}


double localEnergy(const ivec& x, const double& alpha,
                   const sp_mat& v, sp_imat& sMatrix, double& firstNCorr)
{
    /* === CALCULATE THE LOCAL ENERGY ===
     *      This function calculates the local energy.
     *      The calculation is achieved by first getting the
     *      contributions from x'!=x, and then x'=x.
     */
        /* Initialize the required variables */
            double eLocal = 0.0;
            double eBuffer = 0.0;
            unsigned L = x.n_elem;
            unsigned kNext = 0;
            unsigned nPairs = 0;
        /* Add the contribution from x'!=x and calculate nPairs */
            for ( unsigned k = 0; k < L; k++ )
            {
                kNext = nextElementOf(k, L);
                if ( x(k) == x(kNext) )
                {
                    nPairs++;
                } else {
                    eBuffer = wfRatio(x, k, alpha, v, sMatrix, false);
                    eLocal += eBuffer;
                }
            }
            eLocal = 0.5*eLocal;
        /* Add the contribution from x'=x */
            firstNCorr = 0.5*(double)nPairs/(double)L - 0.25;
            eLocal += 0.5*(double)nPairs - 0.25*(double)L; // or 0.25*(double)spinProductSum

        /* Cleanup */

    return eLocal;
}


double staggeredMagnetization(const ivec& x)
{
    /* === CALCULATE THE STAGGERED MAGNETIZATION === */
        unsigned L = x.n_elem;
        int spinSum = 0;
        for (unsigned i = 0; i < L; i+=2)
        {
            spinSum += x(i);
            spinSum -= x(i+1);
        }

    return 0.5*(double)abs(spinSum)/(double)L;
}


int main()
{

    /* === INITIALIZE THE REQUIRED VARIABLES ===
     *      LenSim = length of the simulation
     *      L = number of sites (MUST be EVEN!)
     *      alpha* = the variational parameter settings
     *      x = vector of current configuration s_i (total spin 0)
     *      v = the spin Jastrow factor ln(d^2_ij)
     *      sMatrix = the spin products matrix s_i*s_j for the Jastrow
     *      ELocal = the local energy
     *      mStag = the staggered magnetization
     *      firstNCorr = the first-neighbor correlation function
     *      ofile = output file object
     *      outfilename = name of the output file
     *
     *      CAVEAT: v and sMatrix have been defined as sparse
     *              matrices to spare memory. A vector would also
     *              have been fine (a sp_mat is stored as a vector)
     *              but the matrix notation make it easier to
     *              perform certain operations.
    */
    /* Start timers */
        clock_t cpu0 = clock();
        auto wall0 = chrono::system_clock::now();
    /* System and simulation parameters */
        unsigned LenSim = 1000000;
        unsigned L = 10;
    /* Uncomment for interactive input */
        //cin >> LenSim;
        //cin >> L;
    /* Auxiliary variables for the VMC */
        double alphaMin = 0.20;
        double alphaMax = 0.305;
        double alphaStep = 0.01;
        vec alphaRange = regspace<vec>(alphaMin, alphaStep, alphaMax);
        ivec xInit = zeros<ivec>(L);
        sp_mat v = zeros<sp_mat>(L,L);
        sp_imat sMatrixInit = zeros<sp_imat>(L,L);
        double ELocal = 0.0;
        double mStag = 0.0;
        double firstNCorr = 0.0;
    /* Initialize the file variables */
        ofstream ofile;
        string outfilename = "VarMC.dat";
        ofile.open(outfilename);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setprecision(8);
        ofile << "#alpha" << "\t" << "#LocalEnergy" << "\t" << "#m" << "\t" << "#C(1)"<< "\n";
    /* Initialize the first values */
        /* x initialized to up-down-up-down... */
        for ( unsigned i = 0; i < xInit.n_elem/2; i++ )
        {
            xInit(2*i)   = +1;
            xInit(2*i+1) = -1;
        }
        /* v is basically the distance matrix, doesn't change */
        for ( unsigned i = 1; i < v.n_rows; i++ )
        {
            for ( unsigned j = 0; j < i; j++ )
            {
                v(i,j) = 2.0*log(abs(sin(datum::pi*abs((double)(i-j))/(double)L)));
            }
        }
        /* Just the matrix of the products of the spins of x */
        for ( unsigned i = 1; i < sMatrixInit.n_rows; i++ )
        {
            for ( unsigned j = 0; j < i; j++ )
            {
                sMatrixInit(i,j) = xInit(i)*xInit(j);
            }
        }
    /* Write onscreen info */
        cout << "-------------------------------------" << endl;
        cout << "  VARIATIONAL MONTE CARLO SIMULATOR  " << endl;
        cout << "       - 1D HEISENBERG MODEL -       " << endl;
        cout << "-------------------------------------" << endl << endl;
        cout << "Using the following configuration:" << endl <<
                "     Number of sites         = " << L << endl <<
                "     Number of MC steps      = " << LenSim << endl <<
                "     alpha (min, max, step)  = (" << alphaMin << ", "
                                                   << alphaMax << ", "
                                                   << alphaStep << ")"
                                                   << endl << endl;


        /* === DO THE VARIATIONAL LOOP ===
         *      Here we loop on the variational parameter alpha.
         */
        /* Write onscreen info */
            cout << "Looping on alpha:" << endl;
            cout << "     #alpha"
                 << "\t\t" << "#Energy (per site)"
                 << "\t" << "#m"
                 << "\t\t\t" << "#C(1)" << endl;
            cout << fixed << setprecision(6);
        /* Do the loop */
        for ( unsigned alphaIdx = 0; alphaIdx < alphaRange.n_elem; alphaIdx++ )
        {
            /* === INITIALIZATION === */
            /* Assign alpha; initializing inside the loop should help
               with OpenMP */
                double alpha = alphaRange(alphaIdx);
            /* Make a local copy of the initial x and sMatrix */
                ivec x = zeros<ivec>(L);
                sp_imat sMatrix = zeros<sp_imat>(L,L);
                x = xInit;
                sMatrix = sMatrixInit;
            /* Generate the seed */
                long idum = -time(NULL);
            /* Calculate the local energy for the first time */
                double ELocalCumulative = 0.0;
                ELocal = localEnergy(x, alpha, v, sMatrix, firstNCorr);
                ELocalCumulative += ELocal;
            /* Calculate the magnetization for the first time */
                double mStagCumulative = 0.0;
                mStag = staggeredMagnetization(x);
                mStagCumulative += mStag;
            /* Update the 1st NN correlation cumulant */
                double firstNCorrCumulative = 0.0;
                firstNCorrCumulative += firstNCorr;
            /* Write on file the first configuration */
                ofile << alpha << "\t"
                      << ELocal/(double)L << "\t"
                      << mStag << "\t"
                      << firstNCorr << "\n";


            /* === BEGIN MAIN MC LOOP === */
            /* Initialize Metropolis variables */
                unsigned k = 0;
                unsigned kNext = 0;
                double mRatio = 0.0;
                sp_imat sNew= zeros<sp_imat>(L,L);
            /* Do the loop */
            for ( int cycle = 1; cycle < LenSim; cycle++ )
            {
                /* Select a new move (i.e. flip k and k+1) */
                    k = round(ran1(&idum)*(L-1));
                    kNext = nextElementOf(k,L);

                /* Do Metropolis only if we flip different spins,
                   because of the constraint on the total spin. */
                if ( x(kNext) != x(k) )
                {
                    /* Metropolis test */
                    sNew = sMatrix;
                    mRatio = wfRatio(x, k, alpha, v, sNew, true);
                    mRatio = pow(mRatio, 2.0);
                    if ( (mRatio < 1) ? (ran1(&idum) < mRatio) : 1 )
                    {
                        /* Update x by ASSUMING we flip two different spins.
                           The assumption holds because of the "if" that skips
                           this section if we flip two equal spins. */
                            x(k) = -x(k);
                            x(kNext) = -x(kNext);

                        /* Update the sMatrix */
                            sMatrix = sNew;

                        /* Calculate the local energy */
                            ELocal = localEnergy(x, alpha, v, sMatrix, firstNCorr);

                        /* Calculate the staggered magnetization */
                            mStag = staggeredMagnetization(x);
                    }
                }

                /* Update the cumulant of the local energy */
                    ELocalCumulative += ELocal;

                /* Update the cumulant of the staggered magnetization */
                    mStagCumulative += mStag;

                /* Update the cumulant of the 1st NN correlation function */
                    firstNCorrCumulative += firstNCorr;

                /* Write the sample to file */
                    ofile << alpha << "\t"
                          << ELocal/(double)L << "\t"
                          << mStag << "\t"
                          << firstNCorr << "\n";

            } /* End of main MC loop */

            /* Write onscreen info */
                cout << "     " << alpha
                     << "\t\t" << ELocalCumulative/(double)LenSim/(double)L
                     << "\t\t" << mStagCumulative/(double)LenSim
                     << "\t\t" << firstNCorrCumulative/(double)LenSim
                     << endl;

        } /* End of variational loop */
        cout << endl;


    /* === CLEANING SECTION === */
    /* Close the output file */
        ofile.close();
    /* Stop timers and calculate elapsed time */
        clock_t cpu1 = clock();
        auto wall1 = chrono::system_clock::now();
        double cputime = (cpu1 - cpu0) / (double)CLOCKS_PER_SEC;
        double walltime = chrono::duration<double> (wall1 - wall0).count();
    /* Write onscreen info */
        cout << fixed << setprecision(2);
        cout << "Elapsed time:" << endl;
        cout << "     CPU  time (s) = " << cputime  << endl
             << "     WALL time (s) = " << walltime << endl;
        cout << endl
             << "-------------------------------------" << endl;
        cout << "DONE!" << endl;
        cout << "-------------------------------------" << endl;

    return 0;
}

