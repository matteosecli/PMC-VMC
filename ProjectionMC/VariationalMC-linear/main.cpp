/*    CompPhys@SISSA - Projection Monte Carlo Project
 *               Matteo Seclì,  Spring 2017
 *                <secli.matteo@gmail.com>
 *
 *
 * This is part of the second project for the course
 * "Understanding electron correlation with computer
 * simulation" held at SISSA in the AA 2016-2017 by
 * Sandro Sorella and Federico Becca.
 *
 * It's a VMC simulator for the same system described
 * in the PMC simulator. It's written just to compare
 * results with both methods.
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
#include <ctime>
#include <chrono>
#include "lib.h"


using namespace std;
using namespace arma;


double localEnergy(int& x, double& expAlpha, double& t, double& V, int& L)
{
    /* === CALCULATE THE LOCAL ENERGY ===
     *      This function calculates the local energy.
     */
    if ( x == 1)
    {
        return -t/expAlpha + V*(double)x;
    } else if ( x == L ) {
        return -t*expAlpha + V*(double)x;
    } else {
        return -t*(expAlpha+1.0/expAlpha) + V*(double)x;
    }
}

double MetropolisRatio(int& x, int& xNext, double& alpha, int& L)
{
    /* === CALCULATE THE WF RATIO TIMES THE TRANSITION P ===
     *      This function calculates the ratio ψ(x')/ψ(x),
     *      where x' = xNext.
     *      Then, it multiplies its square by the transition
     *      probability in order to obtain the ratio used
     *      in the Metropolis algorithm.
     */
    double wfRatio = exp(-2.0*alpha*(xNext-x));
    double transP = 1.0;
    if ( x == 1 || x == L ) transP = 0.5;
    if ( ( x == 2 && xNext == 1 ) || ( x == L-1 && xNext == L ) ) transP = 2.0;

    return wfRatio*transP;
}

int main()
{

    /* === INITIALIZE THE REQUIRED VARIABLES ===
     *      LenSim = length of the simulation
     *      L = number of sites
     *      t = hopping amplitude
     *      V = on-site potential amplitude (to be multiplied by site idx)
     *      x = idx of our position
     *      ELocal = the local energy
     *      idum = the usual crap for the random generator
     *      ofile = output file object
     *      outfilename = name of the output file
    */
    /* Start timers */
        clock_t cpu0 = clock();
        auto wall0 = chrono::system_clock::now();
    /* System and simulation parameters */
        int LenSim = 100000000;
        int L = 20;
        double t = 1.0;
        double V = 1.0;
    /* Auxiliary variables for the PMC */
        double alphaMin = 0.65;
        double alphaMax = 0.65;
        double alphaStep = 0.01;
        vec alphaRange = regspace<vec>(alphaMin, alphaStep, alphaMax);
        int x = 0;
        double ELocal = 0.0;
        long idum = -time(NULL);
    /* Initialize the file variables */
        ofstream ofile;
        string outfilename = "VarMC.dat";
        ofile.open(outfilename);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setprecision(8);
    /* Write onscreen info */
        cout << "-------------------------------------" << endl;
        cout << "  VARIATIONAL MONTE CARLO SIMULATOR  " << endl;
        cout << "-------------------------------------" << endl << endl;
        cout << "Using the following configuration:" << endl <<
                "     Number of sites         = " << L << endl <<
                "     Number of MC steps      = " << LenSim << endl <<
                "     alpha (min, max, step)  = (" << alphaMin << ", "
                                                   << alphaMax << ", "
                                                   << alphaStep << ")"
                                                   << endl << endl <<
                "     Hopping amplitude   (t) = " << t << endl <<
                "     Potential amplitude (V) = " << V << endl << endl;


    /* === DO THE VARIATIONAL LOOP ===
     *      Here we loop on the variational parameter alpha.
     */
    /* Write onscreen info */
        cout << "Looping on alpha:" << endl;
        cout << "     #alpha"
             << "\t\t" << "#Energy" << endl;
        cout << fixed << setprecision(6);
    /* Do the loop */
    for ( unsigned alphaIdx = 0; alphaIdx < alphaRange.n_elem; alphaIdx++ )
    {
        /* === INITIALIZATION === */
        /* Assign alpha; initializing inside the loop should help
           with OpenMP */
            double alpha = alphaRange(alphaIdx);
            double expAlpha = exp(alpha);
        /* Generate the seed */
            long idum = -time(NULL);
        /* Calculate local energy and position for the first time */
            double ELocalCumulative = 0.0;
            x = (int) round(ran1(&idum)*(L-1))+1;
            ELocal = localEnergy(x, expAlpha, t, V, L);
            ELocalCumulative += ELocal;
        /* Write on file the first configuration */
            ofile << "#alpha" << "\t" << "#Position" << "\t" << "#LocalEnergy" << "\n";
            ofile << alpha << "\t" << x << "\t" << ELocal << "\n";
        /* === BEGIN MAIN MC LOOP === */
        /* Initialize Metropolis variables */
            int xNext = 0;
            double mRatio = 0.0;
        /* Do the loop */
        for ( int cycle = 1; cycle < LenSim; cycle++ )
        {
            /* Select a new move */
                if ( x == 1)
                {
                    xNext = 2;
                } else if ( x == L ) {
                    xNext = L-1;
                } else {
                    xNext = x + (int) copysign(1.0, ran1(&idum)-0.5);
                }
            /* Metropolis test */
                mRatio = MetropolisRatio(x, xNext, alpha, L);
                if ( (mRatio < 1) ? (ran1(&idum) < mRatio) : 1 )
                {
                     /* Update x */
                        x = xNext;
                     /* Calculate the local energy */
                        ELocal = localEnergy(x, expAlpha, t, V, L);
                }
             /* Update the cumulant of the local energy */
                ELocalCumulative += ELocal;
             /* Write the sample to file */
                ofile << alpha << "\t"
                      << x << "\t"
                      << ELocal << "\n";
         } /* End of main MC loop */
         /* Write onscreen info */
            cout << "     " << alpha
                 << "\t\t" << ELocalCumulative/(double)LenSim
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
        cout << "Elapsed time:" << endl;
        cout << "     CPU  time (s) = " << cputime  << endl
             << "     WALL time (s) = " << walltime << endl;
        cout << endl
             << "-------------------------------------" << endl;
        cout << "DONE!" << endl;
        cout << "-------------------------------------" << endl;

    return 0;
}

