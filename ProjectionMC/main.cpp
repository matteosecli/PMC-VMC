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
 * It's a PMC simulator for a particle on a 1D lattice
 * with L sites, connected by a hopping amplitude and with
 * an on-site potential that grows linearly in the direction
 * of the chain. More details in the presentation attached.
 *
 * At every MC step it prints on file the current position
 * of the particle, the local energy and a weight b_x which
 * is described in the book by Becca&Sorella. The weight
 * allows for the calculation of the GS energy during the
 * post-processing, as detailed in the book; a practical
 * implementation is given in the attached 'Evsp.m', which
 * plots the GS energy as a function of the bin length 'p'.
 * Then, one takes the lowest 'p' that gives a converged
 * energy (in our case, it's basically always 100). The
 * final estimation for the ground state is calculated by
 * 'finalerrors.m' which takes the mean value of the GS energy
 * around p=100 (in order to kill potential fluctuations).
 *
 * In the program, the diagonal matrix Λ required for the
 * projection procedure is chosen in such a way that the
 * diagonal elements are λ = V*L. In this way, the Green's
 * function is always positive and we don't have the sign
 * problem.
 *
 * TODO: (in descending order of importance)
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
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <ctime>
#include <chrono>
#include "lib.h"


using namespace std;


void GUtils( int& x, int& L, double& t, double& V, double& Lambda,
             double& GTransPS, double& GTransOS,
             double& GTransNS, double& bx )
{
    /* === CALCULATE THE PMC UTILITIES ===
     *      This function calculates the G matrix elements
     *      and the normalization bx (their sum).
    */
    /* G matrix elements and their sum (bx) */
        if ( x == 1 )
        {
            GTransPS = 0;
            GTransNS = t;
        } else if ( x == L ) {
            GTransPS = t;
            GTransNS = 0;
        } else {
            GTransPS = t;
            GTransNS = t;
        }
        GTransOS = Lambda - V*x;
        bx = GTransPS + GTransOS + GTransNS;
}

void MoveSelector( int& x, double& bx, long& idum,
                   double& GTransPS, double& GTransOS )
{
    /* === SELECT A NEW MOVE ===
     *      xTrans = the probability of the move times bx
    */
        double xTrans = 0;
        xTrans = ran1(&idum)*bx;
        if ( xTrans < GTransPS )
        {
            x = x-1;
        } else if ( xTrans < GTransPS+GTransOS ) {
        } else {
            x = x+1;
        }
}

int main()
{

    /* === INITIALIZE THE REQUIRED VARIABLES ===
     *      LenSim = length of the simulation
     *      L = number of sites
     *      t = hopping amplitude
     *      V = on-site potential amplitude (to be multiplied by site idx)
     *      Lambda = diagonal element of the Lambda matrix
     *      bx = scale factor (or normalization) at site x
     *      GTransPS, GTransNS, GTansOS = G_{x',x} with x' (respectively)
     *           previous, on and next site
     *      x = idx of our position
     *      ELocal = Lambda - bx, the local energy
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
        double V = 2.0;
        double Lambda = V*L;
    /* Auxiliary variables for the PMC */
        double bx = 0;
        double GTransPS = 0;
        double GTransOS = 0;
        double GTransNS = 0;
        int x = 0;
        double ELocal = 0;
        long idum = -time(NULL);
    /* Initialize the file variables */
        ofstream ofile;
        string outfilename = "ProjMC.dat";
        ofile.open(outfilename);
        ofile << setiosflags(ios::showpoint | ios::uppercase);
        ofile << setprecision(8);
    /* Initialize the first values */
        x = (int) round(ran1(&idum)*(L-1))+1;
        GUtils(x, L, t, V, Lambda, GTransPS, GTransOS, GTransNS, bx);
        ELocal = Lambda - bx ;
    /* Write onscreen info */
        cout << "------------------------------------" << endl;
        cout << "  PROJECTION MONTE CARLO SIMULATOR  " << endl;
        cout << "------------------------------------" << endl << endl;
        cout << "Using the following configuration:" << endl <<
                "     Number of sites         = " << L << endl <<
                "     Number of MC steps      = " << LenSim << endl <<
                "     Hopping amplitude   (t) = " << t << endl <<
                "     Potential amplitude (V) = " << V << endl << endl;
    /* Write on file the first configuration */
        ofile << "#Position" << "\t" << "#LocalEnergy" << "\t" << "#bx" << "\n";
        ofile << x << "\t" << ELocal << "\t" << bx << "\n";


    /* === BEGIN MAIN MC LOOP === */
    for ( int cycle = 1; cycle < LenSim; cycle++ )
    {
        /* Select a new move */
            MoveSelector( x, bx, idum, GTransPS, GTransOS );

        /* Recalculate the G matrix terms and the bx */
            GUtils(x, L, t, V, Lambda, GTransPS, GTransOS, GTransNS, bx);

        /* Calculate the local energy */
            ELocal = Lambda - bx ;

        /* Write position, local energy and bx */
            ofile << x << "\t" << ELocal << "\t" << bx << "\n";

    } /* End of main MC loop */


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
             << "------------------------------------" << endl;
        cout << "DONE!" << endl;
        cout << "------------------------------------" << endl;

    return 0;
}

