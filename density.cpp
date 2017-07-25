#include "density.h"
using namespace std;
#include <gsl/gsl_rng.h>

///////////////////////////////////////////////////////////////////////////////
/// VARIABLE DEFINITION
///////////////////////////////////////////////////////////////////////////////

extern int N_bath; /*!< Size of bath */
extern int Ncut; /*!< Truncation parameter */
double abs_d; /*!< Absolute Value Non-adiabatic Coupling Matrix*/
double Pdotdhat; /*!< Parallel component of momentum*/
double sina;
double cosa;
double de; /*!< */
extern double alpha;

extern double *m; /*!< Mass of particles */
extern double *RR; /*!< */
extern double *PP; /*!< */
double *dhat; /*!< */
extern int  N_slice; /*!< Number of time intervals */
extern double *Pperp; /*!< Perpendicular component of momentum */
extern double TSLICE; /*!< */
extern double *abszsum1;
extern double *argzsum1;
extern double *habszsum1;
extern double *hargzsum1;
extern complex<double> I;


double (*www[2][4][4])(); /*!< Non-adiabatic Coupling Matrix*/
extern double (*phi)(double*, double*); /*!< Density Matrix*/
extern double (*dens_init[4])(double*, double*); /*!< Initial Density Matrix*/
extern double (*obs[4])(double*, double*); /*!< Observable Matrix*/
extern double (*obs1[4])(double*, double*); /*!< Another Observable Matrix*/

extern const gsl_rng_type * TT; /*!< Random Number Generator seed based on Gaussian Distribution */
extern gsl_rng * rr;

complex<double> initd(1,0);


///////////////////////////////////////////////////////////////////////////////
/// PROCESSING TREE
///////////////////////////////////////////////////////////////////////////////

int  density(double *x,double *p){

    int SS0,SS1,SS2,SS3,NNjmp = 0,signPdotdhat;
    double phase0 = 0.0,xx;
    double p0,p1,p2,p3,ap0,ap1,ap2,ap3;
    double dn2;
    complex<double> z = 1.0;
    complex<double> oldz;

    ///////////////////////////////////////////////////////////////////////////////
    /// ALLOCATING MEMORY
    ///////////////////////////////////////////////////////////////////////////////

    dhat = new double[N_bath];


    ///////////////////////////////////////////////////////////////////////////////
    /// INITIALIZATION OF INITIAL SURFACE
    ///////////////////////////////////////////////////////////////////////////////

    gauss_init_W(x, p); /*!< Creates initial random sampling */
    double yy = 4.0*(gsl_rng_uniform (rr));
    /*! Setting initial surface value */
    if (yy < 1.0)
        SS3 = (SS0 = 0);
    else if (yy < 2.0){
        SS0 = 1;
        SS3 = 2;
    }
    else if (yy < 3.0){
        SS0 = 2;
        SS3 = 1;
    }
    else
        SS3 = (SS0 = 3);
    initd = dens_init[SS3](x,p);
    z = 4.0;
    /*! Allocating values for position and momentum from Gaussian */
    for (int l = 0; l < N_bath; ++l){
        RR[l] = x[l];
        PP[l] = p[l];
    }
    SS1 = SS0;
    cout << "Initial Surface " << SS1 << endl;

    ///////////////////////////////////////////////////////////////////////////////
    /// ITERATION FOR TREE: CALCULATE EACH PATH SEGMENT
    ///////////////////////////////////////////////////////////////////////////////
    for (int l = 0; l < N_slice; ++l) {
        cout << "Counter: " << l << endl;
        SS0 = SS1; /*!< Sets beginning surface value */

        ///////////////////////////////////////////////////////////////////////////////
        /// ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        phase0 = U(RR, PP, SS0, TSLICE*0.5); // exp(iLd/2) (before jump)
        z *= exp(I * phase0);

        ///////////////////////////////////////////////////////////////////////////////
        /// NON_ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        dd(RR); // non-adiabatic coupling matrix
        de = dE(RR); /*!< Energy */
        alpha = 0.0;
        Pdotdhat = 0;
        for (int i = 0; i < N_bath; ++i) {
            Pdotdhat += PP[i] * dhat[i];
            alpha += PP[i] * dhat[i] / m[i];
        }
        alpha *= abs_d * TSLICE;
        signPdotdhat = (Pdotdhat < 0 ? -1 : 1); // -1 if neg, 1 if pos
        Pdotdhat = fabs(Pdotdhat);
        for (int i = 0; i < N_bath; ++i)
            Pperp[i] = PP[i] - signPdotdhat * Pdotdhat * dhat[i];
        alpha *= 2.0;
        sina = sin(alpha);
        cosa = cos(alpha);

        /*! Importance Sampling - non-adiabatic coupling matrix gives probabilities */
        ap0 = fabs(p0 = ((www[1][SS0][0]() < -7775.0) ? 0.0 : www[0][SS0][0]()));
        ap1 = fabs(p1 = ((www[1][SS0][1]() < -7775.0) ? 0.0 : www[0][SS0][1]()));
        ap2 = fabs(p2 = ((www[1][SS0][2]() < -7775.0) ? 0.0 : www[0][SS0][2]()));
        ap3 = fabs(p3 = ((www[1][SS0][3]() < -7775.0) ? 0.0 : www[0][SS0][3]()));
        dn2 = ap0 + ap1 + ap2 + ap3;
        xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements
        //alpha goes to 0, pdotdhat very small, matrix becomes identiy and prob of jumping goes to 0
        cout << "Prob:" << "ap0: " << ap0 <<", ap1: "<< ap1 <<", ap2: " << ap2 <<", ap3:"<< ap3 << endl;
        oldz = z;
        SS2 = SS1;
        /*! choose new surface based on probabilities from matrix above
         * calculates new values for the probability weighting
         */
        if (xx < ap0) {
            SS1 = 0;
            z *= p0 * dn2 / ap0;
            cout << SS1 << endl;
        } else if (xx < ap0 + ap1) {
            SS1 = 1;
            z *= p1 * dn2 / ap1;
            cout << SS1 << endl;
        } else if (xx < ap0 + ap1 + ap2) {
            SS1 = 2;
            z *= p2 * dn2 / ap2;
            cout << SS1 << endl;
        } else {
            SS1 = 3;
            z *= p3 * dn2 / ap3;
            cout << SS1 << endl;
        }

        /*! increases jump counter if a jump was undergone,
         * and exiting if jump counter too high (past truncation value)
         */
        if (SS0 != SS1)
            NNjmp++;
        if (NNjmp > Ncut)
            return 0;

        /*! updating momentum values */
        if (www[1][SS0][SS1]() != 9999.0){
            for (int i = 0; i < N_bath; ++i) {
                PP[i] = Pperp[i] + signPdotdhat * www[1][SS0][SS1]() * dhat[i];
            }
        }

        ///////////////////////////////////////////////////////////////////////////////
        /// ADIABATIC PROPAGATOR
        ///////////////////////////////////////////////////////////////////////////////
        phase0 = U(RR,PP,SS1,TSLICE*0.5); // exp(iLd/2) (after jump)
        z *= exp(I*phase0);

        ///////////////////////////////////////////////////////////////////////////////
        /// WRITING SUM VALUES OUT (Solving Integral for Observable and Initial Density)
        ///////////////////////////////////////////////////////////////////////////////
        phi = obs[SS1]; /*!< Observable 1 function */
        abszsum1[l] += real(z*phi(RR,PP)*initd);
        argzsum1[l] += imag(z*phi(RR,PP)*initd);

        phi = obs1[SS1]; /*!< Observable 2 (Hamiltonian) function */
        habszsum1[l] += real(z*phi(RR,PP)*initd);
        hargzsum1[l] += imag(z*phi(RR,PP)*initd);
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// DEALLOCATING MEMEORY
    ///////////////////////////////////////////////////////////////////////////////

    delete [] dhat;

    return 0;
}

