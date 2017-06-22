#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "random.cpp"
#include   "variable.h"
#include   "variable-trajectory.h"
#include   "density.h"
#include   "functions.h"

using namespace std;

#include <gsl/gsl_rng.h>

//need to define separate variables for each trajectory RR[] with same index as index in RR
//each trajectory/phase dependent variable only defined for relevant path
//each path can be calculated on different/individual thread in parallel
//find a method to index each path - doesn't need to be completely static/known to all threads
//# where path coming from global (known to others), but don't need to know anything else
/*inital condition global, everything else private to each thread,
 * at next splitting point that initial condition global, then everything else private as it propagates onward */



int  density(double *x,double *p) { //density evolution

    int l, i, j, SS0, SS1, SS2, SS3, abc, NNjmp = 0, nojumpflag = -1, signPdotdhat, adiabat_flag = -1, zcut = 0;
    double phase0 = 0.0, phase1 = 0.0, phase3 = 0, xx, yy;
    double p0, p1, p2, p3, ap0, ap1, ap2, ap3, wtemp = 1.0;
    double dn1, dn2, dn3, bb, sqrbb, sqrbb2, sqrbbb, BETA, GAMMA, pbb1, pbb2, pbb3, ABSZ, ARGZ;
    complex<double> z = 1.0, oldz;
    complex<double> I(0, 1);
    double (*phi)(double *, double *, int);

    Dt = TSLICE;

    // Initialization of sample

    gauss_init_W(x, p);
    yy = 4.0 * (gsl_rng_uniform(rr));
    if (yy < 1.0)
        SS3 = (SS0 = 0);
    else if (yy < 2.0) {
        SS0 = 1;
        SS3 = 2;
    } else if (yy < 3.0) {
        SS0 = 2;
        SS3 = 1;
    } else
        SS3 = (SS0 = 3);
    // SS[0] = SS0;
    z = 4.0;
    for (l = 0; l < N_bath; l++) {
        RR[0][l] = x[l];
        PP[0][l] = p[l];
    }
    SS1 = SS0;

    // ____________________________________________________________________

    adiabat_flag = -1;
    for (l = 0; l < N_slice; l++) {
        SS0 = SS1;
        // phase0 = 0.0;
        phase0 = U(RR[0], PP[0], SS0, Dt * 0.5);
        z *= exp(I * phase0);

        //  compute Q1 choose
        //  choose ad or Ua or Ub
        // if Ua or Ub choose kick
        // choose s2 for ad, Ua or Ub
        // note that if ad, s2 = s1
        //  propagate with expLd_{s2}
        // update all prob weights
        // update cummulate for current slice
        // go back to begining

        // ****************************************
        //    compute Q1   and choose s1

        dd(RR[0]); //0 corresponds to the # of ensembles (single path) - needs to be defined for every path
        de = dE(RR[0]);
        abc = 0;
        alpha = 0.0;
        bb = de * Dt * abs_d;

        for (Pdotdhat = 0, i = 0; i < N_bath; i++) {
            Pdotdhat += PP[0][i] * dhat[i];
        }
        alpha = Pdotdhat * abs_d * Dt;
        if (NNjmp >= Ncut)
            adiabat_flag = 1;
        if (adiabat_flag == 1) {
            nojumpflag = 1;
        }

        signPdotdhat = (Pdotdhat < 0 ? -1 : 1);
        Pdotdhat = fabs(Pdotdhat);
        for (i = 0; i < N_bath; i++)
            Pperp[i] = PP[0][i] - signPdotdhat * Pdotdhat * dhat[i];
        /* complete a change of state, means abrupt change of energy
             * what happens at one point */
        alpha *= 2.0;
        sina = sin(alpha);
        cosa = cos(alpha);

        //calculating matrix elements

        //avoid index problem if X < # then = 0, otherwise Y
        ap0 = fabs(p0 = ((www[1][SS0][0]() < -7775.0) ? 0.0 : www[0][SS0][0]()));
        ap1 = fabs(p1 = ((www[1][SS0][1]() < -7775.0) ? 0.0 : www[0][SS0][1]()));
        ap2 = fabs(p2 = ((www[1][SS0][2]() < -7775.0) ? 0.0 : www[0][SS0][2]()));
        ap3 = fabs(p3 = ((www[1][SS0][3]() < -7775.0) ? 0.0 : www[0][SS0][3]()));
        dn2 = ap0 + ap1 + ap2 + ap3;

        xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements
        // printf(" oldz  real %lf imag %lf\t", real(z), imag(z));

        //similar to perturbation expansion
        //keeping track of jumps
        oldz = z;
        SS2 = SS1;
        if (xx < ap0) {
            SS1 = 0;
            z *= p0 * dn2 / ap0;
        } else if (xx < ap0 + ap1) {
            SS1 = 1;
            z *= p1 * dn2 / ap1;
        } else if (xx < ap0 + ap1 + ap2) {
            SS1 = 2;
            z *= p2 * dn2 / ap2;
        } else {
            SS1 = 3;
            z *= p3 * dn2 / ap3;
        }
        if (SS0 != SS1)
            NNjmp++;
        if (NNjmp > Ncut)
            return 0;

        //  if  ((abs(z) > ppower) || (abs(z) >  4.5*pow(1.5,Dt*l))){

        if ((abs(z) > ppower)) {
            if (SS0 != SS1)
                NNjmp--;
            SS1 = SS2;
            z = oldz * www[0][SS0][SS1]();
            // z = oldz;
            // goto jmp;
        }
        if (www[1][SS0][SS1]() != 9999.0)
            for (i = 0; i < N_bath; i++)
                PP[0][i] = Pperp[i] + signPdotdhat * www[1][SS0][SS1]() * dhat[i];


        //  printf("zzz %lf %d  %lf\n", l*Dt, NNjmp, abs(z));

        jmp:
        phase0 = U(RR[0], PP[0], SS1, Dt / 2);
        z *= exp(I * phase0);
        phi = obs[SS1];

        abszsum0[l] = real(z * phi(RR[0], PP[0], 0) * dens_init[SS3](x, p, 1));
        argzsum0[l] = imag(z * phi(RR[0], PP[0], 0) * dens_init[SS3](x, p, 1));
        ABSZ = abs(z * phi(RR[0], PP[0], 0) * dens_init[SS3](x, p, 1));
        ARGZ = arg(z * phi(RR[0], PP[0], 0) * dens_init[SS3](x, p, 1));
        printf("zzzz %d %d %d %lf %lf\t %lf %lf\t %lf %lf\t %lf %lf\n", l, NNjmp, SS1, real(z), imag(z), abs(z), arg(z),
               ABSZ, ARGZ, alpha, de);
        realsum[l][NNjmp] += abszsum0[l];
        imagsum[l][NNjmp] += argzsum0[l];
        abszsum1[l] += abszsum0[l];
        argzsum1[l] += argzsum0[l];

        hist[l][NNjmp]++;

        phi = obs1[SS1];
        habszsum0[l] = real(z * phi(RR[0], PP[0], 0) * dens_init[SS3](x, p, 1));
        hargzsum0[l] = imag(z * phi(RR[0], PP[0], 0) * dens_init[SS3](x, p, 1));
        hrealsum[l][NNjmp] += habszsum0[l];
        himagsum[l][NNjmp] += hargzsum0[l];
        habszsum1[l] += habszsum0[l];
        hargzsum1[l] += hargzsum0[l];
    }
return 0;
}
