#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "variable.h"
#include   "variable-trajectory.h"
#include   "functions.h"
using namespace std;

#include <gsl/gsl_rng.h>

double gam(double *R){ //depends on coords. of bath

    double asyEps = 0.4;
    double x;
    int i;

    for (i = 0, x = 0.0; i < N_bath; i++)
        x += c[i]*R[i]; // - coupling potential energy ('02 after eq32)
    x += asyEps;    // asymmetric spin boson
    return -x;
}


double Hb(double *R, double *P){ /* Bath Hamiltonian */ //pure harmonic oscillator

    double x;
    int i;

    for (i = 0, x = 0.0; i < N_bath; i++)
        x += P[i]*P[i] - mww[i]*R[i]*R[i];
    return x/2.0;
}


void dgamma(double *R){ //derivative of gam - may not want to assume bath is harmonic

    int i;

    for (i = 0; i < N_bath; i++)
        dgam[i] = -c[i];
}

void Fb(double *R){ /* Pure Bath Force Field */ // '02 eq 38

    int i;
    double x;

    for (i= 0; i < N_bath; i++)
        f[i] = mww[i]*R[i];
}

void F1(double *R){ /* 00 force field  */  //'02 eq 39.1

    int i;
    double g,h;

    g = gam(R);
    h = g/sqrt(ddd4 + g*g);
    for (i = 0; i< N_bath; i++){
        f[i]  = mww[i]*R[i] -  h*c[i];
    }
    /* ok  note mww[i] = - m*w[i]*w[i] as calculated in main */
}

void F2(double *R){ /* 11 force field */ // '02 eq 39.2

    int i;
    double g,h;

    g = gam(R);
    h = g/sqrt(ddd4 + g*g);

    for (i = 0; i< N_bath; i++)
        f[i] = mww[i]*R[i] + h*c[i];
}

double dE(double *R){ /* Energy difference between adibiatic surface (E1 - E0) */

    double g;

    g = gam(R);
    g *= 4.0*g;
    return (sqrt(ddd + g));
}


double G(double *R){ // '02 eq 40

    double x,g;

    g = gam(R);
    if (fabs(g/delta) < 1.0e-7)
        return (g/delta);
    x = (-delta + sqrt(ddd + 4*g*g))/(2*g);
    return x;
}

void dd(double*R){ //non-adiabatic coupling matrix

    double x1,x2,x3;
    int i;

    /* Computing d = <1|d/dr|0>  */
    x2 = gam(R);
    dgamma(R);
    if ( fabs(x2) < 1.0e-4)
        x3 = 1/delta;
    else {
        x1 = G(R);
        x3 = -x1/x2 + 2.0/(delta + 2.0*x2*x1);
        x3 = x3/(1.0 + x1*x1);
    }
    for (i = 0, abs_d = 0.0; i < N_bath; i++){
        dhat[i] = -dgam[i]*x3;
        abs_d += dhat[i]*dhat[i];
    }
    abs_d = sqrt(abs_d);
    for (i = 0; i < N_bath; i++)
        dhat[i] /= abs_d;
    /* Note may need to change the sign of dhat here

    /* dhat now is not normalized, and equal d_{10}(R) */
}

void integ_step(double *r, double *v, double dt, int Sa){ //classical propagator

    /* ********* Velocity Verlet ************************** */

    double y;
    int i;

    y = 0.5*dt*dt;
    for (i = 0; i < N_bath; i++)
        r[i] += dt*v[i] + y*f[i];
    y = 0.5*dt;
    for (i = 0; i < N_bath; i++)
        v[i] += y*f[i];
    force[Sa](r); //Hellman-Feynman Forces
    for (i = 0; i < N_bath; i++)
        v[i] += y*f[i];
}


void bath_para(double eta, double w_max){

    /* Parameters for bath (corresponding to an ohmic spectral density) */

    int i;
    double w_0;

    w_0 = (1 - exp(-w_max))/N_bath;

    for (i = 0; i < N_bath; i++){
        m[i] = 1.0;
        w[i] = -log( 1-(i+1)*w_0 );
        c[i] = sqrt(eta*w_0*m[i])*w[i];
    }
}


double U( double *r,double *v, int Sa, double t){ //propagator //exp(+-iLdel)

    double  dE0, phase, dt, x1, x2, x3, v1, v2, v3;
    int Nsteps, i;

    /* ******** Adiabatic Propagator *********************  */


    force[Sa](r);
    dt = timestep;

    if (t <= timestep){
        dt = t;
        Nsteps = 1;
    }
    else{
        Nsteps = t/dt +1;
        dt = t/Nsteps;
    }

    if ((Sa == 0) || (Sa == 3)){
        for (i = 0; i < Nsteps; i++){
            integ_step(r , v,  dt, Sa);
        }
        return 0.0;
    }
    //collecting phase info from H-F forces
    phase = dE(r)*0.5;
    for (i = 0; i < Nsteps; i++){
        integ_step(r , v,  dt, Sa);
        phase += dE(r);
    }
    phase -=dE(r)*0.5;
    phase*= dt;

    if (Sa == 1)
        phase *= -1.0;

    return phase;

    /* For some reason before I was returing  - phase  */

    /*   note dE = E1 - E0 > 0 , w_10 = dE , corresponding to the 1 0 pair
     of surfaces, equivalent to Sa = 2    */


}