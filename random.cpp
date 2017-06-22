#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
//#include   "random.h"
#include   "variable.h"

using namespace std;

#include <gsl/gsl_rng.h>

#define PI 3.141592653589793

double ranVector[10001];

/* Random number stuff ***************************************************************  */


double  gauss1(double mean_x, double sigma_x, int i){

    double x1,y1,y2,y3;

    y1 = ranVector[i];

    while (fabs(y1) < 1.0e-200){
        y1 = gsl_rng_uniform (rr);
    }

    y2 = ranVector[i+N_bath];
    y3 = sqrt(-2*log(y1));
    x1 = y3*cos(2*PI*y2);
    return (mean_x = sigma_x*x1);
}

void randnums(int rand_dim, double *rand_vec){

    int i;
    for (i = 0; i < rand_dim; i++){
        rand_vec[i] = gsl_rng_uniform (rr); //Random number between 0 and 1 inputed into rand_vec
    }
}


void gauss_init_W(double *R, double *v){ /* Gaussian number generator  for (R,P) */
//sample from gaussian given mean and variance
//takes into account zero-point motion (in x and p) given a variance //assuming harmonic oscillators
//used in Density part of the code
    double sigma_x, sigma_v, mean_x, mean_v;
    int i;

    randnums(4*N_bath, ranVector);
    for (i = 0; i< N_bath; i++){
        mean_x = meann[i];
        sigma_x = sig[i];
        // sigma_x = 1.0;
        mean_v = 0;
        sigma_v = sig[i+N_bath];
        // sigma_v = 1.0;
        R[i] = gauss1(mean_x, sigma_x,i);
        v[i] = gauss1(mean_v, sigma_v,i + 2*N_bath);
    }
}