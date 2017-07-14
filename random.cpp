#include "random.h"
using namespace std;
#include <gsl/gsl_rng.h>

#define PI 3.141592653589793

extern int N_bath;
extern const gsl_rng_type * TT;
extern gsl_rng * rr;
extern double *sig;

double ranVector[10001];





/* Random number stuff ***************************************************************  */


double  gauss1(double sigma_x, int i){
    double x1,y1,y2,y3;

    y1 = ranVector[i];

    while (fabs(y1) < 1.0e-200){
        y1 = gsl_rng_uniform (rr);
    }

    y2 = ranVector[i+N_bath];
    y3 = sqrt(-2*log(y1));
    x1 = y3*cos(2*PI*y2);
    return (sigma_x*x1);
}

void randnums(int rand_dim, double *rand_vec){
    for (int i = 0; i < rand_dim; ++i){
        rand_vec[i] = gsl_rng_uniform (rr);
    }
}

void gauss_init_W(double *R, double *v){
    double sigma_x, sigma_v;
    /* Gaussian number generator  for (R,P) */

    randnums(4*N_bath, ranVector);
    for (int i = 0; i < N_bath; ++i){
        sigma_x = sig[i];
        sigma_v = sig[i+N_bath];
        R[i] = gauss1(sigma_x,i);
        v[i] = gauss1(sigma_v,i + 2*N_bath);
    }
}
