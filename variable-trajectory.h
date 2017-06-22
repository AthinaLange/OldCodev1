#ifndef TRAJECTORY
#define TRAJECTORY

#include <gsl/gsl_rng.h>

double ddd4, ddd, ppower, de; //indexed by trajectory index
double **RR, **PP, *dhat, *f, *dgam, *Pperp /*perp element of momentum to d - position dependent*/; //trajectory dependent
double abs_d /*trajectory dependent*/, Pdotdhat /*component of momentum parallel to d (function) trajectory-dependent*/;

double sina, cosa, alpha; //trajectory dependent

double  *abszsum0, *abszsum1, *argzsum0, *argzsum1, **realsum, **imagsum; //0 - summand trajectory dependent ; 1 - sum
double  *habszsum0, *habszsum1, *hargzsum0, *hargzsum1, **hrealsum, **himagsum;


#endif
