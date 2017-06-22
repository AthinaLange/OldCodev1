#ifndef VARIABLE
#define VARIABLE

#include <gsl/gsl_rng.h>

//separate out trajectory-dependent variables

int N_bath;
int  N_slice, Nsample, Ncut;

const gsl_rng_type * TT;
gsl_rng * rr;

double (*dens_init[4])(double *, double *, int );
double (*obs[4])(double *, double *, int );
double (*obs1[4])(double *, double *, int );
double (* www[6][4][4])();
void (*force[4])(double *);


int *SS, **hist; //could use local hist and accumulate into a global hist
double *mww;
double *meann, *sig;
double *m /*mass - set to 1.0 in bath_para*/, *c/*system para */, *w /*frequency 0 system para*/, *d /*phase-space dependent */, delta /*timestep*/;
double timestep, TSLICE, Dt /*integrating timestep */;


#endif

//may want to try and separate out physical system