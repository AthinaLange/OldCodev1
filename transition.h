#ifndef TRANSITION
#define TRANSITION

#include <complex>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>


void setwww();

double dens_init_0(double *x,double *p);

double dens_init_1(double *x,double *p);

double dens_init_2(double *x,double *p);

double dens_init_3(double *x,double *p);

double obs_0(double *x,double *p);

double obs_1(double *x,double *p);

double obs_2(double *x,double *p);

double obs_3(double *x,double *p);

double H_0(double *x,double *p);

double H_1(double *x,double *p);

double H_2(double *x,double *p);

double H_3(double *x,double *p);

#endif 
