#ifndef TRANSITION
#define TRANSITION

#include <complex>
#include <vector>
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include <iostream>

double gam(double *R);

double Hb(double *R, double *P);

void dgamma(double *R);

void Fb(double *R);

void F1(double *R);

void F2(double *R);

double dE(double *R);

double G(double *R);

void dd(double*R);

void integ_step(double *r, double *v, double dt, int Sa);

void bath_para(double eta, double w_max);

double U(double *r,double *v, int Sa, double t);

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
