#include "functions.h"

using namespace std;

// =========================================================================
// VARIABLES
// =========================================================================

extern int N_bath;
extern double ddd;
extern double ddd4;
extern double delta;
extern double abs_d;
extern double timestep;
extern double *dgam;
extern double *dhat;
extern double *mww;
extern double *f;
extern double *c;
extern double *m;
extern double *w;

extern double cosa;
extern double sina;
extern double Pdotdhat;
extern double de;

extern double (* www[2][4][4])();

extern void (*force[4])(double *);

// =========================================================================
// FUNCTIONS
// =========================================================================

double gam(double *R){
    double x = 0.0;
    double asyEps = 0.4;

    for (int i = 0; i < N_bath; ++i)
        x += c[i]*R[i];
    x += asyEps;    // asymmetric spin boson
    return -x;
}

double Hb(double *R, double *P){ /* Bath Hamiltonian */
    double x = 0.0;
    for (int i = 0; i < N_bath; ++i)
        x += P[i]*P[i] - mww[i]*R[i]*R[i];
    return x*0.5;
}

void dgamma(double *R){
    for (int i = 0; i < N_bath; ++i)
        dgam[i] = -c[i];
}

void Fb(double *R){ /* Pure Bath Force Field */
    double x;
    for (int i= 0; i < N_bath; ++i)
        f[i] = mww[i]*R[i];
}

void F1(double *R){ /* 00 force field   */
    double g,h;
    g = gam(R);
    h = g/sqrt(ddd4 + g*g);
    for (int i = 0; i < N_bath; ++i){
        f[i]  = mww[i]*R[i] -  h*c[i];
    }
}

void F2(double *R){ /* 11 force field */
    double g,h;
    g = gam(R);
    h = g/sqrt(ddd4 + g*g);
    for (int i = 0; i < N_bath; ++i)
        f[i] = mww[i]*R[i] + h*c[i];
}

double dE(double *R){ /* Energy difference between adibiatic surface (E1 - E0) */
    double g;
    g = gam(R);
    g *= 4.0*g;
    return (sqrt(ddd + g));
}

double G(double *R){
    double x,g;
    g = gam(R);
    if (fabs(g/delta) < 1.0e-7)
        return (g/delta);
    x = (-delta + sqrt(ddd + 4*g*g))/(2*g);
    return x;
}

void dd(double*R){
    double x1,x2,x3;
    int i;
    x2 = gam(R);
    dgamma(R);
    if (fabs(x2) < 1.0e-4)
        x3 = 1/delta;
    else {
        x1 = G(R);
        x3 = -x1/x2 + 2.0/(delta + 2.0*x2*x1);
        x3 = x3/(1.0 + x1*x1);
    }
    for (i = 0,abs_d = 0.0; i < N_bath; ++i){
        dhat[i] = -dgam[i]*x3;   // 1-3-05 put - sign in
        abs_d += dhat[i]*dhat[i];
    }
    abs_d = sqrt(abs_d);
    for (i = 0; i < N_bath; ++i)
        dhat[i] /= abs_d;
}

void integ_step(double *r, double *v, double dt, int Sa){ // Velocity Verlet
    double y;
    y = 0.5*dt*dt;
    for (int i = 0; i < N_bath; ++i)
        r[i] += dt*v[i] + y*f[i];
    y = 0.5*dt;
    for (int i = 0; i < N_bath; ++i)
        v[i] += y*f[i];
    force[Sa](r);
    for (int i = 0; i < N_bath; ++i)
        v[i] += y*f[i];
}

void bath_para(double eta, double w_max){ /* Parameters for bath (corresponding to an ohmic spectral density) */
    double w_0;
    w_0 = (1 - exp(-w_max))/N_bath;
    for (int i = 0; i < N_bath; ++i){
        m[i] = 1.0;
        w[i] = -log( 1-(i+1)*w_0 );
        c[i] = sqrt(eta*w_0*m[i])*w[i];
    }
}

double U( double *r,double *v, int Sa, double t){ // Adiabatic Propagator
    double  dE0, phase,dt,x1,x2,x3,v1,v2,v3;
    int Nsteps;

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
        for (int i = 0; i < Nsteps; ++i){
            integ_step(r , v,  dt, Sa);
        }
        return 0.0;
    }
    phase = dE(r)*0.5;
    for (int i = 0; i < Nsteps; ++i){
        integ_step(r , v,  dt, Sa);
        phase += dE(r);
    }
    phase -=dE(r)*0.5;
    phase*= dt;

    if (Sa == 1)
        phase *= -1.0;

    return phase;
}

// TRANSITION MATRIX ===========================================================


/* Q1 */


double wwa0_00(){
    double x;
    x = 1 + cosa;
    return x*0.5;
}

double wwa0_01(){
    double x;
    x = -sina;
    return x*0.5;
}

double wwa0_02(){
    double x;
    x = -sina;
    return x*0.5;
}

double wwa0_03(){
    double x;
    x = 1.0 - cosa;
    return x*0.5;
}

double wwa0_10(){
    double x;
    x = sina;
    return x*0.5;
}

double wwa0_11(){
    double x;
    x = 1.0 + cosa;
    return x*0.5;
}

double wwa0_12(){
    double x;
    x = -1.0 + cosa;
    return x*0.5;
}

double wwa0_13(){
    double x;
    x = -sina;
    return x*0.5;
}

double wwa0_20(){
    double x;
    x = sina;
    return x*0.5;
}

double wwa0_21(){
    double x;
    x = -1.0 + cosa ;
    return x*0.5;
}

double wwa0_22(){
    double x;
    x = 1.0 + cosa;
    return x*0.5;
}

double wwa0_23(){
    double x;
    x = -sina ;
    return x*0.5;
}

double wwa0_30(){
    double x;
    x = 1 - cosa;
    return x*0.5;
}

double wwa0_31(){
    double x;
    x =  sina;
    return x*0.5;
}

double wwa0_32(){
    double x;
    x = sina;
    return x*0.5;
}

double wwa0_33(){
    double x;
    x = 1 + cosa;
    return x*0.5;
}

//         W_{a 1}

/* _____________________________________________  */

double wwa1_00(){
    return 9999.0;
}

double wwa1_01(){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_02(){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_03(){
    double x;
    x = Pdotdhat*Pdotdhat - 2.0*de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_10(){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_11(){
    return 9999.0;
}

double wwa1_12(){
    return 9999.0;
}

double wwa1_13(){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_20(){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_21(){
    return 9999.0;

}

double wwa1_22(){
    return 9999.0;
}

double wwa1_23(){
    double x;
    x = Pdotdhat*Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_30(){
    double x;
    x = Pdotdhat*Pdotdhat + 2.0*de;
    return sqrt(x);
}

double wwa1_31(){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_32(){
    double x;
    x = Pdotdhat*Pdotdhat + de;
    return sqrt(x);
}

double wwa1_33(){
    return 9999.0;
}

void setwww(){

    // W_{a0}

    www[0][0][0] = wwa0_00;
    www[0][0][1] = wwa0_01;
    www[0][0][2] = wwa0_02;
    www[0][0][3] = wwa0_03;
    www[0][1][0] = wwa0_10;
    www[0][1][1] = wwa0_11;
    www[0][1][2] = wwa0_12;
    www[0][1][3] = wwa0_13;
    www[0][2][0] = wwa0_20;
    www[0][2][1] = wwa0_21;
    www[0][2][2] = wwa0_22;
    www[0][2][3] = wwa0_23;
    www[0][3][0] = wwa0_30;
    www[0][3][1] = wwa0_31;
    www[0][3][2] = wwa0_32;
    www[0][3][3] = wwa0_33;

    // W_{a1}

    www[1][0][0] = wwa1_00;
    www[1][0][1] = wwa1_01;
    www[1][0][2] = wwa1_02;
    www[1][0][3] = wwa1_03;
    www[1][1][0] = wwa1_10;
    www[1][1][1] = wwa1_11;
    www[1][1][2] = wwa1_12;
    www[1][1][3] = wwa1_13;
    www[1][2][0] = wwa1_20;
    www[1][2][1] = wwa1_21;
    www[1][2][2] = wwa1_22;
    www[1][2][3] = wwa1_23;
    www[1][3][0] = wwa1_30;
    www[1][3][1] = wwa1_31;
    www[1][3][2] = wwa1_32;
    www[1][3][3] = wwa1_33;

}


/* ______________________________________   */

/* ************** Observables and initial density Matrices *********  */


/* Definition of initial density matrix element */


double wigner_harm_osc(double *x, double *p){
    return 1.0;
}

double dens_init_0(double *x,double *p){
    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = 0.5*(1.0 + 2.0*g + gg)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_1(double *x,double *p){
    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_2(double *x,double *p){
    double z;
    double g, gg;
    g = G(x); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_3(double *x,double *p){
    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = 0.5*(gg - 2*g  + 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double obs_0(double *x,double *p){
    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = 2.0*g/(1 + gg);
    return z;
}

double obs_1(double *x,double *p){
    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = (gg-1)/(1 + gg);
    return z;
}

double obs_2(double *x,double *p){
    double z;
    double g, gg;
    g = G(x); gg = g*g;
    z =  (gg-1)/(1 + gg);
    return z;
}

double obs_3(double *x,double *p){
    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = -2.0*g/(1 + gg);
    return z;
}

/* These are the matrix elements of the Hamiltonian */

double H_0(double *x,double *p){
    double z;
    z = Hb(x,p) - dE(x)*0.5;
    return z;
}

double H_1(double *x,double *p){
    return 0.0;
}

double H_2(double *x,double *p){
    return 0.0;
}

double H_3(double *x,double *p){
    double z;
    z = Hb(x,p) + dE(x)*0.5;
    return z;
}

