#include "functions.h"
//#include "transition.h"


using namespace std;

// VARIABLES ================================================================

extern double cosa;
extern double sina;
extern double Pdotdhat;
extern double de;

extern double (* www[2][4][4])();


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
