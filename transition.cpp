#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "variable.h"
#include   "variable-trajectory.h"
#include   "functions.h"
#include   "transition.h"
using namespace std;

#include <gsl/gsl_rng.h>

double cosb1,cosb2, sinb1,sinb2,cosg1,cosg2,sing1,sing2;

//    Transition Matrices ___________________________________________________________________________-

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
    return  x*0.5;
}

double wwa0_20(){
    double x;
    x =  sina;
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
    x = 1.0 - cosa;
    return x*0.5;
}

double wwa0_31(){
    double x;
    x =  sina;
    return x*0.5;
}

double wwa0_32(){
    double x;
    x =  sina;
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
    x = Pdotdhat * Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_02(){
    double x;
    x = Pdotdhat * Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_03(){
    double x;
    x = Pdotdhat * Pdotdhat - 2.0*de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_10() {
    double x;
    x = Pdotdhat * Pdotdhat + de;
    return sqrt(x);
}

double wwa1_11(){
    return 9999.0;
}

double wwa1_12(){
    return 9999.0;
}

double wwa1_13() {
    double x;
    x = Pdotdhat * Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_20(){
    double x;
    x = Pdotdhat * Pdotdhat + de;
    return sqrt(x);
}

double wwa1_21(){
    double x;
    return 9999.0;
}

double wwa1_22(){
    double x;
    return 9999.0;
}

double wwa1_23(){
    double x;
    x = Pdotdhat * Pdotdhat - de;
    if (x <= 0)
        return -7777.0;
    else
        return sqrt(x);
}

double wwa1_30(){
    double x;
    x = Pdotdhat * Pdotdhat + 2.0*de;
    return sqrt(x);
}

double wwa1_31(){
    double x;
    x = Pdotdhat * Pdotdhat + de;
    return sqrt(x);
}

double wwa1_32(){
    double x;
    x = Pdotdhat * Pdotdhat + de;
    return sqrt(x);
}

double wwa1_33(){
    return 9999.0;
}

//   W_{a2}

/* _____________________________________________  */

double wwa2_00(){
    double x;
    x = sina - sinb2;
    return x;
}

double wwa2_01(){
    double x;
    x = cosb2;
    return x;
}

double wwa2_02(){
    double x;
    x = cosb2;
    return x;
}

double wwa2_03(){
    double x;
    x = sina  + sinb2;
    return  x;
}

double wwa2_10(){
    double x;
    x = cosa;
    return x;
}

double wwa2_11(){
    double x;
    return 0.0;
}

double wwa2_12(){
    double x;
    return 0.0;
}

double wwa2_13(){
    double x;
    x = cosa;
    return x;
}

double wwa2_20(){
    double x;
    x = cosa;
    return x;
}

double wwa2_21(){
    double x;
    return 0.0;
}

double wwa2_22(){
    double x;
    return 0.0;
}

double wwa2_23(){
    double x;
    x = cosa;
    return x;
}

double wwa2_30(){
    double x;
    x = sina + sinb2;
    return -x;
}

double wwa2_31(){
    double x;
    x = cosb2;
    return x;
}

double wwa2_32(){
    double x;
    x = cosb2;
    return x;
}

double wwa2_33(){
    double x;
    x = -sina + sinb2;
    return x;
}

//   Wb0

double wwb0_00(){
    double x;
    x = cosa;
    return x*x;
}

double wwb0_01(){
    double x;
    x = cosa * sina;
    return x;
}

double wwb0_02(){
    double x;
    x = cosa * sina;
    return x;
}

double wwb0_03(){
    double x;
    x = sina;
    return 1.0 + x*x;
}

double wwb0_10(){
    double x;
    x = -cosa * sina;
    return x;
}

double wwb0_11(){
    double x;
    x = cosa;
    return x*x;
}

double wwb0_12(){
    double x;
    x = cosa;
    return x*x;
}

double wwb0_13(){
    double x;
    x = cosa * sina;
    return -x;
}

double wwb0_20() {
    double x;
    x = cosa * sina;
    return -x;
}

double wwb0_21(){
    double x;
    x = cosa;
    return x*x;
}

double wwb0_22(){
    double x;
    x = cosa;
    return x*x;
}

double wwb0_23(){
    double x;
    x = cosa * sina;
    return x;
}

double wwb0_30(){
    double x;
    x = sina;
    return 1.0 + x*x;
}

double wwb0_31(){
    double x;
    x = cosa * sina;
    return -x;
}

double wwb0_32(){
    double x;
    x = -cosa * sina;
    return x;
}

double wwb0_33(){
    double x;
    x = cosa;
    return x*x;
}


//  W_{b1}

double wwb1_00(){
    double x;
    x = 1.0 - sina*sing1;
    return x;
}

double wwb1_01(){
    double x;
    x = cosg1 * sina;
    return x;
}

double wwb1_02(){
    double x;
    x = cosg1 * sina;
    return x;
}

double wwb1_03(){
    double x;
    x = 1.0 + sina*sing1;
    return x;
}

double wwb1_10(){
    double x;
    x = cosa * sing1;
    return -x;
}

double wwb1_11(){
    double x;
    x = cosa * cosg1;
    return x;
}

double wwb1_12(){
    double x;
    x = cosa * cosg1;
    return x;
}

double wwb1_13(){
    double x;
    x = cosa * sing1;
    return x;
}

double wwb1_20(){
    double x;
    x = cosa * sing1;
    return -x;
}

double wwb1_21(){
    double x;
    x = cosa * cosg1;
    return x;
}

double wwb1_22(){
    double x;
    x = cosa * cosg1;
    return x;
}

double wwb1_23(){
    double x;
    x = cosa * sing1;
    return x;
}

double wwb1_30(){
    double x;
    x = 1.0 + sina*sing1;
    return x;
}

double wwb1_31(){
    double x;
    x = -cosg1 * sina;
    return x;
}

double wwb1_32(){
    double x;
    x = -cosg1 * sina;
    return x;
}

double wwb1_33(){
    double x;
    x = 1.0 - sina*sing1;
    return x;
}


//  W_{b2}

double wwb2_00(){
    double x;
    x = 1.0 - sina*sing2;
    return x;
}

double wwb2_01(){
    double x;
    x = cosg2 * sina;
    return x;
}

double wwb2_02(){
    double x;
    x = cosg2 * sina;
    return x;
}

double wwb2_03(){
    double x;
    x = 1.0 + sina*sing2;
    return x;
}

double wwb2_10(){
    double x;
    x = cosa * sing2;
    return -x;
}

double wwb2_11(){
    double x;
    x = cosa * cosg2;
    return x;
}

double wwb2_12(){
    double x;
    x = cosa * cosg2;
    return x;
}

double wwb2_13(){
    double x;
    x = cosa * sing2;
    return x;
}

double wwb2_20(){
    double x;
    x = cosa * sing2;
    return -x;
}

double wwb2_21(){
    double x;
    x = cosa * cosg2;
    return x;
}

double wwb2_22(){
    double x;
    x = cosa * cosg2;
    return x;
}

double wwb2_23(){
    double x;
    x = cosa * sing2;
    return x;
}

double wwb2_30(){
    double x;
    x = 1.0 + sina*sing2;
    return x;
}

double wwb2_31(){
    double x;
    x = cosg2 * sina;
    return -x;
}

double wwb2_32(){
    double x;
    x = cosg2 * sina;
    return -x;
}

double wwb2_33(){
    double x;
    x = 1.0 - sina*sing2;
    return x;
}

/*---------------------------------------------------------*/

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

    // W_{a2}

    www[2][0][0] = wwa2_00;
    www[2][0][1] = wwa2_01;
    www[2][0][2] = wwa2_02;
    www[2][0][3] = wwa2_03;
    www[2][1][0] = wwa2_10;
    www[2][1][1] = wwa2_11;
    www[2][1][2] = wwa2_12;
    www[2][1][3] = wwa2_13;
    www[2][2][0] = wwa2_20;
    www[2][2][1] = wwa2_21;
    www[2][2][2] = wwa2_22;
    www[2][2][3] = wwa2_23;
    www[2][3][0] = wwa2_30;
    www[2][3][1] = wwa2_31;
    www[2][3][2] = wwa2_32;
    www[2][3][3] = wwa2_33;

    //  wb0

    www[3][0][0] = wwb0_00;
    www[3][0][1] = wwb0_01;
    www[3][0][2] = wwb0_02;
    www[3][0][3] = wwb0_03;
    www[3][1][0] = wwb0_10;
    www[3][1][1] = wwb0_11;
    www[3][1][2] = wwb0_12;
    www[3][1][3] = wwb0_13;
    www[3][2][0] = wwb0_20;
    www[3][2][1] = wwb0_21;
    www[3][2][2] = wwb0_22;
    www[3][2][3] = wwb0_23;
    www[3][3][0] = wwb0_30;
    www[3][3][1] = wwb0_31;
    www[3][3][2] = wwb0_32;
    www[3][3][3] = wwb0_33;

    // W_{b1}

    www[4][0][0] = wwb1_00;
    www[4][0][1] = wwb1_01;
    www[4][0][2] = wwb1_02;
    www[4][0][3] = wwb1_03;
    www[4][1][0] = wwb1_10;
    www[4][1][1] = wwb1_11;
    www[4][1][2] = wwb1_12;
    www[4][1][3] = wwb1_13;
    www[4][2][0] = wwb1_20;
    www[4][2][1] = wwb1_21;
    www[4][2][2] = wwb1_22;
    www[4][2][3] = wwb1_23;
    www[4][3][0] = wwb1_30;
    www[4][3][1] = wwb1_31;
    www[4][3][2] = wwb1_32;
    www[4][3][3] = wwb1_33;

    // W_{b2}

    www[5][0][0] = wwb2_00;
    www[5][0][1] = wwb2_01;
    www[5][0][2] = wwb2_02;
    www[5][0][3] = wwb2_03;
    www[5][1][0] = wwb2_10;
    www[5][1][1] = wwb2_11;
    www[5][1][2] = wwb2_12;
    www[5][1][3] = wwb2_13;
    www[5][2][0] = wwb2_20;
    www[5][2][1] = wwb2_21;
    www[5][2][2] = wwb2_22;
    www[5][2][3] = wwb2_23;
    www[5][3][0] = wwb2_30;
    www[5][3][1] = wwb2_31;
    www[5][3][2] = wwb2_32;
    www[5][3][3] = wwb2_33;

}


/* ************** Observables and initial density Matrices *********  */


/* Definition of initial density matrix element */


double wigner_harm_osc(double *x, double *p){
    double prod,y1,y2;
    int i;
    /* set wigner to one has it drops out of calculation for our initial
    density */
    return 1.0;
}


double dens_init_0(double *x,double *p, int k){
    double z;
    double g,gg;
    g = G(x); gg = g*g;
    z = 0.5*(1.0 + 2.0*g + gg)/(1 + gg);

    return (z*wigner_harm_osc(x,p));


}

double dens_init_1(double *x,double *p, int k){
    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_2(double *x,double *p, int k){
    double z;
    double g, gg;

    g = G(x); gg = g*g;
    z = 0.5*(gg - 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double dens_init_3(double *x,double *p, int k){
    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = 0.5*(gg - 2*g  + 1.0)/(1 + gg);
    return (z*wigner_harm_osc(x,p));
}

double obs_0(double *x,double *p, int k){
    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = 2.0*g/(1 + gg);
    return z;
}

double obs_1(double *x,double *p, int k){
    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = (gg-1)/(1 + gg);
    return z;
}

double obs_2(double *x,double *p, int k){
    double z;
    double g, gg;

    g = G(x); gg = g*g;
    z =  (gg-1)/(1 + gg);
    return z;
}

double obs_3(double *x,double *p, int k){
    double z;
    double g,gg;

    g = G(x); gg = g*g;
    z = -2.0*g/(1 + gg);
    return z;
}

/* These are the matrix elements of the Hamiltonian */

double H_0(double *x,double *p, int k){
    double z;
    double g,gg;

    z = Hb(x,p) - dE(x)/2.0;
    return z;
}

double H_1(double *x,double *p, int k){
    return 0.0;
}

double H_2(double *x,double *p, int k){
    return 0.0;
}

double H_3(double *x,double *p, int k){
    double z;
    double g,gg;

    z = Hb(x,p) + dE(x)/2.0;
    return z;
}