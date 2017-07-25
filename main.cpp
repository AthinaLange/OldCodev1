#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "density.h"
#include   "functions.h"

using namespace std;
#include <gsl/gsl_rng.h>

/*!
 *  \brief Using a Trotter Approximation to calculate a quantum-classical non-adiabatic approximation to
 * separate the quantum and the quasi-classical degrees of freedom to allow a surface-hopping scheme to
 * be implemented that creates a tree-like structure.
 *
 * This code follows one of the possible paths through the tree structure. The non-adiabatic propagator
 * determines the likelihood of a jump based on importance sampling
 *
 * \author Donal MacKernan and Athina Lange.
 * \date 24.7.17
 */





///////////////////////////////////////////////////////////////////////////////
/// SYSTEM INPUT
///////////////////////////////////////////////////////////////////////////////

char  datafilename[80];
FILE   *stream;

const gsl_rng_type * TT; /*!< Random Number Generator seed based on Gaussian Distribution */
gsl_rng * rr;

int N_bath; /*!< Size of bath */
int N_slice; /*!< Number of time intervals */
int Ncut; /*!< Truncation parameter */
double timestep; /*!< Size of time interval */
int Nsample; /*!< Sample Size (No. of trees calculated) */
double beta; /*!< Inverse Temperature */
double delta; /*!< MD Integrating Timestep*/
double ppower; /*!< */

double ddd4; /*!< */
double ddd; /*!< */
double *m; /*!< Mass of particles */
double *c; /*!< */
double *w; /*!< */
double *f; /*!< Force on particles */

double *dgam; /*!< */
double *mww; /*!< */
double *sig; /*!< Sigma/Variance */
double *mu; /*!< */
double TSLICE; /*!< */

double *abszsum1;
double *argzsum1;
double *habszsum1;
double *hargzsum1;

double (*dens_init[4])(double*, double*); /*!< Initial Density Matrix*/
double (*obs[4])(double*, double*); /*!< Observable Matrix*/
double (*obs1[4])(double*, double*); /*!< Another Observable Matrix*/
void (*force[4])(double *); /*!< Hellman-Feynman Forces*/


// ================================================================
// MAIN
// ================================================================


int main(int argc, char *argv[]){

    double  *R1,  *v;
    double w_max, eta, T;
    int  i,init_seed;

    int Nblock = 1024;
    int t_strobe = 1;

    /*! Sets up the use of a Gaussian Random Number Generator from GSL */
    gsl_rng_env_setup();
    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);


    ///////////////////////////////////////////////////////////////////////////////
    /// SYSTEM INPUT
    ///////////////////////////////////////////////////////////////////////////////
    //cout << " Print information about new stream: " << endl;
    //cout << "Input datafilename, N_bath, N_slice, Ncut\n timestep, T, init_seed, Nsample\n w_max, eta, beta, delta,power " << endl;
    //cin >> datafilename >> N_bath >> N_slice >> Ncut >> timestep >> T >> init_seed >> Nsample >> w_max >> eta >> beta >> delta >> ppower;
    cout << " Print information about new stream:" << endl;
    cout << "Input datafilename" << endl;
    cin >> datafilename;
    N_bath = 200;
    N_slice = 60;
    Ncut = 10;
    timestep = 0.05;
    T = 15;
    init_seed = 0;
    Nsample = 10000;
    w_max = 3;
    eta = 0.13;
    beta = 25;
    delta = 0.8;
    ppower = 100000;

    ///////////////////////////////////////////////////////////////////////////////
    /// ALLOCATING MEMORY
    ///////////////////////////////////////////////////////////////////////////////

    mww = new double[N_bath];
    mu = new double[N_bath];
    sig =  new double[2*N_bath];
    dgam = new double[N_bath];
    R1 = new double[N_bath];
    v = new double[N_bath];
    f = new double[N_bath];
    c = new double[N_bath];
    m = new double[N_bath];
    w = new double[N_bath];
    abszsum1  = new double[N_slice];
    argzsum1  = new double[N_slice];
    habszsum1  = new double[N_slice];
    hargzsum1  = new double[N_slice];

    ///////////////////////////////////////////////////////////////////////////////
    /// INITIALIZATION OF SYSTEM
    ///////////////////////////////////////////////////////////////////////////////


    dens_init[0] = dens_init_0; dens_init[1] = dens_init_1;
    dens_init[2] = dens_init_2; dens_init[3] = dens_init_3;
    obs[0] = obs_0; obs[1] = obs_1; obs[2] = obs_2; obs[3] = obs_3;
    obs1[0] = H_0; obs1[1] = H_1; obs1[2] = H_2; obs1[3] = H_3;

    ddd4 = delta*delta*0.25;
    ddd =  delta*delta;
    TSLICE  = T/N_slice;

    bath_para(eta, w_max);       /*!< Defining Bath Parameters */

    for (i = 0; i < N_bath; ++i)
        mu[i] = beta*w[i]*0.5;
    for (i = 0; i < N_bath; ++i){
        sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
        mww[i] = -m[i]*w[i]*w[i];
    }
    for (i = 0; i < N_bath; ++i)
        sig[i+N_bath] = 1.0*sqrt(w[i]/(2.0*tanh(mu[i])));

    /*! Defining force field */
    force[0] = F1;
    force[1] = Fb;
    force[2] = Fb;
    force[3] = F2;
    /*! Defining non-adiabatic coupling matrix */
    setwww();


    stream = fopen(datafilename,"w");
    fprintf(stream,"%s\n w_max %lf eta %lf beta %lf delta %lf killz %lf N_bath %d N_slice %d\n", argv[0], w_max, eta, beta, delta, ppower, N_bath, N_slice);
    fclose(stream);

    /*! Initializing sum1 counters*/
    for (int i = 0; i < N_slice; ++i){
        abszsum1[i] = 0.0;
        argzsum1[i]  = 0.0;
        habszsum1[i] = 0.0;
        hargzsum1[i] = 0.0;
    }

    ///////////////////////////////////////////////////////////////////////////////
    /// PROCESSING TREE
    ///////////////////////////////////////////////////////////////////////////////
    density(R1, v);
/*#pragma omp parallel for num_threads(8)
    {
        for (int i = 0; i < NN_sample; ++i){
            density(x,p);
        /*    if (((i+1) % Nblock) == 0){
                l  = 1.0/(i+1);
                stream = fopen(datafilename,"a");
                for (k = 0; k < N_slice; k++)
                    if ( ((k+1) % t_strobe) == 0) {
                        for (j = 0; j <=(Ncut+1); j++){
                            fprintf(stream,"%d %lf %d %lf %lf %lf %lf  %lf %lf %lf %lf%.6lf\n", i+1, Dt*(k+1), j, (abszsum1[k])*l, (argzsum1[k])*l, realsum[k][j]*l, imagsum[k][j]*l,(habszsum1[k])*l, (hargzsum1[k])*l, hrealsum[k][j]*l, himagsum[k][j]*l, hist[k][j]*l );
                            printf("%d %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %.6lf\n", i+1, Dt*(k+1), j, (abszsum1[k])*l, (argzsum1[k])*l,realsum[k][j]*l ,imagsum[k][j]*l,(habszsum1[k])*l, (hargzsum1[k])*l, hrealsum[k][j]*l, himagsum[k][j]*l, hist[k][j]*l);              }
                    }
                fclose(stream);
            }
        }
    }*/


    stream = fopen(datafilename,"a");
    fprintf(stream,"dt %lf T %lf Nsample %d\n", timestep, T, Nsample);
    for (int i = 0; i < N_slice; ++i){
        //for (int j =0; j < Ncut;++j){
        printf("%lf   %d  %lf %lf  %lf  %lf\n", TSLICE*(i+1),i, (abszsum1[i]/Nsample), (argzsum1[i]/Nsample),  (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample));
        fprintf(stream,"%lf   %lf  %lf  %lf  %lf\n", TSLICE*(i+1), (abszsum1[i]/Nsample), (argzsum1[i]/Nsample),  (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample));
    }
    fclose(stream);

    ///////////////////////////////////////////////////////////////////////////////
    /// DEALLOCATING MEMEORY
    ///////////////////////////////////////////////////////////////////////////////

    delete [] abszsum1; delete [] argzsum1; delete [] habszsum1; delete [] hargzsum1;
    delete [] mww; delete [] mu; delete [] sig;

    delete [] dgam; delete [] R1; delete [] v; delete [] f;
    delete [] c; delete [] m; delete [] w;


    return 0;

}