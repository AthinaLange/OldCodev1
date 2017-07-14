
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "random.h"
#include   "functions.h"



using namespace std;

#include <gsl/gsl_rng.h>

// VARIABLES ======================================================================

char  datafilename[80];
FILE   *stream;

const gsl_rng_type * TT;
gsl_rng * rr;


int N_bath;
double ddd4;
double ddd;
double delta;
double abs_d;
double timestep;
double Pdotdhat;
double sina;
double cosa;
double de;

double *m;
double *c;
double *w;
double *f;
double *dhat;
double *dgam;
double *mww;
double *sig;
double *RR;
double *PP;
double *SS;


void (*force[4])(double *);
double (*www[2][4][4])();

// =======================================

int  N_slice;
int Nsample;
int Ncut;

double beta;

double *Pperp;
double *mu;
double *dtdtm;
double TSLICE;
double ppower;
double abszsum0;
double *abszsum1;
double argzsum0;
double *argzsum1;
double habszsum0;
double *habszsum1;
double hargzsum0;
double *hargzsum1;
complex<double> I(0,1);

int t_strobe, Nblock = 1024; /* t_strobe is the frequency at which results for slices are printed,
                                 Nblock is the size of sub-ensembles */

double alpha;

double (*phi)(double*, double*);
double (*dens_init[4])(double*, double*);
double (*obs[4])(double*, double*);
double (*obs1[4])(double*, double*);



int  density(double *x,double *p){

    int SS0,SS1,SS2,SS3,NNjmp =0,signPdotdhat;
    double phase0 = 0.0,xx;
    double p0,p1,p2,p3,ap0,ap1,ap2,ap3;
    double dn2;
    complex<double> z = 1.0;
    complex<double> oldz;

    // Initialization of sample

    gauss_init_W(x, p);
    double yy = 4.0*(gsl_rng_uniform (rr));
    if (yy < 1.0)
        SS3 = (SS0 = 0);
    else if (yy < 2.0){
        SS0 = 1;
        SS3 = 2;
    }
    else if (yy < 3.0){
        SS0 = 2;
        SS3 = 1;
    }
    else
        SS3 = (SS0 = 3);
    SS[0] = SS0;
    z = 4.0;
    for (int l = 0; l < N_bath; ++l){
        RR[l] = x[l];
        PP[l] = p[l];
    }
    SS1 = SS0;
    cout << "rand" << yy << endl;
    cout << "surface" << SS1 << endl;
    // ____________________________________________________________________

    for (int l = 0; l < N_slice; ++l) {
        SS0 = SS1;
        phase0 = U(RR, PP, SS0, TSLICE*0.5); // exp(iLd/2) (before jump)
        z *= exp(I * phase0);

        dd(RR); // non-adiabatic coupling matrix
        de = dE(RR); // energy
        alpha = 0.0;
        Pdotdhat = 0;
        for (int i = 0; i < N_bath; ++i) {
            Pdotdhat += PP[i] * dhat[i]; //parall component of dhat to momentum
        }
        alpha = Pdotdhat * abs_d * TSLICE;

        signPdotdhat = (Pdotdhat < 0 ? -1 : 1); // -1 if neg, 1 if pos
        Pdotdhat = fabs(Pdotdhat);
        for (int i = 0; i < N_bath; ++i)
            Pperp[i] = PP[i] - signPdotdhat * Pdotdhat * dhat[i]; // perp component of dhat
        alpha *= 2.0;
        sina = sin(alpha);
        cosa = cos(alpha);

        ap0 = fabs(p0 = ((www[1][SS0][0]() < -7775.0) ? 0.0 : www[0][SS0][0]()));
        ap1 = fabs(p1 = ((www[1][SS0][1]() < -7775.0) ? 0.0 : www[0][SS0][1]()));
        ap2 = fabs(p2 = ((www[1][SS0][2]() < -7775.0) ? 0.0 : www[0][SS0][2]()));
        ap3 = fabs(p3 = ((www[1][SS0][3]() < -7775.0) ? 0.0 : www[0][SS0][3]()));
        dn2 = ap0 + ap1 + ap2 + ap3;
        xx = dn2 * (gsl_rng_uniform(rr));   // choosing matrix elements
        //alpha goes to 0, pdotdhat very small, matrix becomes identiy and prob of jumping goes to 0

        cout << "Prob:" << "ap0" << ap0 <<" ap1"<< ap1 <<" ap2" << ap2 <<" ap3"<< ap3 << endl;
        oldz = z;
        SS2 = SS1;



        if (xx < ap0) {
            SS1 = 0;
            z *= p0 * dn2 / ap0;
            cout << SS1 << endl;
        } else if (xx < ap0 + ap1) {
            SS1 = 1;
            z *= p1 * dn2 / ap1;
            cout << SS1 << endl;
        } else if (xx < ap0 + ap1 + ap2) {
            SS1 = 2;
            z *= p2 * dn2 / ap2;
            cout << SS1 << endl;
        } else {
            SS1 = 3;
            z *= p3 * dn2 / ap3;
            cout << SS1 << endl;
        }



/*        if (SS0 != SS1)
            NNjmp++;
        if (NNjmp > Ncut)
            return 0;
*/
        if (www[1][SS0][SS1]() != 9999.0)
            for (int i = 0; i < N_bath; ++i)
                PP[i] = Pperp[i] + signPdotdhat * www[1][SS0][SS1]() * dhat[i];


        phase0 = U(RR,PP,SS1,TSLICE*0.5); // exp(iLd/2) (after jump)
        z *= exp(I*phase0);

        phi = obs[SS1];
        abszsum0  = real(z*phi(RR,PP)*dens_init[SS3](x,p));
        argzsum0  = imag(z*phi(RR,PP)*dens_init[SS3](x,p));
        abszsum1[l] += abszsum0;
        argzsum1[l] += argzsum0;

        phi = obs1[SS1];
        habszsum0  = real(z*phi(RR,PP)*dens_init[SS3](x,p));
        hargzsum0  = imag(z*phi(RR,PP)*dens_init[SS3](x,p));
        habszsum1[l] += habszsum0;
        hargzsum1[l] += hargzsum0;
        cout << "lth " <<l << endl;
    }

    return 0;
}



int monte(int NN_sample,double *x, double *p){

    for (int i = 0; i < N_slice; ++i){
        abszsum1[i] = 0.0;
        argzsum1[i]  = 0.0;
        habszsum1[i] = 0.0;
        hargzsum1[i] = 0.0;
    }
#pragma omp parallel for num_threads(8)
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
            }*/
        }
    }
    return 0;

}


int main(int argc, char *argv[]){

    double  *R1,  *v;
    double w_max, eta,T;
    int  i,init_seed;

    t_strobe = 1;   /* Print out the results for every t_strobe slice  */

    gsl_rng_env_setup();

    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);

    //printf(" Print information about new stream:\n");
    //fputs("Input datafilename, N_bath, N_slice, Ncut\n timestep, T, init_seed, Nsample\n w_max, eta, beta, delta,power\n ",stderr);
    //scanf("%s%d%d%d%lf%lf%d%d%lf%lf%lf%lf%lf", datafilename, &N_bath, &N_slice, &Ncut, &timestep, &T, &init_seed, &Nsample, &w_max, &eta, &beta, &delta, &ppower);
    /* initialize stream  - scope 0 */
    cout << " Print information about new stream:" << endl;
    cout << "Input datafilename" << endl;
    cin >> datafilename;
    N_bath = 200;
    N_slice = 20;
    Ncut = 10;
    timestep = 0.05;
    T = 15;
    init_seed = 0;
    Nsample = 100;
    w_max = 3;
    eta = 0.13;
    beta = 25;
    delta = 0.8;
    ppower = 100000;

    /* Allocate memory ----------------------------------------------------- */

    mww = new double[N_bath];
    mu = new double[N_bath];
    sig =  new double[2*N_bath];
    dtdtm = new double[N_bath];
    dgam = new double[N_bath];
    dhat = new double[N_bath];
    R1 = new double[N_bath];
    v = new double[N_bath];
    f = new double[N_bath];
    c = new double[N_bath];
    m = new double[N_bath];
    w = new double[N_bath];
    RR = new double[N_bath];
    PP = new double[N_bath];
    SS = new double[N_slice];
    Pperp = new double[N_bath];
    abszsum1  = new double[N_slice];
    argzsum1  = new double[N_slice];
    habszsum1  = new double[N_slice];
    hargzsum1  = new double[N_slice];

    dens_init[0] = dens_init_0; dens_init[1] = dens_init_1;
    dens_init[2] = dens_init_2; dens_init[3] = dens_init_3;
    obs[0] = obs_0; obs[1] = obs_1; obs[2] = obs_2; obs[3] = obs_3;
    obs1[0] = H_0; obs1[1] = H_1; obs1[2] = H_2; obs1[3] = H_3;

    ddd4 = delta*delta*0.25;
    ddd =  delta*delta;
    TSLICE  = T/N_slice;

    bath_para(eta,w_max);       /* compute system parameters etc */

    //  bath corresponds to eq. 53
    for (i = 0; i < N_bath; i++)
        mu[i] = beta*w[i]*0.5;
    for (i = 0; i < N_bath; i++){
        sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
        mww[i] = -m[i]*w[i]*w[i];
        dtdtm[i] = -0.5*timestep*timestep/m[i];
    }

    for (i = 0; i < N_bath; i++)
        sig[i+N_bath] = 1.0*sqrt(w[i]/(2.0*tanh(mu[i])));
    force[0] = F1;         /* assign pointers to force fields */
    force[1] = Fb;
    force[2] = Fb;
    force[3] = F2;
    setwww();


/*    stream = fopen(datafilename,"w");
    fprintf(stream,"%s\n w_max %lf eta %lf beta %lf delta %lf killz %lf N_bath %d N_slice %d\n", argv[0], w_max, eta, beta, delta, ppower, N_bath, N_slice);
    fprintf(stream,"Nens\t time\t j\t O_j\t O_tot\t En_j\t En_tot\n");
    fclose(stream);*/

    monte(Nsample,R1,v);

/*    stream = fopen(datafilename,"a");
    fprintf(stream,"dt %lf T %lf Nsample\n", timestep, T, Nsample);
    for (i = 0; i < N_slice; i++)
        if (((i+1)% t_strobe) == 0)
            for (j =0; j <= (Ncut+1);j++){
                printf("%lf   %lf  %lf  %lf  %lf  %lf  %lf    %.6lf\n", Dt*(i+1), (abszsum1[i]/Nsample), (argzsum1[i]/Nsample), realsum[i][j]/Nsample,  (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample), hrealsum[i][j]/Nsample,  1.0*hist[i][j]/Nsample);
                fprintf(stream,"%lf   %lf  %lf  %lf  %lf  %lf  %lf %.6lf\n", Dt*(i+1), (abszsum1[i]/Nsample), (argzsum1[i]/Nsample), realsum[i][j]/Nsample, (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample), hrealsum[i][j]/Nsample, 1.0*hist[i][j]/Nsample);


            }
    fclose(stream);*/


    delete [] abszsum1; delete [] argzsum1; delete [] habszsum1; delete [] hargzsum1;
    delete [] Pperp; delete [] mww; delete [] mu; delete [] sig; delete [] dtdtm;
    delete [] dgam; delete [] dhat; delete [] R1; delete [] v; delete [] f;
    delete [] c; delete [] m; delete [] w; delete [] RR; delete [] PP; delete [] SS;

    return 0;

}