#include   <stdlib.h>
#include   <stdio.h>
#include   <math.h>
#include   <iostream>
#include   <complex>
#include   "dmatrix.c"
#include   "imatrix.c"
#include   "omp.h"
#include   "variable.h"
#include   "variable-trajectory.h"
#include   "density.cpp"
#include   "functions.cpp"
#include   "transition.cpp"
using namespace std;

#include <gsl/gsl_rng.h>

char  datafilename[80];
FILE   *stream;

/* For the definition of gamma, dgamma etc, please
   see our JCP paper  */


/* Monte Carlo Sampling ----------------------------------------------------- */

int t_strobe, Nblock = 1024; /* t_strobe is the frequency at which results for slices are printed,
                                 Nblock is the size of sub-ensembles */


int monte(int NN_sample,double *x, double *p){

    int i, j, k, flag, sl;
    double y, l;

    /* This section calls slice algorithm, and prints results  */
    for (i = 0; i < N_slice; i++){
        abszsum1[i]  = 0.0;
        argzsum1[i]  = 0.0;
        habszsum1[i] = 0.0;
        hargzsum1[i] = 0.0;
    }
#pragma omp parallel for num_threads(8)
    {
        for (i = 0; i < NN_sample; i++){
            density(x,p);
            if (((i+1) % Nblock) == 0){
                l  = 1.0/(i+1);
                stream = fopen(datafilename,"a");
                for (k = 0; k < N_slice; k++)
                    if (((k+1) % t_strobe) == 0) {
                        for (j = 0; j <=(Ncut+1); j++){
                            fprintf(stream,"%d %lf %d %lf %lf %lf %lf  %lf %lf %lf %lf%.6lf\n", i+1, Dt*(k+1), j, (abszsum1[k])*l, (argzsum1[k])*l, realsum[k][j]*l, imagsum[k][j]*l,(habszsum1[k])*l, (hargzsum1[k])*l, hrealsum[k][j]*l, himagsum[k][j]*l, hist[k][j]*l );
                            printf("%d %lf %d %lf %lf %lf %lf %lf %lf %lf %lf %.6lf\n", i+1, Dt*(k+1), j, (abszsum1[k])*l, (argzsum1[k])*l,realsum[k][j]*l ,imagsum[k][j]*l,(habszsum1[k])*l, (hargzsum1[k])*l, hrealsum[k][j]*l, himagsum[k][j]*l, hist[k][j]*l);              }
                    }
                fclose(stream);
            }
        }
    }
    return 0;
}

/* Main ------------------------------------------------- */


int MAX_tr;
int *jump;
double *b, beta, *mu, *dtdtm;
//int  rngstreamnum, nrngstreams, *rngstream;
//int rng_seed;
//int rho_top_init[4], obs_top[4];
//double BETA1,BETA2, GAMMA1, GAMMA2;

int main(int argc, char *argv[]){

    double z, sum1, sum2, z2, sum3, sum5, xx, Dt;
    double tw[4][4];
    double  dt, t1, *R1, *v, x;
    double w_max, eta, reS, imS, prod, Ndev, yy, T;
    int  i, Sa, j, Nsteps, N_space, ii, k, init_seed, ppp, N0, stepnum, error_flag = 0;
    int i1, i2;
    t_strobe = 1;   /* Print out the results for every t_strobe slice  */

    gsl_rng_env_setup();

    TT = gsl_rng_default;
    rr = gsl_rng_alloc(TT);

    // rngstreamnum = 0;
    // nrngstreams = 1;
    // rng_seed = make_sprng_seed();
    // rngstream = init_sprng(rngstreamnum,nrngstreams,rng_seed,SPRNG_DEFAULT);
    /* initialize stream   */
    printf(" Print information about new stream:\n");
    // print_sprng(rngstream);
    fputs("Input datafilename, N_bath, N_slice, Ncut\n timestep, T, init_seed, Nsample\n w_max, eta, beta, delta,power\n ",stderr);
    scanf("%s%d%d%d%lf%lf%d%d%lf%lf%lf%lf%lf", datafilename, &N_bath, &N_slice, &Ncut, &timestep, &T, &init_seed, &Nsample, &w_max, &eta, &beta, &delta, &ppower);
    reS = 0.0; imS = 0.0;
    MAX_tr = 2;


    /* Allocate memory ----------------------------------------------------- */
    R1 = (double *)malloc((N_bath)*sizeof(double));
    v =  (double *)malloc((N_bath)*sizeof(double));
    f =  (double *)malloc((N_bath)*sizeof(double));
    m =  (double *)malloc((N_bath)*sizeof(double));
    w =  (double *)malloc((N_bath)*sizeof(double));
    c =  (double *)malloc((N_bath)*sizeof(double));
    dhat =  (double *)malloc((N_bath)*sizeof(double));
    dgam =  (double *)malloc((N_bath)*sizeof(double));
    b  =  (double *)malloc((N_bath)*sizeof(double));
    Pperp  =  (double *)malloc((N_bath)*sizeof(double));
    SS = (int *)malloc((2)*sizeof(double));
    RR = dmatrix(0,MAX_tr-1,0,N_bath-1);
    PP = dmatrix(0,MAX_tr-1,0,N_bath-1);
    realsum = dmatrix(0,N_slice,0,N_slice+1);
    imagsum = dmatrix(0,N_slice,0,N_slice+1);
    hrealsum = dmatrix(0,N_slice,0,N_slice+1);
    himagsum = dmatrix(0,N_slice,0,N_slice+1);

    hist = imatrix(N_slice,N_slice+1);

    mu =  (double *)malloc((N_bath)*sizeof(double));
    mww = (double *)malloc((N_bath)*sizeof(double));
    dtdtm = (double *)malloc((N_bath)*sizeof(double));
    sig =  (double *)malloc((2*N_bath)*sizeof(double));
    meann = (double *)malloc((2*N_bath)*sizeof(double));
    abszsum0  = (double  *)malloc((N_slice)*sizeof(double));
    abszsum1  = (double  *)malloc((N_slice)*sizeof(double));
    argzsum0  = (double  *)malloc((N_slice)*sizeof(double));
    argzsum1  = (double  *)malloc((N_slice)*sizeof(double));

    habszsum0  = (double  *)malloc((N_slice)*sizeof(double));
    habszsum1  = (double  *)malloc((N_slice)*sizeof(double));
    hargzsum0  = (double  *)malloc((N_slice)*sizeof(double));
    hargzsum1  = (double  *)malloc((N_slice)*sizeof(double));

    dens_init[0] = dens_init_0; dens_init[1] = dens_init_1;
    dens_init[2] = dens_init_2; dens_init[3] = dens_init_3;
    obs[0] = obs_0; obs[1] = obs_1; obs[2] = obs_2; obs[3] = obs_3;
    obs1[0] = H_0; obs1[1] = H_1; obs1[2] = H_2; obs1[3] = H_3;
    // obs[0] = H_0; obs[1] = H_1; obs[2] = H_2; obs[3] = H_3;
    ddd4 = delta*delta*0.25;
    ddd =  delta*delta;
    TSLICE  = T/N_slice;
    Dt = TSLICE;
    Ndev = 4.0;
    bath_para(eta,w_max);       /* compute system parameters etc */
    /*
    bath corresponds to eq. 54 -- asymmetric boson
    for (i = 0; i < N_bath; i++)
       mu[i] = beta*w[i]/2.0;
    for (i = 0; i < N_bath; i++){
       sig[i] = 1.0/sqrt(w[i]*2.0*tanh(mu[i]));
       meann[i] = c[i]/(w[i]*w[i]);
       mww[i] = -m[i]*w[i]*w[i];
       dtdtm[i] = -0.5*timestep*timestep/m[i];
    }

    for (i = 0; i < N_bath; i++){
       sig[i+N_bath] = 1.0*sqrt(w[i]/(2.0*tanh(mu[i])));
       meann[i+N_bath] = 0.0;
    }
    */
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
    /*
    alpha = 0.0; BETA = 1.5*alpha; GAMMA = 2.0*alpha;


      for (i1 = 0; i1 < 4 ;i1++)
	for(i2 = 0; i2 < 4; i2++){
           tw[i1][i2] = 0.0;
	   for (j=0;j <3;j++)
	      tw[i1][i2] +=  www[j][i1][i2]();
        }
     for (i1 = 0; i1 < 4 ;i1++){
        printf("\n");
        for(i2 = 0; i2 < 4; i2++)
	    printf("%lf\t", tw[i1][i2]);
        printf("\n");
    }
     printf("OK\n");

    Checked for alpha = beta = gamma = 0 and gives unity as it should

    */


    stream = fopen(datafilename,"w");
    fprintf(stream,"%s\n w_max %lf eta %lf beta %lf delta %lf killz %lf N_bath %d N_slice %d\n", argv[0], w_max, eta, beta, delta, ppower, N_bath, N_slice);
    fprintf(stream,"Nens\t time\t j\t O_j\t O_tot\t En_j\t En_tot\n");
    fclose(stream);
    for (i = 0; i < N_slice; i++)
        for (j = 0; j <= N_slice; j++){
            realsum[i][j] = 0.0;
            imagsum[i][j] = 0.0;
            hrealsum[i][j] = 0.0;
            himagsum[i][j] = 0.0;

            hist[i][j] = 0;
        }
    monte(Nsample,R1,v);
    stream = fopen(datafilename,"a");
    fprintf(stream,"dt %lf T %lf Nsample %d\n", timestep, T, Nsample);
    for (i = 0; i < N_slice; i++)
        if (((i+1)% t_strobe) == 0)
            for (j =0; j <= (Ncut+1);j++){
                printf("%lf   %lf  %lf  %lf  %lf  %lf  %lf    %.6lf\n", Dt*(i+1), (abszsum1[i]/Nsample), (argzsum1[i]/Nsample), realsum[i][j]/Nsample,  (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample), hrealsum[i][j]/Nsample,  1.0*hist[i][j]/Nsample);
                fprintf(stream,"%lf   %lf  %lf  %lf  %lf  %lf  %lf %.6lf\n", Dt*(i+1), (abszsum1[i]/Nsample), (argzsum1[i]/Nsample), realsum[i][j]/Nsample, (habszsum1[i]/Nsample), (hargzsum1[i]/Nsample), hrealsum[i][j]/Nsample, 1.0*hist[i][j]/Nsample);
            }
    fclose(stream);

    free(R1);  free(v); free(f); free(b);
    free(m), free(c); free(w); free(dhat); free(dgam);
    free(SS); free(jump); free(Pperp);
    free_dmatrix(RR,0,MAX_tr-1,0,N_bath-1);
    free_dmatrix(PP,0,MAX_tr-1,0,N_bath-1);
    free_dmatrix(realsum,0,N_slice,0,N_slice+1);
    free_dmatrix(imagsum,0,N_slice,0,N_slice+1);
    free_dmatrix(hrealsum,0,N_slice,0,N_slice+1);
    free_dmatrix(himagsum,0,N_slice,0,N_slice+1);
    free_imatrix(hist,N_slice,N_slice+1);
    free(abszsum0);  free(argzsum0);    free(abszsum1);  free(argzsum1);
    free(habszsum0);  free(hargzsum0);    free(habszsum1);free(hargzsum1);
    free(mu);free(sig); free(mww);free(dtdtm); free(meann);

    return 0;
}