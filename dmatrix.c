/* > c.dmatrix */ 


#include <stdlib.h>



double **dmatrix(int nrl,int nrh,int ncl,int nch) { 

   int i; 

   double **m; 

   m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*)); 

   m -= nrl; 

   for(i=nrl;i<=nrh;i++) { 

      m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double)); 

      m[i] -= ncl; 

   } 

   return m; 

}


void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch) { 

   int i;

   for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl)); 

   free((char*) (m+nrl)); 

} 
