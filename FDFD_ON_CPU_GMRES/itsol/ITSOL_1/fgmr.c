#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "./LIB/globheads.h"
#include "./LIB/protos.h"


#define  epsmac  1.0e-16

int fgmr(SMatptr Amat, SPreptr lu, double *rhs, double *sol, double tol,
	    int im, int *itmax, FILE *fits){ 
/*----------------------------------------------------------------------
|                 *** Preconditioned FGMRES ***                  
+-----------------------------------------------------------------------
| This is a simple version of the ARMS preconditioned FGMRES algorithm. 
+-----------------------------------------------------------------------
| Y. S. Dec. 2000. -- Apr. 2008  
+-----------------------------------------------------------------------
| on entry:
|---------- 
|
|(Amat)   = matrix struct. the matvec operation is Amat->matvec.
|(lu)     = preconditioner struct.. the preconditioner is lu->precon
|           if (lu == NULL) the no-preconditioning option is invoked.
| rhs     = real vector of length n containing the right hand side.
| sol     = real vector of length n containing an initial guess to the
|           solution on input.
| tol     = tolerance for stopping iteration
| im      = Krylov subspace dimension 
| (itmax) = max number of iterations allowed. 
| fits    = NULL: no output
|        != NULL: file handle to output " resid vs time and its" 
|
| on return:
|---------- 
| fgmr      int =  0 --> successful return.
|           int =  1 --> convergence not achieved in itmax iterations.
| sol     = contains an approximate solution (upon successful return).
| itmax   = has changed. It now contains the number of steps required
|           to converge -- 
+-----------------------------------------------------------------------
| internal work arrays:
|----------       
| vv      = work array of length [im+1][n] (used to store the Arnoldi
|           basis)
| hh      = work array of length [im][im+1] (Householder matrix)
| z       = work array of length [im][n] to store preconditioned vectors
+-----------------------------------------------------------------------
| subroutines called :
| lu->precon - preconditionning operation 
+---------------------------------------------------------------------*/
  int n=Amat->n, maxits = *itmax; 
  int i, i1, ii, j, k, k1, its, retval, tmp1 = 1;
  double **hh, *c, *s, *rs, t;
  double negt, beta, eps1=0.0, gam, **vv, **z; 
  double t1=0.0, t2=0.0;
  its = 0;
  vv = (double **)Malloc((im+1)*sizeof(double *), "fgmres");
  z  = (double **)Malloc(im*sizeof(double *), "fgmres");
  for (i=0; i<=im; i++) {
    vv[i] = (double *)Malloc(n*sizeof(double), "fgmres");
  }
  hh = (double **)Malloc(im*sizeof(double *), "fgmres");
  for (i=0; i<im; i++) {
    hh[i] = (double *)Malloc((i+2)*sizeof(double), "fgmres");
    z[i]  = (double *)Malloc(n*sizeof(double), "fgmres");
  }
  c  = (double *)Malloc(im*sizeof(double), "fgmres");
  s  = (double *)Malloc(im*sizeof(double), "fgmres");
  rs = (double *)Malloc((im+1)*sizeof(double), "fgmres");
/*-------------------- outer loop starts here */
label20:
/*-------------------- compute initial residual vector */
  Amat->matvec(Amat, sol, vv[0]); 
  for (j=0; j<n; j++)
    vv[0][j] = rhs[j] - vv[0][j];    /*  vv[0]= initial residual */
  beta = DNRM2(n, vv[0], tmp1);
/*-------------------- print info if fits != null */
  if (fits != NULL && its == 0) {
    t1 = sys_timer();
//    fprintf(fits,"%8d   %10.2e   %10.2e\n",its,t1-t1, beta) ;     
    fprintf(fits,"%8d   %10.2e\n",its, beta) ;     
  }
  if (beta == 0.0) goto label990;
  t = 1.0 / beta;
/*--------------------   normalize:  vv[0] = vv[0] / beta */
  for (j=0; j<n; j++)
    vv[0][j] = vv[0][j]*t;
  if (its == 0) eps1 = tol*beta;
/*    ** initialize 1-st term  of rhs of hessenberg system   */
  rs[0] = beta;
  i = -1;
 label4:
  i++;
  its++;
  i1 = i+1; 
/*------------------------------------------------------------
|  (Right) Preconditioning Operation   z_{j} = M^{-1} v_{j}
+-----------------------------------------------------------*/
  if (lu == NULL) 
    memcpy(z[i],vv[i],n*sizeof(double));
  else{
  /*-------------------- ## debug for no precon, just comment next line! */
    lu->precon(vv[i], z[i],lu) ;
    }
/*-------------------- matvec operation w = A z_{j} = A M^{-1} v_{j} */
   Amat->matvec(Amat, z[i], vv[i1]); 
/*------------------------------------------------------------
|     modified gram - schmidt...
|     h_{i,j} = (w,v_{i})
|     w  = w - h_{i,j} v_{i}
+------------------------------------------------------------*/
   for (j=0; j<=i; j++) {
     t = DDOT(n, vv[j], tmp1, vv[i1], tmp1);
     hh[i][j] = t;
     negt = -t;
     DAXPY(n, negt, vv[j], tmp1, vv[i1], tmp1);
   }
/*-------------------- h_{j+1,j} = ||w||_{2}    */
   t = DNRM2(n, vv[i1], tmp1);
   hh[i][i1] = t;
   if (t == 0.0) goto label58;
   t = 1.0/t;
/*-------------------- v_{j+1} = w / h_{j+1,j}  */
   for (k=0; k<n; k++)
     vv[i1][k] = vv[i1][k]*t;
/*---------------------------------------------------
|     done with modified gram schimdt and arnoldi step
|     now  update factorization of hh
+--------------------------------------------------*/
 label58:
/*--------------------------------------------------------
|   perform previous transformations  on i-th column of h
+-------------------------------------------------------*/
   for (k=1; k<=i; k++) {
     k1 = k-1;
     t = hh[i][k1];
     hh[i][k1] = c[k1]*t + s[k1]*hh[i][k];
     hh[i][k] = -s[k1]*t + c[k1]*hh[i][k];
   }
   gam = sqrt( pow(hh[i][i],2) + pow(hh[i][i1],2) );
/*---------------------------------------------------
|     if gamma is zero then any small value will do...
|     will affect only residual estimate
+--------------------------------------------------*/
   if (gam == 0.0) gam = epsmac;
/*-------------------- get  next plane rotation    */
   c[i] = hh[i][i]/gam;
   s[i] = hh[i][i1]/gam;
   rs[i1] = -s[i]*rs[i];
   rs[i] =  c[i]*rs[i];
/*----------------------------------------------------
|   determine residual norm and test for convergence
+---------------------------------------------------*/
   hh[i][i] = c[i]*hh[i][i] + s[i]*hh[i][i1];
   beta = fabs(rs[i1]);
   if( fits != NULL ) {
     t2 = sys_timer();
//     fprintf(fits,"%8d   %10.2e   %10.2e\n",its,t2-t1, beta) ; 
    fprintf(fits,"%8d   %10.2e\n",its, beta) ;          
   }
   if ( (i < im-1) && (beta > eps1) && (its < maxits) )  goto label4;
/*-------------------- now compute solution. 1st, solve upper 
                       triangular system*/
   rs[i] = rs[i]/hh[i][i];
   for (ii=1; ii<=i; ii++) {
     k=i-ii;
     k1 = k+1;
     t=rs[k];
     for (j=k1; j<=i; j++)
       t = t - hh[j][k]*rs[j];
     rs[k] = t/hh[k][k];
   }
/*--------------------  linear combination of v[i]'s to get sol. */
   for (j=0; j<=i; j++) {
     t = rs[j];
     for (k=0; k<n; k++)
       sol[k] += t*z[j][k];
   }
/*--------------------  restart outer loop if needed */
   if (beta <= eps1) goto label990;
   if (its >= maxits) goto label991;
   goto label20;
label990:
   retval = 0;
   goto label888;
label991:
   retval = 1;
label888:
   for (i=0; i<=im; i++)
     free(vv[i]);
   free(vv);
   for (i=0; i<im; i++){
     free(hh[i]);
     free(z[i]);
   }
   free(hh);
   free(z);
   free(c);
   free(s);
   free(rs);
   *itmax = its; 
   return retval;
}
/*-----------------end of fgmr ---------------------------------------
----------------------------------------------------------------------*/
