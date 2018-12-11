/*
 * Written by Philippe Casgrain on Dec 23, 2016
 *
 * e-mail: p.casgrain@mail.utoronto.ca
 *
 * *=================================================================
 *
 * compile with
 * mex -v  HMMfilter.c
 * optimized: mex COPTIMFLAGS='-O2' -v HMMfilter.c
 *
 *=================================================================*/

#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif


/* Actual function that computes the transition probs */

void compute_transition(int M, int n, int k, double* mu, double* kappa, double* theta,
                        double Delta, double* X, double* DX, double* PXTrans){
    
    int m,i,t; /* Define time and theta dimension to loop over */
    double p1, p_1; /* Define some temporary storage */
    
    for (m=0;m<M;m++) {
        for (t=0;t<n;t++) {
            for (i=0;i<k;i++){
                /* Compute jump probability for jumps of the + and - Poisson processes */
                p1 = 1 - exp(-Delta*(mu[i] + kappa[i]*MAX(0,theta[i]-X[m+M*t])));
                p_1 =  1 - exp(-Delta*(mu[i] + kappa[i]*MAX(0,X[m+M*t]-theta[i])));
                
                /* Depending on the sign of the observed jump, plug in appropriate probability */
                if (DX[m+M*t]>0) PXTrans[m+M*(i+k*t)] = p1*(1-p_1);
                else if (DX[m+M*t]<0) PXTrans[m+M*(i+k*t)] = p_1*(1-p1);
                else PXTrans[m+M*(i+k*t)] = (1-p1)*(1-p_1);
            }
        }
    }
}


/*---------------------------------
    MEX Function stuff is below
 ---------------------------------*/


/*  grad = Pois_Transition_Prob(mu, kappa, theta, Delta, X, DX)
 in:   mu,kappa = real number parameters for poisson transition probability
 theta = real array of size 1xk containing mean parameters for each state
 Delta = real positive number representing time step for Poisson process
 X = real array of size 1xn containing observed values of Poisson Process at (m,t)=(0:M-1,0:n-1)
 DX = real array of size 1xn containing the jump size of X at (m,t)=(0:M-1,1:n)
 
 out:  PXTrans = real array of size kxn containing the posterior for the HMC at t=0:(n-1)
 */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ){
    int n, k, M , dims[3];
    double *mu, *kappa, *theta, *Delta, *X, *DX, *PXTrans;
    
    if (nrhs != 6) mexErrMsgTxt("Six input arguments required.");
    
    /* Read arguments into proper C variables */
    n=mxGetN(prhs[4]);
    k=mxGetN(prhs[2]);
    M = mxGetM(prhs[4]);
    
    mu = mxGetPr(prhs[0]);
    kappa = mxGetPr(prhs[1]);
    theta = mxGetPr(prhs[2]);
    Delta = mxGetPr(prhs[3]);
    X = mxGetPr(prhs[4]);
    DX = mxGetPr(prhs[5]);
    
    /* Fill array with output Dimensions */
    dims[0]=M; dims[1]=k; dims[2]=n;
    
    /* Create an output matrix */
    plhs[0] = mxCreateNumericArray(3,dims,mxDOUBLE_CLASS,mxREAL);
    PXTrans = mxGetPr(plhs[0]);
    
    
    /* Compute the gradient */
    compute_transition(M, n, k, mu, kappa, theta, *Delta, X, DX, PXTrans);
    
    return;    
}




