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
#include <string.h>


#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

/* Prototypes */

void grad_element(double* mu, double* kappa, double* theta, double* Delta, double* X,
                  double Jump , double* temp_grad, double* temp_prob);

void lik_element(double* mu, double* kappa, double* theta, double* Delta, double* X,
                 double Jump , double* temp_grad, double* temp_prob);

int index_(int i, int j, int N);

/* Function that actually computes the gradient.
 The gradient ends up being a big sum of terms.
 We loop through and compute each element.
 */


void compute( int M,  int n,  int k, double* mu, double* kappa,
             double* theta, double* Delta, double* X, double* DX, double* Post,
             double* lik, double* grad , double* hess, short int numArgs){
    
    int m,i,t;
    int temp_hess_size, ind1, ind2; /* Define time and theta dimension to loop over */
    double temp_lik, temp_grad[3], currPost; /* Temporary place to store some variables */
    double temp_prob[2];
    
    
    /* Compute Log-Likelihood */
    if (numArgs==1){
        for (m=0;m<M;m++) {
            for(t=0;t<n;t++){
                for(i=0;i<k;i++){
                    
                    /* Get current posterior value */
                    currPost = Post[m+M*(i+k*t)];
                    
                    /* Compute partial log-likelihood */
                    lik_element(&mu[i],&kappa[i],&theta[i],Delta,&X[m+M*t],DX[m+M*t],&temp_lik,temp_prob);
                    
                    /* Add new element to likelihood */
                    *lik -= temp_lik*currPost;
                    
                }
            }
        }
    }
    
    /* Compute Log-Likelihood and Gradient */
    else if (numArgs==2) {
        for (m=0;m<M;m++) {
            for(t=0;t<n;t++){
                for(i=0;i<k;i++){
                    
                    /* Get current posterior value */
                    currPost = Post[m+M*(i+k*t)];
                    
                    /* Compute partial log-likelihood */
                    lik_element(&mu[i],&kappa[i],&theta[i],Delta,&X[m+M*t],DX[m+M*t], &temp_lik, temp_prob);
                    
                    /* Compute partial gradient. Stores in temp_grad variable */
                    grad_element(&mu[i],&kappa[i],&theta[i],Delta,&X[m+M*t],DX[m+M*t], temp_grad, temp_prob);
                    
                    /* Add new element to likelihood */
                    *lik -= temp_lik*currPost;
                    
                    /* Properly re-assign elements for gradient */
                    grad[i] -= temp_grad[0]*currPost; // mu_i derivative
                    grad[k+i] -= temp_grad[1]*currPost; // kappa_i derivative
                    grad[2*k+i] -= temp_grad[2]*currPost; // theta_i derivative

                }
            }
        }
    }
}



/*---------------------------------
 MEX Function is below
 ---------------------------------*/


/*  grad = Pois_Likelihood_Objective_Merged(mu, kappa, theta, Delta, X, DX , Post)
 in:   mu,kappa = real number parameters for poisson transition probability
 theta = real array of size 1xk containing mean parameters for each state
 Delta = real positive number representing time step for Poisson process
 X = real array of size 1xn containing observed values of Poisson Process at t=0:(n-1)
 DX = real array of size 1xn containing the jump size of X at t=0:n
 Post = real array of size kxn containing the posterior for the HMC at t=0:(n-1)
 
 out: lik = scalar representing negative of log-likelihood of model      
      grad = 1xk array representing negative gradient of portion of the E-step log likelihood
 */

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ){
    int n, k, M;
    double *mu, *kappa, *theta, *Delta, *X, *DX, *Post, *lik, *grad, *hess;
    
    if (nrhs != 7) mexErrMsgTxt("Seven input arguments required.");
    if (nlhs>2) mexErrMsgTxt("There is a maximum of two output variables.");
    
    /* Get Dimension Sizes */
    n=mxGetN(prhs[5]);
    k=mxGetN(prhs[2]);
    M=mxGetM(prhs[4]);
    
    /* Read arguments into proper C variable */
    mu = mxGetPr(prhs[0]);
    kappa = mxGetPr(prhs[1]);
    theta = mxGetPr(prhs[2]);
    Delta = mxGetPr(prhs[3]);
    X = mxGetPr(prhs[4]);
    DX = mxGetPr(prhs[5]);
    Post = mxGetPr(prhs[6]);
    
    plhs[0] = mxCreateDoubleScalar(0);
    lik = mxGetPr(plhs[0]);
    
    plhs[1] = mxCreateDoubleMatrix(3*k, 1, mxREAL);
    grad = mxGetPr(plhs[1]);
    
    /* Compute the gradient */
    compute(M, n, k, mu, kappa, theta, Delta, X, DX, Post, lik, grad, hess, nlhs);
    
    return;
}




/*---------------------------------
 Helper Functions
 
 These functions are used for the
 computation of the log-likelihood,
 the gradient and the hessian.
 ---------------------------------*/


void lik_element(double* mu, double* kappa, double* theta, double* Delta, double* X,
                 double Jump , double* temp_log, double* temp_prob ) {
    
    double p1, p_1; /* Temporary place to store some variables */
    
    /* Compute jump probability for jumps of the + and - Poisson processes */
    p1 = exp(-(*Delta)*(*mu + *kappa*MAX(0,*theta-*X)));
    p_1 = exp(-(*Delta)*(*mu + *kappa*MAX(0,*X-*theta)));
    
    temp_prob[1] = p1;
    temp_prob[2] = p_1;
    
    /* Depending on the sign of the observed jump, plug in appropriate probability */
    if (Jump>0) *temp_log = log((1-p1)*p_1);
    else if (Jump<0) *temp_log = log((1-p_1)*p1);
    else *temp_log = log(p1*p_1);
    
    
}


void grad_element(double* mu, double* kappa, double* theta, double* Delta, double* X,
                  double Jump , double* temp_grad, double* temp_prob) {
    
    double C1; /* Extra Storage */
    
    /* When there is a positive jump */
    if (Jump>0) {
        C1 = (1/(1-temp_prob[1])-1);
        temp_grad[0] = (*Delta)*(C1-1);
        
        if (*X < *theta) {
            temp_grad[1] = (*Delta)*C1*(*theta-*X);
            temp_grad[2] = (*Delta)*C1*(*kappa);
        }
        else {
            temp_grad[1] = -(*Delta)*(*X - *theta);
            temp_grad[2] = -(*Delta)*(*kappa);
        }
    }
    
    /* When there is a negative jump */
    else if (Jump<0){
        C1 = (1/(1-temp_prob[2])-1);
        temp_grad[0] = (*Delta)*(C1-1) ;
        
        if (*X < *theta) {
            temp_grad[1] = (*Delta)*C1*(*theta-*X);
            temp_grad[2] = (*Delta)*C1*(*kappa);
        }
        else {
            temp_grad[1] = -(*Delta)*(*X-*theta);
            temp_grad[2] = -(*Delta)*(*kappa);
        }
    }
    
    /* When there is no jump */
    else {
        temp_grad[1] = -2*(*Delta);
        temp_grad[2] = -(*Delta)*(*theta - *X);
        temp_grad[3] = temp_grad[1]*(*kappa);
    }
}



