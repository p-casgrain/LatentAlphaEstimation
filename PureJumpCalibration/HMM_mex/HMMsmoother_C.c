/*  
 *  Based on files written by Aurelien Garivier, CNRS & Telecom Paristech
 *  January 2012
 *
 * Baum-Welch algorithm for discrete Hidden Markov models
 * see http://www.telecom-paristech.fr/~garivier/code/index.html
 * 
 * *=================================================================
 *
 * 
 * compile with "make", or:
 * mex -v  HMMsmoother.c 
 * optimized: mex COPTIMFLAGS='-O2' -v HMMsmoother.c 
 *
 *=================================================================*/

#include <math.h>
#include "mex.h"
#include <stdlib.h>
#include <stdio.h>

#if !defined(MAX)
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif

void smoother(int N, int n, double* Q, double* g, double* c, double* beta){
	int i,j,t;
	for(j=0; j<N; j++) beta[j+N*(n-1)] = 1;
	for(t=n-2; t>=0; t--){
		for(i=0; i<N; i++){
			double z=0;
			for(j=0; j<N; j++)
			  z += Q[i+N*j] * g[j+N*(t+1)] * beta[j+N*(t+1)];
			beta[i+N*t] = z / c[t+1];
		}
	}
}

/*function beta = smoother(Q, g, c)
in:  Q = transition matrix of size k
     g = pre-computed emission probability of size k x n
     c = vector computed by function filter: c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))
out: beta = smoothing factors: 
     beta(x,t) = P(Y(t+1:n)=y(t+1:n) | X(t)=x) / P(Y(t+1:n)=y(t+1:n) | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
     permits to compute the posterior distribution of the hidden states as post = phi .* beta:
     post(x,t) = P(X(t)=x | Y(1:n)=y(1:n))
*/
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ){
  int N,n;
  double *Q, *g, *c, *beta;
  if (nrhs != 3) mexErrMsgTxt("Four input arguments required."); 
  else if (nlhs > 1) mexErrMsgTxt("One output arguments provided."); 

/* Read arguments into proper C variable */
  N = mxGetN(prhs[0]);
  Q = mxGetPr(prhs[0]);
  g = mxGetPr(prhs[1]);
  n = mxGetN(prhs[1]);
  c = mxGetPr(prhs[2]);


  plhs[0] = mxCreateDoubleMatrix(N, n, mxREAL);
  beta = mxGetPr(plhs[0]);

  smoother(N, n, Q, g, c, beta);
  return;    
}
