/*
 *  Based on files written by Aurelien Garivier, CNRS & Telecom Paristech
 *  January 2012
 *
 * Baum-Welch algorithm for discrete Hidden Markov models
 * see http://www.telecom-paristech.fr/~garivier/code/index.html
 * 
 * *=================================================================
 *
 * compile with "make", or:
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

void filter(int N, int n, double* nu, double* Q, double* g, double* phi, double* c){
	int i,j,t;
	double z[N];
	c[0]=0;
	for(j=0; j<N; j++){
	  z[j] = nu[j]*g[j]; /* g(:,1) in matlab */
		c[0] += z[j];
	}
	for(j=0; j<N; j++) phi[j+0*N] = z[j]/c[0];
	for(t=1; t<n; t++){
		c[t]=0;
		for(j=0; j<N; j++){
			z[j]=0;
			for(i=0; i<N; i++) z[j]+=phi[i+(t-1)*N]*Q[i+N*j]*g[j+N*t];
			c[t] += z[j];
		}
		for(j=0; j<N; j++) phi[j+N*t] = z[j]/c[t];
	}
}

/*  [phi, c] = filter(nu, Q, g)
in:   nu = initial distribution as vector of size k
      Q = transition matrix of size k
      g = pre-computed emission probability of size k x n
out:  phi = filter(x,t) = P(X(t)=x | Y(1:t)=y(1:t)) for 1<=x<=k and 1<=t<=n
c(t) = P(Y(t) = y(t)| Y(1:t-1)=y(1:t-1))*/

void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[] ){
  int N,n;
  double *nu, *Q, *g, *phi, *c;

  if (nrhs != 3) mexErrMsgTxt("Four input arguments required."); 
  else if (nlhs > 2) mexErrMsgTxt("Two output arguments provided."); 

/* Read arguments into proper C variable */
  N = mxGetN(prhs[0]);
  nu = mxGetPr(prhs[0]);
  Q = mxGetPr(prhs[1]);
  g = mxGetPr(prhs[2]);
  n = mxGetN(prhs[2]);
  
  plhs[0] = mxCreateDoubleMatrix(N, n, mxREAL);
  phi = mxGetPr(plhs[0]);
  plhs[1] = mxCreateDoubleMatrix(1, n, mxREAL);
  c = mxGetPr(plhs[1]);

  filter(N, n, nu, Q, g, phi, c);
 
  return;    
}
