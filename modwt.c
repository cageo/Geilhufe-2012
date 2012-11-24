#include <stdio.h>
#include "math.h"
#include "mex.h"

int mod_1_N(int number, int divisor) {                 /* GMA */
  /* Compute positive number between 1 and divisor, inclusive */
  int Mod = number%divisor;
  Mod = Mod < 0 ? Mod+divisor : Mod;
  return (Mod);
}

void modwt(double *Vin, int N, int k, int *L, double *ht, double *gt, 
	   double *Wout, double *Vout)
{

  int j, l, t;

  for(t = 0; t < N; t++) {
    j = t;
    Wout[t] = ht[0] * Vin[j];
    Vout[t] = gt[0] * Vin[j];
    for(l = 1; l < *L; l++) {
      j -= pow(2.0, k - 1.0);
      j = mod_1_N(j, N);                 /* GMA bug fix */
      Wout[t] += ht[l] * Vin[j];
      Vout[t] += gt[l] * Vin[j];
    }
  }

}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  int M, J, L;
  int m, n;
  double *Wout, *Vout;
  double *Vin, *ht, *gt;
  
  /* Check for proper number of arguments */
  
  if (nrhs != 4) {
    mexErrMsgTxt("DWT requires four input arguments.");
  } else if (nlhs > 2) {
    mexErrMsgTxt("DWT requires two output arguments.");
  }
  
  Vin = mxGetPr(prhs[0]);
  ht = mxGetPr(prhs[1]);
  gt = mxGetPr(prhs[2]);
  J = (int) mxGetScalar(prhs[3]);
  /* mexPrintf("J = %d\n", J); */

  m = mxGetM(prhs[0]);
  /* mexPrintf("m = %d\n", m); */
  n = mxGetN(prhs[0]);
  /* mexPrintf("n = %d\n", n); */
  
  /* Create matrices for the return arguments */
  
  plhs[0] = mxCreateDoubleMatrix(m, n, mxREAL);
  plhs[1] = mxCreateDoubleMatrix(m, n, mxREAL);
  
  /* Assign pointers to the various parameters */
  
  Wout = mxGetPr(plhs[0]);
  Vout = mxGetPr(plhs[1]);
  
  M = mxGetNumberOfElements(prhs[0]);
  L = mxGetNumberOfElements(prhs[1]);
  
  /* Do the actual computations in a subroutine */
  
  modwt(Vin, M, J, &L, ht, gt, Wout, Vout);
  return;

}
