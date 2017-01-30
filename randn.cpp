/*
 * For unit testing with ARMADILLO random matrices
 * it is nice to compare with MATLAB code with some existing
 * test functions.
 *
 *
 */

#include <armadillo>
#include "armaMex.hpp"

using namespace arma;

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check the number of input arguments.
  if (nrhs > 2)
    mexErrMsgTxt("Currently supporting only random scalars, vectors, matrices!");
  
  // Check type of input.
  int i;
  for (i = 0; i < nrhs; i++) {
      if (mxGetClassID(prhs[i]) != mxDOUBLE_CLASS) {
          mexErrMsgTxt("Input must be of type double!");
      }
      // TODO: check for integer! 
  }
  
  int m = 1;
  int n = 1;

  if (nrhs == 1) {
      // create a random vector
      m = (int)armaGetDouble(prhs[0]);
      n = m;
  }
  if (nrhs == 2) {
      // create a random matrix
      m = (int)armaGetDouble(prhs[0]);
      n = (int)armaGetDouble(prhs[1]);
  }
  mat M(m,n,fill::randn);
  
  // Create the output argument plhs[0] to return matrix
  plhs[0] = armaCreateMxMatrix(m, n);
  
  // Return the matrix as plhs[0] in Matlab  
  armaSetPr(plhs[0], M);
  return;
}



