/*
 * For unit testing with ARMADILLO random matrices
 * it is nice to compare with MATLAB code with some existing
 * test functions.
 *
 * Note: We always use the same seed: 1!
 *
 */

#include <armadillo>
#include "armaMex.hpp"

using namespace arma;

void
mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  // Check the number of input arguments.
  if (nrhs != 1)
    mexErrMsgTxt("Currently supporting only one input!");
  
  // Check type of input.
  int i;
  for (i = 0; i < nrhs; i++) {
      if (mxGetClassID(prhs[i]) != mxDOUBLE_CLASS) {
          mexErrMsgTxt("Input must be of type double!");
      }
      // TODO: check for integer! 
  }
  
  int seed_value = (int)armaGetDouble(prhs[0]);
  arma_rng::set_seed(seed_value);
  
  return;
}