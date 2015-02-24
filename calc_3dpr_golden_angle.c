/* calc_3dpr_golden_angle.c - For compiling a MEX-function in Matlab
 * 
 * Calculates a 3D projection reconstruction golden angle trajectory ordering.
 * Currently only supports isotropic FOV
 * Compile in Matlab with 'mex calc_3dpr_golden_angle.c’
 *
 * Original authors: Wenwen Jiang, Peder Larson
 * (c) 2014, The Regents of the University of California.
 */

#include "mex.h" 
#include "math.h"

#define MAX_3D_ANGLES 1000000
#define PI 3.14159265358979323846
#define GA_PHI1 0.465571231876768
#define GA_PHI2 0.682327803828019

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


{
int num_projs, n;
double *theta; double *phi; double *kmax; double *dcf;

/* read inputs */
num_projs = (int) mxGetScalar(prhs[0]);

plhs[0] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
theta = mxGetPr(plhs[0]);
plhs[1] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
phi = mxGetPr(plhs[1]);
plhs[2] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
kmax = mxGetPr(plhs[2]);
plhs[3] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
dcf = mxGetPr(plhs[3]);
	
      double mphi1[MAX_3D_ANGLES], mphi2[MAX_3D_ANGLES];

      for (n=0; n<num_projs; n++) {
	mphi1[n] = (n > 0) ? fmodf(GA_PHI1 + mphi1[n-1],1.0) : GA_PHI1; 
	theta[n] = acos(mphi1[n]*2.0 - 1.0);

	mphi2[n] = (n > 0) ? fmodf(GA_PHI2 + mphi2[n-1],1.0) : GA_PHI2; 
	phi[n] = 2.0*PI*mphi2[n];

	kmax[n] = 1.0 / (2.0);
	dcf[n] = 1.0;
      }

mxSetN(plhs[0], num_projs);
mxSetN(plhs[1], num_projs);
mxSetN(plhs[2], num_projs);
mxSetN(plhs[3], num_projs);

}





