/* calc_3dpr_ellipse.c - For compiling a MEX-function in Matlab
 * 
 * Utilizes radial_fov.c functions for designing 3D PR trajectories.
 * Compile in Matlab with 'mex calc_3dpr_ellipse.c'
 * This also serves as a good example of how to use the radial_fov.c functions.
 *
 * Peder Larson, 10/19/2006
 * (c) 2006, Board of Trustees, Leland Stanford Junior University
 */

#include "mex.h" 
#include "radial_fov.c"

#define MAX_3D_ANGLES 1000000

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


{
double *fovptr;
double fovx; double fovy; double fovz;
double resxy; double resz;
int *nangles;
double *theta; double *phi; double *kmax; double *dcf;

/* read inputs */
fovptr = mxGetPr(prhs[0]);
fovx = fovptr[0];
fovy = fovptr[1];
fovz = fovptr[2];
resxy = mxGetScalar(prhs[1]);
resz = mxGetScalar(prhs[2]);

/* unused variable */
nangles = mxMalloc(sizeof(int));

plhs[0] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
theta = mxGetPr(plhs[0]);
plhs[1] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
phi = mxGetPr(plhs[1]);
plhs[2] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
kmax = mxGetPr(plhs[2]);
plhs[3] = mxCreateDoubleMatrix(1,MAX_3D_ANGLES, mxREAL);
dcf = mxGetPr(plhs[3]);
	
calc_3d_angles(nangles, theta, phi, kmax, dcf, fovx, fovy, fovz, resxy, resz);

mxSetN(plhs[0], *nangles);
mxSetN(plhs[1], *nangles);
mxSetN(plhs[2], *nangles);
mxSetN(plhs[3], *nangles);

}





