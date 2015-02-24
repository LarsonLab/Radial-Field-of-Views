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

#define MAX_2D_ANGLES 4096
#ifndef PI
#define PI 3.1415926535897932384626433832795028841971
#endif

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])


{
double *fovptr;
double fovx; double fovy;
double *resptr;
double resx; double resy;
int *nangles;
double *theta; double *kmax; double *dcf;

/* read inputs */
fovptr = mxGetPr(prhs[0]);
fovx = fovptr[0];
fovy = fovptr[1];
resptr = mxGetPr(prhs[1]);
resx = resptr[0];
resy = resptr[1];

/* unused variable */
nangles = mxMalloc(sizeof(int));

plhs[0] = mxCreateDoubleMatrix(1,MAX_2D_ANGLES, mxREAL);
theta = mxGetPr(plhs[0]);
plhs[1] = mxCreateDoubleMatrix(1,MAX_2D_ANGLES, mxREAL);
kmax = mxGetPr(plhs[1]);
plhs[2] = mxCreateDoubleMatrix(1,MAX_2D_ANGLES, mxREAL);
dcf = mxGetPr(plhs[2]);
	
calc_2d_angles(nangles, theta, kmax, dcf, fovx, fovy, resx, resy, 2.0*PI);

mxSetN(plhs[0], *nangles);
mxSetN(plhs[1], *nangles);
mxSetN(plhs[2], *nangles);

}





