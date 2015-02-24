/* radial_fov.c - Radial FOV code, in C
 * 
 * Some support functions written in C for designing anisotropic FOV 
 * radial trajectories.
 * Currently includes support for both 2D and 3D PR.
 * See calc_3dpr_ellipse.c for an example of how to use these functions.
 *
 * Note: Only elliptical FOV shape is currently implemented.
 * The reason is that this is likely the most commonly desired shape,
 * and thus will be flexible for the most applications.
 *
 * Peder Larson, 10/19/2006
 * (c) 2006, Board of Trustees, Leland Stanford Junior University
 */

#include <math.h>

#define MAX_2D_ANGLES 4096
#ifndef PI
#define PI 3.1415926535897932384626433832795028841971
#endif
#define HPI  PI/2.0

void calc_3d_angles(int *nangles, double *theta, double *phi, double *kmax, double *dcf, double fov1, double fov2, double fov3, double res12, double res3);

void calc_interpolation_params(double *t, int *ncones, double *theta_cones, double *kmax_cones, double *dcf_cones, double fov1, double fov2, double fov3, double res12, double res3);

void calc_2d_angles(int *nangles, double *theta, double *kmax, double *dcf, double fov1, double fov2, double res1, double res2, double theta_width);

double interp_lin(double *x, double *Sx, int y);

double ellipse(double angle, double x, double y);



void calc_3d_angles(int *nangles, double *theta, double *phi, double *kmax, double *dcf, double fov1, double fov2, double fov3, double res12, double res3)
{
  int ncones, count;
  double t[MAX_2D_ANGLES];
  double theta_cones[MAX_2D_ANGLES], kmax_cones[MAX_2D_ANGLES], dcf_cones[MAX_2D_ANGLES];
  double del_phi_est, del_phi, fov12;


  fov12 = (fov1 > fov2) ? fov1 : fov2;

    /*
  mexPrintf("fovxy = %f, xfov = %f, yfov = %f\n", fov12, fov1, fov2);
  mexPrintf("resxy = %f, resz = %f\n", res12, res3);
  */

  calc_interpolation_params(t, &ncones, theta_cones, kmax_cones, dcf_cones, fov1, fov2, fov3, res12, res3);

  *nangles = (int) floor(t[ncones-1]) + 1;

  /*
  mexPrintf("Nangles = %d (Ncones = %d) \n", *nangles, ncones);
  mexPrintf("t[0, 1, ncones-2,ncones-1] = [%f,%f,%f,%f]\n", t[0], t[1], t[ncones-2], t[ncones-1]);
  */
  
  /* problem is before this point - dealing with the theta, phi*/

  for (count = 0; count < *nangles; count++) {
    
    theta[count] = interp_lin(t, theta_cones, count);

    kmax[count] = interp_lin(t, kmax_cones, count);

    if (count == 0)
      phi[count] = 0;
    else {
      del_phi_est = 1.0 / (kmax[count] * sin(theta[count]) 
			       * ellipse(phi[count-1] + HPI, fov1, fov2));
      del_phi = 1.0 / (kmax[count] * sin(theta[count]) 
			   * ellipse(phi[count-1] + HPI + del_phi_est/2.0, fov1, fov2));
      phi[count] = phi[count-1] + del_phi;

      /* move into 0 to 2*pi range */
      phi[count] -= floor(phi[count]/(2.0*PI)) * 2.0*PI; 

    }

    dcf[count] = kmax[count] / (ellipse(theta[count]+HPI, fov3, fov12) *
				ellipse(phi[count] + HPI, fov1, fov2));

  }

  return;
}


double interp_lin(double *x, double *Sx, int y)
{
  int count = 0;
  double frac;
  double Sy;

  while (1) {
    if (x[count] >= (double)y)
      break;
    count++;
  }

  if (count == 0)
    Sy = Sx[0];
  else {
    frac = (x[count] - (double)y) / (x[count]-x[count-1]);
    Sy = Sx[count] * (1.0 - frac) + Sx[count-1] * frac;
  }

  return (Sy);
}


void calc_interpolation_params(double *t, int *ncones, double *theta_cones, double *kmax_cones, double *dcf_cones, double fov1, double fov2, double fov3, double res12, double res3)
{
  int nphiest, count;
  double theta_temp[MAX_2D_ANGLES], kmax_temp[MAX_2D_ANGLES], dcf_temp[MAX_2D_ANGLES];
  double fov12;

  fov12 = (fov1 > fov2) ? fov1 : fov2;

  calc_2d_angles(ncones, theta_cones, kmax_cones, dcf_cones, 
		 fov3, fov12, res3, res12, PI);
  theta_cones[*ncones] = PI;
  kmax_cones[*ncones] = 1.0 / (2.0*res3);
  (*ncones)++;


  calc_2d_angles(&nphiest, theta_temp, kmax_temp, dcf_temp, 
		 fov1, fov2, res12, res12, 2.0*PI);

  t[0] = 0.0;

  for (count = 1; count < *ncones; count++)
    t[count] = t[count-1] + 
      (double)nphiest * sin( (theta_cones[count]+theta_cones[count-1])/2.0 ) *
      (kmax_cones[count]+kmax_cones[count-1]) * res12;

  return;
}

void calc_2d_angles(int *nangles, double *theta, double *kmax, double *dcf, double fov1, double fov2, double res1, double res2, double theta_width)
{
  int n = 0;
  double S;
  double del_theta_est, del_theta;

  theta[n] = 0.0;

  while (theta[n] < theta_width) {

    del_theta_est = 1.0 / ( ellipse(theta[n], 1.0 / (2.0*res1), 1.0 / (2.0*res2)) * 
			    ellipse(theta[n]+ HPI, fov1, fov2) );

    del_theta = 1.0 / ( ellipse(theta[n] + del_theta_est/2.0 , 1.0 / (2.0*res1), 1.0 / (2.0*res2)) * 
			ellipse(theta[n] + del_theta_est/2.0 + HPI, fov1, fov2) );
 
    theta[n+1] = theta[n] + del_theta;

    n++;

    /*
    mexPrintf("theta[%d] = %f, dtheta = %f\n", n, theta[n], del_theta);
    */
  }

  /* adjust angles for symmetry based on which spoke is closer to theta_width */

  if ( theta[n] - theta_width > theta_width - theta[n-1] ) {
    *nangles = n-1;
    S = theta_width / theta[n-1];
  } else {
    *nangles = n;
    S = theta_width / theta[n];
  }

  for (n = 0; n < *nangles; n++) {
    theta[n] *= S;
    kmax[n] = ellipse(theta[n], 1.0/(2.0*res1), 1.0/(2.0*res2));
    dcf[n] = kmax[n] / ellipse(theta[n]+ HPI, fov1, fov2);
  }

  return;
}



double ellipse(double angle, double x, double y)
{
  return ( 1.0 / sqrt( pow(cos(angle)/x,2.0) + pow(sin(angle)/y,2.0) ) );
}
