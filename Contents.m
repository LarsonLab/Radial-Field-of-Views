% RADIAL_FOVS - Tools to design FOVs for radial imaging
%
%    This package is used to design field-of-view's (FOVs) for radial
%    imaging techniques, commonly used in Computed Tomography (CT)
%    and as well as Magnetic Resonance Imaging (MRI).  Non-circular
%    and non-spherical shapes are supported by this package.
%    
%    Radial imaging acquires data on projections, or spokes, in
%    frequency space.  These functions will return the projection
%    acquisition angles that will support the input FOV shape and
%    size, and there are functions for both 2D and 3D imaging.
%    This package also designs acquisitions for the PROPELLER MRI
%    trajectory.
%
% Design functions:
%
%    calc_angles       - Computes 2D projection angles for the
%                        desired FOV (required by all other functions)
%    calc_cones        - 3D cones design
%    calc_3Dangles     - 3D projection angles with cones-based design
%    calc_3Dspiral     - 3D projection angles with spiral-based design
%    design_blades     - PROPELLER blade angles and line spacing design
%    design_blades_epi - PROPELLER-EPI trajectory design
%
% Additional C code (not required for above functions):
%
%    radial_fov.c        - Basic functions for desiging anisotropic FOV shapes
%                          Currently only includes elliptical FOVs
%    calc_3dpr_ellipse.c - Code that can be compiled as a MEX-function for
%                          designing 3dpr with elliptical FOVs, and uses the
%                          radial_fov.c functions
%
%
% The package also includes FOV shape functions, see
%    help shape_fcns
%
% Updated versions of this package are available at 
%    http://www-mrsrl.stanford.edu/~peder/radial_fovs/
%
% This package accompanies the journal article and conference abstract:
%    Larson PEZ, Gurney PT, Nishimura DG. "Anisotropic
%    Field-of-Views in Radial Imaging." IEEE - Transactions on
%    Medical Imaging. 27(1): 47-57 (2008).
%
%    Larson PEZ, Nishimura DG. "Anisotropic
%    Field-of-Views for PROPELLER MRI."
%    15th Annual Meeting of the International Society for Magnetic 
%    Resonance in Medicine. #1726 (2007).
%    Magnetic Resonance in Medicine, In Press, (2009).
%
% Written by Paul Gurney and Peder Larson, 6/28/2006
%   Last updated Jan. 29, 2009
% (c) 2007, 2008 Board of Trustees, Leland Stanford Junior University
