% SHAPE_FCNS - Shaping functions used in RADIAL_FOVS package
%
%    These functions return the polar coorindates for a shape,
%    where the distance returned is the diameter, measured
%    through the center at a given angle.
%
%    They are intended for use with the RADIAL_FOVS package as
%    shapes for field-of-views and also acquisition shapes.
%    See "help radial_fovs" for more information.
%
%
% Shape functions: (input arguments)
%   const: (D) - circle of diameter D
%   ellipse: (X, Y) - ellipse with axes diameters of X and Y
%   rect: (X, Y) - rectangle with side lengths of X and Y
%   oval: (X, Y, R) - oval (rounded rectangle) with side lengths 
%        of X and Y, and corner radius of R
%   diamond: (X, Y) - diamond shape with diameters of X and Y
%   star: (D, S) - star/squished circle of diameter D with a "squish" factor S
%        (S = 0 will yeild a circle)
%   polygon: (D, N) - N-sided polygon with a diameter D measured
%        from the middle of the sides
%
% The resulting FOV shape can be demonstrated using:
%    theta = [0:.01:1]*2*pi; polar(theta, FCN(theta, input arguments)/2);
%
% When used with calc_XXX functions for determining projection angles,
% function handles (@fcn_name) followed by the input arguments shown above
% should be passed in.  See the design functions for examples.
% 
% Written by Paul Gurney and Peder Larson, 6/28/2006
% (c) 2006, Board of Trustees, Leland Stanford Junior University