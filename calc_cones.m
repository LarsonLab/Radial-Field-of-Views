function [theta, kmax, dcf] = calc_cones(FOV, varargin)
% [theta, kmax, dcf] = calc_cones(FOV, F1, F2, ..., [KFCN, K1, K2, ...])
%
% Calculates the angles and k-space extents of a set of cones for desired 
% field-of-view (FOV).
% Requires calc_angles.m. See "help shape_fcns" for some FOV/KFCN shapes.  
%
% Inputs:
%   FOV - desired FOV shape, which will be circularly symmetric about the z-axis.
%         Passed in as a function handle (@fcn_name).
%   F1, F2, ... - Input parameters to FOV function
%   KFCN (optional) - function of desired kmax (defaults to constant)
%   K1, K2, ... - Inputs to KFCN function
%
% Outputs:
%   theta - resulting cone deflections from the kz axis
%   kmax - trajectory extents for corresponding theta's
%   dcf - angular density compensation weighting function
%         assumes equal number of points acquired on each spoke
%
% Examples:
%
%   % Cylindrical FOV (Z is the FOV height)
%   Z = 160; D = 80;
%   [theta, kmax, dcf] = calc_cones(@rect, Z, D);
%
%   % Ellipsoid FOV with "star" kmax pattern
%   Z = 60; XY = 140; S = 0.25;
%   [theta, kmax, dcf] = calc_cones(@ellipse, Z, XY, @star, 1, S);
%
% Paul Gurney and Peder Larson, 6/1/2006
% (c) 2006, Board of Trustees, Leland Stanford Junior University

% Make sure shape_fcns is in the path
currpath = path;
if (isempty(strfind(currpath, 'shape_fcns')))
  error('shape_fcns directory must be in the path')
end

for k = 1:length(varargin)
  if isa(varargin{k}, 'function_handle')
    KFCN = varargin{k};
    K = varargin(k+1:end);
    break
  else
    F(k) = varargin(k);
  end
end

if (~exist('KFCN'))
  KFCN = @const;
  K = {1};
end

% Choose theta0 such that no cone will lie on the z-axis (theta = 0)
theta0 = 1/(2*feval(KFCN, 0, K{:})/2 * feval(FOV, pi/2, F{:}));

[theta, kmax, dcf] = calc_angles(theta0, FOV, F{:}, KFCN, K{:});