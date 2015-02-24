function [theta, kmax, dcf_theta] = calc_angles(theta0, theta_width, FOV, varargin)
% [theta, kmax, dcf_theta] = calc_angles([theta0, theta_width] FOV, F1, F2, ..., 
%                                        [KFCN, K1, K2, ...])
%   
% Calculates a set of angles and k-space extents for desired 
% radial imaging field-of-view (FOV)
% Some FOV/KFCN functions are provided in the shape_fcns directory.
% Use "help shape_fcns" for more information on these functions.
% 
% Inputs:
%   theta0 (optional) - initial starting angle, in radians.  Default is 0.
%   theta_width (optional, requires theta0 when used) - size of range of 
%         angles to design.  Used in 3D PR design.  Default is pi. 
%   FOV - function handle of the desired FOV (@fcn_name).
%         This function must be pi-periodic and convex
%   F1, F2, ... - Input parameters to FOV function
%   KFCN (optional) - function of desired kmax (defaults to constant)
%   K1, K2, ... - Inputs to KFCN function, see comments below on values here
% 
% Outputs:
%   theta - resulting angles for desired FOV/KFCN
%           angles are for a full-spoke (0 < theta < pi) 
%   kmax - trajectory extents for corresponding theta's
%   dcf_theta - angular density compensation weighting function
%
% For FOV, the inputs to the external functions correspond to the pixel sizes.
% For KFCN, they correspond to the inverse of the resolution pixel size.
%  (KFCN = @const, K1 = 0.5 will give a resolution size of 2 pixels, 
%   while K1 = 1 will give a resolution of 1 pixel)
%
% Examples:
%   % Add shape_fcns directory to the path
%   addpath([pwd '/shape_fcns'])
%
%   % Design a circular FOV
%   X = 80;  
%   [theta, kmax, dcf_theta] = calc_angles(@const, X);
%   % theta will consist of X * pi / 2 equally spaced angles
%
%   % Design a 2D elliptical FOV
%   X = 80; Y = 120;
%   [theta, kmax, dcf_theta] = calc_angles(@ellipse, X, Y);
%
%   % Rectangular FOV with rectangular kmax pattern
%   [theta, kmax, dcf_theta] = calc_angles(@rect, X, Y, @rect, 1, X/Y);
%   % NOTE: resolution will be 1 pixel in x, and Y/X in y
%
%   % Oval FOV with initial angle of pi/2
%   theta0 = pi/2; R = 40;
%   [theta, kmax, dcf_theta] = calc_angles(theta0, @oval, X, Y, R);
%
%   % Diamond FOV using half projection acquisitions
%   [theta, kmax, dcf_theta] = calc_angles(0, 2*pi, @diamond, X, Y);
%
% Paul Gurney and Peder Larson, 12/12/2005, updated 5/30/2006
% (c) 2006, Board of Trustees, Leland Stanford Junior University

% Make sure shape_fcns is in the path
currpath = path;
if (isempty(strfind(currpath, 'shape_fcns')))
  error('shape_fcns directory must be in the path')
end

% Make theta0, theta_width optional
if isa(theta0, 'function_handle')
  varargin(3:end+2) = varargin;
  varargin{1} = theta_width;
  if nargin > 2
    varargin{2} = FOV;
  end
  FOV = theta0;
  theta0 = 0;
  theta_width = pi;
elseif isa(theta_width, 'function_handle')
  varargin(2:end+1) = varargin;
  varargin{1} = FOV;
  FOV = theta_width;
  theta_width = pi;
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

n = 1;
theta(n) = theta0;
theta_cutoff = theta0 + theta_width;

while (theta(n) < theta_cutoff)
  dtheta_approx = 1 / ( feval(KFCN, theta(n), K{:})/2 * feval(FOV, theta(n) + pi/2, F{:}));
  kmax_mid = feval(KFCN, theta(n) + dtheta_approx/2, K{:})/2;

  dtheta = 1 / ( kmax_mid * feval(FOV, theta(n) + dtheta_approx/2 + pi/2, F{:}));
  theta(n+1) = theta(n) + dtheta;
  n = n+1;
end

% adjust theta for symmetry
% choose adjustment based on which spoke is closest to pi
% NOTE: could also choose adjustment based on whether an even or odd number 
% of full-spokes is desired
if (theta(end) - (theta_cutoff) > (theta_cutoff) - theta(end-1))
  theta = (theta(1:end-2) - theta0)*theta_width/(theta(end-1) - theta0) + theta0;
else
  theta = (theta(1:end-1) - theta0)*theta_width/(theta(end) - theta0) + theta0;
end

% Calculate kmax function
kmax = feval(KFCN, theta, K{:})/2;

% DCF = projection length * projection spacing
dcf_theta = kmax./feval(FOV, theta + pi/2, F{:});