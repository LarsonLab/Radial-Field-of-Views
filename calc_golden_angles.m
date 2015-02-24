function [theta] = calc_golden_angles(FOV, Nprojections)
% [theta] = calc_golden_angles(FOV [, Nprojections])
%   
% Calculates a set of angles with "golden angle" ordering, supporting
% anisotropic field-of-views (FOVs).  Currently supports an elliptical
% FOV shape.
% 
% Inputs:
%   FOV = [FOVx, FOVy] - FOV dimensions, in pixels, in x,y directions on values here
%   Nprojections (optional) - total number of projections, where FOV ratio
%       defined is maintained
%   
% Outputs:
%   theta - resulting angles for desired FOV/KFCN
%           angles are for a full-spoke (0 < theta < pi) 
%
% Original authors: Wenwen Jiang, Peder Larson
% (c) 2014, The Regents of the University of California.

if length(FOV) == 1
    FOV = [FOV, FOV];
end

if nargin < 2
    % estimate total number of projections
    theta_sequential  = calc_angles(0, 2*pi, @ellipse, FOV(1), FOV(2));
    Nprojections = length(theta_sequential);
end

% Create parameterization curve from equally spaced, isotropic FOV angles
% to anisotropic FOV angles
theta_map_orig = linspace(0, 2*pi, 1e4);
% elliptical FOV formula
dtheta_map = sqrt( (cos(theta_map_orig).^2 / (FOV(2)/2)^2 + sin(theta_map_orig).^2 / (FOV(1)/2)^2));

theta_map_param = cumsum(dtheta_map) / sum(dtheta_map) * 2*pi;

d_golden_angle = pi*(3-sqrt(5));

theta(1) = 0;
theta_isotropic(1) = 0;
for n = 2:Nprojections
    theta_isotropic(n) = mod(theta_isotropic(n-1) + d_golden_angle, 2*pi);
    
    [temp Ith] = min(abs(theta_map_orig - theta_isotropic(n)));
    theta(n) = theta_map_param(Ith);
end

