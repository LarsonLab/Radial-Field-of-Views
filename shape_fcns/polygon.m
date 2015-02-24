function out = polygon(angle, D, N)
% out = polygon(angle, D, N)
%   Returns the diameter of an N-sided polygon with a diameter,
%   measured between the middle of the sides) of D.
%
%   N should be even if being used with the radial_fov package
%   design functions.
%
% Peder Larson, 6/28/2006
angle = abs(rem(angle, 2*pi/N)) - pi/N; % puts angles to be between -pi/N to pi/N
out = D*sec(angle);