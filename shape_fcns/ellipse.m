function d = ellipse(angle, X, Y)
% d = ellipse(angle, X, Y)
%   Ellipse with diameters of X and Y.
%   Returns diameter for input polar coordinates.
%
% Peder Larson and Paul Gurney, 6/28/2006
d = (cos(angle).^2/X^2 + sin(angle).^2/Y^2).^(-1/2);

