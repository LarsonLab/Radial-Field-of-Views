function d = diamond(angle, X, Y)
% d = diamond(angle, X, Y)
%   Diamond with vertical and horizontal diameters of X and Y, 
%   returned for polar coordinates.
%
% Paul Gurney and Peder Larson, 6/28/2006
d = sqrt((Y./((Y/X)+abs(tan(angle)))).^2 + ...
           (Y*abs(tan(angle))./((Y/X)+abs(tan(angle)))).^2);

