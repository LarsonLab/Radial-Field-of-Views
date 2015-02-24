function d = rect(angle, X, Y)
% d = rect(angle, X, Y)
%   Rectangle with side lengths of X and Y, returned for polar coordinates.
%
% Peder Larson, 6/28/2006
corangle = atan2(Y, X);
angle = abs(rem(angle, pi)); % forces angle to be between 0 and pi 
for k = 1:length(angle)
  if ((angle(k) < corangle) | (angle(k) > pi - corangle)) 
    d(k) = X/abs(cos(angle(k)));
  else
    d(k) = Y/abs(sin(angle(k)));
  end
end
