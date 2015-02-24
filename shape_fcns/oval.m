function d = oval(angle, X, Y, R)
% d = oval(angle, X, Y, R)
%   Returns the radial diameter for input angle(s) of an oval with
%   side lengths of X and Y, and a corner radii of R
%
% Peder Larson, 6/28/2006

if (R > min(X/2,Y/2))
  error('Corner radius is too large (must be less than half the minimum side length)')
end

corangle1 = atan2(Y/2-R, X/2);
corangle2 = atan2(Y/2, X/2-R);

corangle = atan2(Y/2-R,X/2-R);
cordist = sqrt((X/2-R)^2 + (Y/2-R)^2);

angle = abs(rem(angle, pi)); % forces angle to be between -pi to pi 

for k = 1:length(angle)
  if ((angle(k) < corangle1) | (angle(k) > pi - corangle1)) 
    d(k) = X/abs(cos(angle(k)));
  elseif ((angle(k) > corangle2) & (angle(k) < pi - corangle2)) 
    d(k) = Y/abs(sin(angle(k)));
  else
    if (angle(k) < pi/2)
      tmpangle = abs(angle(k) - corangle);
    else
      tmpangle = abs(angle(k) - (pi - corangle));
    end
    d(k) = 2*(cordist*cos(tmpangle) + sqrt(R^2 - cordist^2*sin(tmpangle)^2));
  end
end
