function d = star(angle, D, S)
% d = star(angle, D, S)
%   Star/squished circle shape with a diameter of D and a
%   "squish factor" S (0 <= S < 1).  S = 0 results in a circle.
%
% Peder Larson and Paul Gurney, 6/28/2006
if ((S < 0) | (S > 1))
  error('S must be between 0 and 1')
end

d = D*((cos(4*angle)+1)/2*S+(1-S));