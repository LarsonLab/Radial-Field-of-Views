function [alpha, dk] = design_blades(alpha0, L, rot_est, FOV, varargin);
% [alpha, dk] = design_blades([alpha0,] L, rot_est, FOV, F1, F2, ...)
%
% Designs a set of PROPELLER blades for the desired field-of-view
% (FOV) shape and size.
% Some FOV functions are provided in the shape_fcns directory - use
% "help shape_fcns" for more information.
% The FOV size is specified in pixels, where trajectory has a
% resolution of 1 pixel.
%
% Inputs:
%   alpha0 (optional) - angle of initial blade, defaults to 0
%   L - number of lines acquired in each blade
%   rot_est - anticipated rotation estimate, in radians
%   FOV - function handle of the minimum required FOV shape
%         The resulting FOV will be this size or larger
%   F1, F2, ... - Input parameters to FOV function
%
% Outputs:
%   alpha - resulting blade angles
%   dk - line spacings on each blade 
%
% Examples:
%   % Design a 150x250 pixel elliptical FOV
%   L = 20; X = 150; Y = 250; rot_est = 10*pi/180;
%   [alpha, dk] = design_blades(L, rot_est, @ellipse, X, Y);
%   
%   % 180x200 pixel rectangular FOV
%   X = 180; Y = 200;
%   [alpha, dk] = design_blades(L, rot_est, @rect, X, Y);
%
% Peder Larson 8/11/2006, last updated 1/23/2008
% (c) 2007, 2008 Board of Trustees, Leland Stanford Junior University

% make alpha0 optional
if isa(rot_est, 'function_handle')
  F = {FOV, varargin{:}};
  FOV = rot_est;
  rot_est = L;
  L = alpha0;
  alpha0 = 0;
else
  F = varargin;
end

% variable kmax also possible
kmax = 1/2;

[alpha, dk, sc] = do_design(alpha0, rot_est, L/2 / kmax, 1, FOV, F);

sc_tol = 1e-4;
Fscale = 1;
Fdel = .01;
Fdel_min = 1e-4;
Nblades = length(alpha);

while (abs(1 - sc) > sc_tol)
  Fscalenew = Fscale*(1 + Fdel); % small increases in FOV size
  [alphanew, dknew, sc] = do_design(alpha0, rot_est, L/2 / kmax, Fscalenew, FOV, F);
  
  % protects to insure more blades are not used
  if (length(alphanew) > Nblades)
    if (Fdel < Fdel_min)
      break;
    else
      % shrink to converge on solution
      Fdel = Fdel/2;
    end
  else
    alpha = alphanew; dk = dknew;
    Fscale = Fscalenew;
  end
end


% internal function:
function [alpha, dk, sc] = do_design(alpha0, rot_est, ascale, Fscale, FOV, F);
% [alpha, dk, sc] = do_design(alpha0, rot_est, ascale, Fscale, FOV, F1, F2, ...)
%
% Does one iteration of a PROPELLER design.
%
% Inputs:
%   alpha0 - angle of initial blade
%   rot_est - estimate of anticipated rotation
%   ascale - geometrical scaling of blade angles (L/2 / kmax)
%   Fscale - scales FOV size
%   FOV - function handle of the FOV shape
%   F1, F2, ... - Input parameters to FOV function
%
% Outputs:
%   alpha - blade angles
%   dk - line spacings on each blade 
%   sc - scaling factor required for symmetry
%
% Peder Larson 11/10/2006
% (c) 2007 Board of Trustees, Leland Stanford Junior University

tol = .01;
alphawid = pi;

alpha(1) = alpha0;
dk(1) = calc_dk(alpha0, rot_est, Fscale, FOV, F);
end_angle = alpha0 + alphawid;

n = 1;

while (alpha(n) < end_angle)
  da1 = atan( ascale * dk(n) );

  % solving non-linear equation for line spacing of next blade
  [da2 Fda2] = fzero(@(da) F_next_blade(da, alpha(n) + da1, rot_est, ...
                                          ascale, Fscale, FOV, F), da1);
  
  alpha(n+1) = alpha(n) + da1 + da2;
  dk(n+1) = calc_dk(alpha(n+1), rot_est, Fscale, FOV, F);
  n = n+1;
end

% scale for symmetry
sc = alphawid / (alpha(end) - alpha0);
alpha = (alpha(1:end-1) - alpha0) * sc + alpha0;
dk = dk(1:end-1) * sc; % reduce oversampling of blades

%disp(['sc = ' num2str(sc) ', Nblades = ' int2str(length(alpha))])

% used to solve non-linear equation calculating line spacing
function Fval = F_next_blade(da, a0, rot_est, ascale, Fscale, FOV, F);
Fval = ascale * calc_dk(a0 + da, rot_est, Fscale, FOV, F) - tan(da);

% calculates line spacing, including anticipated rotation if specified
function dk = calc_dk(a, rot_est, Fscale, FOV, F);

if rot_est == 0
  dk = 1/ (feval(FOV, a + pi/2, F{:})*Fscale);
else
  Neval = 100;
  da = linspace(-rot_est, rot_est, Neval);
  dk = min( 1 ./ (feval(FOV, a + da + pi/2, F{:})*Fscale) );
end