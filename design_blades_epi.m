function [alpha, L, dka, Nr, Nramp] = design_blades_epi(T, dt, norm_slew, rot_est, FOV, varargin);
% [alpha, L, dka, Nr, Nramp] = 
%        design_blades_epi(N, rot_est, FOV, F1, F2, ...)
%    OR
% [alpha, L, dka, Nr, Nramp] = 
%        design_blades_epi(T, dt, norm_slew, rot_est, FOV, F1, F2, ...)
%
% Designs a set of PROPELLER-EPI blades that are tailored to the desired
% field-of-view (FOV) size and shape.
% Some FOV functions are provided in the shape_fcns directory - use
% "help shape_fcns" for more information.
% The FOV size is specified in pixels, where trajectory has a
% resolution of 1 pixel.
%
% Inputs:
%   N - number of samples per blade.  
%    ** Designs the blades for an equal number of samples per
%       blade, not including gradient ramps.
%      OR
%   T - readout duration (s)
%   dt - sampling interval (s)
%   norm_slew - normalized slew-rate = slew-rate * resolution  (G/s)
%    ** Designs the blades for equal readout durations and includes
%       gradient ramps
%
%   rot_est - anticipated rotation estimate, in radians
%   FOV - function handle of the desired FOV (@fcn_name).
%         This function must be pi-periodic and convex
%   F1, F2, ... - Input parameters to FOV function
%
% Outputs:
%   alpha - blade angles, in radians
%   L - lines per blade
%   dka - line spacing in normalized k-space units
%   Nr - number of samples per line
%   Nramp - number of samples on gradient ramps
%   
% Examples:
%   % Add shape_fcns directory to the path
%   addpath([pwd '/shape_fcns'])
%
%   % Design an elliptical FOV with 1400 samples per blade
%   X = 80;  Y = 120;  N = 1400; rot_est = 10*pi/180;
%   [alpha, L, dka, Nr] = design_blades_epi(N, rot_est, @ellipse, X, Y);
%
%   % Rectangular FOV with 24ms readout duration and 1mm resolution
%   X = 80;  Y = 50;  rot_est = 0;
%   T = 18e-3; dt = 16e-6; % in sec
%   res = 0.1; norm_slew = 15000 * res; 
%   [alpha, L, dka, Nr, Nramp] = ...
%           design_blades_epi(T, dt, norm_slew, rot_est, @rect, X, Y);
%
% Peder Larson, 12/14/2007, last modified 1/23/2008
% (c) 2007, 2008, Board of Trustees, Leland Stanford Junior University

    
  % Can be modified, but usually not important
  alpha0 = 0;
  
  % parse inputs
  if isa(norm_slew, 'function_handle')
    if exist('FOV')
      F = {rot_est, FOV, varargin{:}};
    else
      F = {rot_est};
    end
    FOV = norm_slew;
    rot_est = dt;
    norm_slew = Inf;
    dt = 1;
  else
    F = varargin;
  end
  
  dbg = 0;
  
  % variable kmax also possible
  kmax = 1/2;

  sc_tol = 1e-4;
  Fdel_min = 1e-4;
  Fscale = 1;
  Fdel = .01;

  alphawid = pi;
  end_angle = alpha0 + alphawid;

  % initialize variables so they can be used by all inlined functions
  alpha = []; sc = [];
  dka = [];
  L = []; Nr = [];
  Nramp = [];
  
  alpha_best = []; sc_best = [];
  dka_best = [];  L_best = []; Nr_best = [];
  Nramp_best = [];

  
  do_design;
  
  function do_design;

    alpha(1) = alpha0;
    dka(1) = calc_dka(alpha(1), rot_est, Fscale, FOV, F);

    % must satisfy dkr constraint and baseline samples per line
    L(1) = getL(alpha(1), T, dt, norm_slew, rot_est, Fscale, FOV, F);
    
    Nr(1) = calc_Nr( L(1), T, dt, norm_slew );
    Nramp(1) = calc_ramp(Nr(1), dt, norm_slew);
    
    n = 1;
    
    while (alpha(n) < end_angle)
      da1 = atan( L(n) * dka(n) / (2*kmax) );
      
      [da2l, Fda2l] = fzero(@(da) F_next_blade(da, alpha(n)+da1, T, dt, ...
          norm_slew, rot_est, kmax, Fscale, FOV, F), 0);
      [da2h, Fda2h] = fzero(@(da) F_next_blade(da, alpha(n)+da1, T, dt, ...
          norm_slew, rot_est, kmax, Fscale, FOV, F), pi/4);
      Ll = getL(alpha(n)+da1+da2l, T, dt, norm_slew, rot_est, Fscale, FOV, F);
      Lh = getL(alpha(n)+da1+da2h, T, dt, norm_slew, rot_est, Fscale, FOV, F);
    
      if Ll < Lh
        da2 = da2l;
        L(n+1) = Ll;
      else
        da2 = da2h;
        L(n+1) = Lh;
      end
      
      alpha(n+1) = alpha(n) + da1 + da2;
      dka(n+1) = calc_dka(alpha(n+1), rot_est, Fscale, FOV, F);
      Nr(n+1) = calc_Nr(L(n+1), T, dt, norm_slew );
      Nramp(n+1) = calc_ramp(Nr(n+1), dt, norm_slew);
      
      n = n+1;
    end
  end % end of do_design inline function

  Idup = length(alpha);
  get_sc(Idup);

  function get_sc(I);
    sc = alphawid / (alpha(I) - alpha0);
  end
  
  Nblades = length(alpha)-1;

  set_best;
  Idup_best = Idup_best;

  function set_best
    alpha_best = alpha;
    dka_best = dka;
    L_best = L;
    Nr_best = Nr;
    Nramp_best = Nramp;
    sc_best = sc;
    Idup_best = Idup;
  end
    
  
  while (abs(1 - sc)/sc > sc_tol)
    Fscale = Fscale*(1+Fdel); % small increases in FOV size

    clear alpha dka L Nr Nramp
    
    do_design;

    if (length(alpha)-1 == Nblades)
      Idup = length(alpha);
    else
      Idup = length(alpha)-1;
      % Roll back Fscale
      Fscale = Fscale / (1+Fdel);
      Fdel = Fdel*.5;
      if Fdel < Fdel_min
        break;
      end
    end

    get_sc(Idup);
  
    if abs(1-sc)/sc < abs(1-sc_best)/sc_best
      set_best;
    end
    
    if dbg >= 1
      disp(['sc = ' num2str(sc) ', (1-sc)/sc = ' num2str(1/sc-1) ', Nblades = ' int2str(length(alpha))])
      disp(['L = ' int2str(L)])
      disp(['Nr = ' int2str(Nr)])
      disp(['Nramp = ' int2str(Nramp)])
      disp(['Tblades = ' num2str(L.*dt.*(Nr+Nramp))])
      
    end
  
  end

  scale_params(Idup_best);

  function scale_params(I);
    alpha = (alpha_best(1:I-1) - alpha0) * sc_best + alpha0;
    dka = dka_best(1:I-1) * sc;
    L = L_best(1:I-1);
    Nr = Nr_best(1:I-1);
    Nramp = Nramp_best(1:I-1);
  end
  
end

% for calculating line spacing, including anticipated rotation if specified
function dka = calc_dka(a, rot_est, Fscale, FOV, F);

  if rot_est == 0
    dka = 1 ./ (feval(FOV, a + pi/2, F{:})*Fscale);
  else
    Neval = 100;
    da = linspace(-rot_est, rot_est, Neval);
    dka = min( 1 ./ (feval(FOV, a + da + pi/2, F{:})*Fscale) );
  end

end

% for solving non-linear equations to determine the number of lines annd
% line spacing of the next blade
function Fval = F_next_blade(da, a0, T, dt, norm_slew, rot_est, kmax, Fscale, FOV, F)
  
  for n = 1:length(da)
    L(n) = getL(a0 + da(n), T, dt, norm_slew, rot_est, Fscale, FOV, F);
  end

  dka = calc_dka(a0 + da, rot_est, Fscale, FOV, F);

  Fval = L .* dka / (2*kmax)  - tan( da );
  
end

% calculate required number of lines
function [L, Nrmin, Nramp] = getL(a, T, dt, norm_slew, rot_est, Fscale, FOV, F)
  Nrmin = ceil( 1./ calc_dka(a-pi/2, rot_est, Fscale, FOV, F) );
  
  dka = calc_dka(a, rot_est, Fscale, FOV, F);
  
  da_corner = atan( 1 / (dka*Nrmin) );

  % Ensure the "corner" aliasing lobes are not inside the FOV
  while ( sqrt( 1/dka^2 + Nrmin^2) <= ...
		  1./calc_dka(a + da_corner -pi/2, rot_est, Fscale, FOV, F) ) | ...
        ( sqrt( 1/dka^2 + Nrmin^2) <= ...
          1./calc_dka(a - da_corner -pi/2, rot_est, Fscale, FOV, F) )
    Nrmin = Nrmin + 1;
    da_corner = atan( 1 / (dka*Nrmin) );
  end

  Nramp = calc_ramp(Nrmin, dt, norm_slew);
  
  L = floor( T/dt ./ (Nrmin + 2*Nramp) );
end

% calculate number of ramp samples
function Nramp = calc_ramp(Nr, dt, norm_slew)
  gr = 1 ./ (Nr * 4257 * dt); % *2*kmax
  Nramp = ceil(gr / (norm_slew*dt) );
end

% calculate number of readout samples
function Nr = calc_Nr(L, T, dt, norm_slew)
   Nr = floor( T./(2*L*dt) + ...
               sqrt( (T./L/dt).^2 - 4 * 2 /(norm_slew*4257*dt^2) )/2 );
end
