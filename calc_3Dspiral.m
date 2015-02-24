function [theta, phi, kmax, dcf] = calc_3Dspiral(proj_type, FOVtheta, varargin)
% [theta, phi, kmax, dcf] = calc_3Dspiral([proj_type,] FOVtheta, Ft1, Ft2, ... 
%                                         FOVphi, Fp1, Fp2, ...
%                                         [KFCNtheta, Kt1, Kt2])
%
% Calculates a set of angles and k-space extents for desired 3D PR
% field-of-view (FOV) using a spiral on sphere method.
% Spherical coordinates are returned, where theta is the deflection
% from the z-axis (polar angle), and phi is the deflection from the
% x-axis (azimuthal angle).
% The resulting FOV = min(FOVtheta, FOVphi), where FOVtheta is
% the polar FOV (circularly symmetric about the z-axis) and FOVphi
% is the azimuthal FOV which is invariant in z.
%
% Requires calc_angles.m. Some basic FOV/KFCN shapes are provided
% in shape_fcns directory - use "help shape_fcns" for more
% information on these functions.
%
% Inputs:
%   proj_type (optional) - Type of projections (or spokes, echoes) desired:
%                          0 - full projections/spokes
%                          1 - half projections (default)
%   FOVtheta - desired FOV shape in theta, input as a function
%         handle (@fcn_name).
%   Ft1, Ft2, ... - Input parameters to FOVtheta function
%   FOVphi - desired FOV in phi, input as a function handle
%   Fp1, Fp2, ... - Input parameters to FOVphi function
%   KFCNtheta (optional) - function of desired kmax in theta
%         (defaults to constant)
%   Kt1, Kt2, ... - Inputs to KFCNtheta function
%         Variation of kmax in phi is not supported
%
% Outputs:
%   theta, phi - resulting spherical angle coordinates of projections
%                theta is deflection from the kz axis, 
%                and phi is deflection from the kx axis
%   kmax - trajectory extents for corresponding theta/phi's
%   dcf - angular density compensation weighting function
%         assumes equal number of points acquired on each spoke
%
% Examples:
%   % Cylindrical FOV with a diameter D and height Z
%   D = 40; Z = 80;
%   [theta, phi, kmax, dcf] = calc_3Dspiral(@rect, Z, D, @const, D);
%
%   % Cylinder with and ellipsoid base with diameters of X and Y
%   X = 20; Y = 40;
%   [theta, phi, kmax, dcf] = calc_3Dspiral(@rect, Z, max(X,Y), @ellipse, X, Y);
%
%   % Ellipsoid FOV using full-projection acquisitions
%   X = 60; Y = 60; Z = 30;
%   [theta, phi, kmax, dcf] = calc_3Dspiral(0, @ellipse, Z, max(X,Y), @ellipse, X, Y);
%
%   % Ellipsoid FOV with elliptical kmax shape in theta
%   D = 30; Z = 60;
%   [theta, phi, kmax, dcf] = calc_3Dspiral(@ellipse, Z, D, @const, D, ...
%                                           @ellipse, D/Z, 1);
%
%   % Cuboid
%   X = 30; Y = 30; Z = 60;
%   [theta, phi, kmax, dcf] = calc_3Dspiral(@rect, Z, sqrt(X^2+Y^2), @rect, X, Y);
%
% Paul Gurney and Peder Larson, 6/19/2006
% (c) 2006, Board of Trustees, Leland Stanford Junior University
  
% Parse input to functions:
% makes proj_type optional
if isa(proj_type, 'function_handle')
  varargin(2:end+1) = varargin;
  varargin{1} = FOVtheta;
  FOVtheta = proj_type;
  proj_type = 1;
end

% get FOV/KFCN inputs
Ifcn = [];
for k = 1:length(varargin)
  if isa(varargin{k}, 'function_handle')
    Ifcn = [Ifcn k];
  end
end

if (length(Ifcn) == 0)
  error('FOVphi not specified')
end

% FOVphi
Ft = varargin([1:Ifcn(1)-1]);
FOVphi = varargin{Ifcn(1)};

if (length(Ifcn) == 1)
  Fp = varargin([Ifcn(1)+1:end]);
  KFCNtheta = @const;
  Kt = {1};
else
  Fp = varargin([Ifcn(1)+1:Ifcn(2)-1]);
  KFCNtheta = varargin{Ifcn(2)};
  Kt = varargin([Ifcn(2)+1:end]);
end

% Wong & Roos style 3D PR design

if proj_type % half projections
  thwid = pi;
else % full projections
  thwid = pi/2;
end
  
[thcones, kmaxcones, dcones] = calc_angles(0, thwid, FOVtheta, Ft{:}, KFCNtheta, Kt{:});
thcones(end+1) = thwid;
kmaxcones(end+1) = feval(KFCNtheta, thwid, Kt{:})/2;
%dcones(end+1) = kmaxcones(end)/feval(FOVtheta, pi, Ft{:});

t(1) = 1;

kmaxest = feval(KFCNtheta, pi/2, Kt{:});
phiest = calc_angles(0, 2*pi, FOVphi, Fp{:}, @const, kmaxest);
Nphiest = length(phiest);
  
ncones = length(thcones);
for k = 2:ncones
  t(k) = t(k-1) + Nphiest * sin((thcones(k) + thcones(k-1))/2) * ...
                                (kmaxcones(k) + kmaxcones(k-1)) / kmaxest;  
end

if (proj_type == 0)
  % add extra quarter turn of spiral
  thcones(end+1) = pi/2 + 1/(kmaxest/2 * feval(FOVtheta, pi, Ft{:})) / 4;

  t(ncones+1) = t(ncones) + Nphiest/4;

  kmaxcones(end+1) = feval(KFCNtheta, thcones(end), Kt{:})/2;
end
  
n = 1:floor(t(end));
theta = interp1(t, thcones, n, 'linear');
kmax = interp1(t, kmaxcones, n, 'linear');

phi(1) = 0;
for k = 2:length(n)
  rk = kmax(k) * sin(theta(k));
    
  dphi_approx = 1 / (rk * feval(FOVphi, phi(k-1) + pi/2, Fp{:}));

  dphi = 1 / (rk * feval(FOVphi, phi(k-1) + dphi_approx/2  + pi/2, Fp{:}));
  phi(k) = phi(k-1) + dphi;
end

dcf = (kmax./feval(FOVtheta, theta + pi/2, Ft{:})) ./ feval(FOVphi, phi + pi/2, Fp{:});

if (proj_type == 0) % adjust dcf for full projections near k_xy
  I1 = find((phi > phi(end) - 2*pi) & (phi <= phi(end) - pi));
  dcf(I1) = dcf(I1) .* (1 - .5 * (phi(I1)- (phi(end) - 2*pi)) / pi );
  I2 = find(phi > phi(end) - pi);
  dcf(I2) = dcf(I2) .* (1 - .5 * (phi(I2)- (phi(end) - pi)) / pi);
end

phi = mod(phi, 2*pi);