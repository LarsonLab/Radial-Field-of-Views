function [theta, phi, kmax, dcf] = calc_3Dangles(proj_type, FOVtheta, varargin)
% [theta, phi, kmax, dcf] = calc_3Dangles(proj_type, FOVtheta, Ft1, Ft2, ... 
%                                         FOVphi, Fp1, Fp2, ...
%                                         [KFCNtheta, Kt1, Kt2])
%
% Calculates a set of angles and k-space extents for desired 3D PR
% field-of-view (FOV) using a cones-based method. 
% Spherical coordinates are used, where theta is the deflection
% from the z-axis (polar angle), and phi is the deflection from the
% x-axis (aziumathal angle).
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
%                          0 - full projections/spokes (default)
%                          1 - half projections
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
%   [theta, phi, kmax, dcf] = calc_3Dangles(@rect, Z, D, @const, D);
%
%   % Cylinder with and ellipsoid base with diameters of X and Y
%   X = 20; Y = 40;
%   [theta, phi, kmax, dcf] = calc_3Dangles(@rect, Z, max(X,Y), @ellipse, X, Y);
%
%   % Ellipsoid FOV using half projection acquisitions
%   X = 60; Y = 60; Z = 30;
%   [theta, phi, kmax, dcf] = calc_3Dangles(1, @ellipse, Z, max(X,Y), @ellipse, X, Y);
%
%   % Ellipsoid FOV with elliptical kmax shape in theta
%   D = 30; Z = 60;
%   [theta, phi, kmax, dcf] = calc_3Dangles(@ellipse, Z, D, @const, Z, D, ...
%                                           @ellipse, D/Z, 1);
%
%   % Cuboid
%   X = 30; Y = 30; Z = 60;
%   [theta, phi, kmax, dcf] = calc_3Dangles(@rect, Z, sqrt(X^2+Y^2), @rect, X, Y);
%
% Paul Gurney and Peder Larson, 6/1/2006
% (c) 2006, Board of Trustees, Leland Stanford Junior University

% Parse input to functions:
% makes proj_type optional
if isa(proj_type, 'function_handle')
  varargin(2:end+1) = varargin;
  varargin{1} = FOVtheta;
  FOVtheta = proj_type;
  proj_type = 0;
end

% Parse input to functions:
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

if proj_type % half projections
  thwid = pi;
else
  thwid = pi/2;
end

% thcones is angle deflected from z-axis
[thcones, kmaxcones, dcones] = calc_angles(0, thwid, FOVtheta, Ft{:}, KFCNtheta, Kt{:});

% first cone is just a single spoke
theta = [thcones(1)];  phi = [0]; kmax = [kmaxcones(1)];
dcf = [dcones(1)/(2*feval(FOVphi, phi(1), Fp{:}))];
  
phi0_base = 1/(1/2 * feval(FOVphi, pi/2, Fp{:}));
for k = 2:length(thcones)
  rk = 2*kmaxcones(k)*sin(thcones(k));
  phi0 = phi0_base/rk * rand;
  
  [phiring, kmaxring, dring] = calc_angles(phi0, 2*pi, FOVphi, Fp{:}, @const, rk);
    
  theta = [theta, thcones(k)*ones(1,length(phiring))];
  phi = [phi, phiring];
  kmax = [kmax, kmaxcones(k)*ones(1,length(phiring))];
  dcf = [dcf, dcones(k) .* (dring/rk)];
end

if (proj_type == 0)
  % For full projection design
  % last cone is at theta = pi/2, and only goes halfway around in phi
  rk = feval(KFCNtheta, pi/2, Kt{:});
  phi0 = phi0_base/rk * rand;
  
  [phiring, kmaxring, dring] = calc_angles(phi0, pi, FOVphi, Fp{:}, @const, rk);
  
  theta = [theta, pi/2*ones(1,length(phiring))];
  phi = [phi, phiring];
  kmax = [kmax, kmaxring];
  dcf = [dcf, (rk/2)/feval(FOVtheta, pi, Ft{:}) .* (dring/rk)];
end