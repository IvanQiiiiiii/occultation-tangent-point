% Compute COSMIC-style tangent ("impact") point from two ECEF positions.
% Units: km throughout. Requires Mapping or Aerospace Toolbox for ecef2geodetic (if method='geodetic').
% Copyright Owner: Qi Yifan (Ivan Qi) e-mail:qyfivan@cug.edu.cn
% version 1.0: 2025-10-24

function [p, h, glat, glon, meta] = impact_TangentPoint(x1, y1, z1, x2, y2, z2, opts)
% Why: provide two latitude/height definitions. 'radial' matches legacy; 'geodetic' gives standard geodetic height.
%
% Inputs:
%   x1,y1,z1  : ECEF of LEO (km)
%   x2,y2,z2  : ECEF of GPS (km)
%   opts.method : 'radial'  (default, matches original math)
%                 'geodetic' (recommended if you want standard geodetic height)
%   opts.ellipsoid : referenceEllipsoid/wgs84Ellipsoid in kilometers (default WGS-84 km)
%
% Outputs:
%   p     : impact parameter = min distance to origin along LOS (km)
%   h     : height (km) per chosen method
%   glat  : geodetic latitude (deg)
%   glon  : longitude (deg)
%   meta  : struct (fields: method, betweenSegment, s1, s2, n, ellipsoid)

arguments
    x1 (:,:) double; y1 (:,:) double; z1 (:,:) double
    x2 (:,:) double; y2 (:,:) double; z2 (:,:) double
    opts.method (1,1) string {mustBeMember(opts.method,["radial","geodetic"])} = "radial"
    opts.ellipsoid = wgs84Ellipsoid("kilometer")
end

% basic shapes
x1 = double(x1); y1 = double(y1); z1 = double(z1);
x2 = double(x2); y2 = double(y2); z2 = double(z2);

% line-of-sight unit vector
xn = x2 - x1; yn = y2 - y1; zn = z2 - z1;
rabs = hypot(hypot(xn,yn),zn);
if any(rabs(:) == 0), error('ro:impact_tangent:degenerateLOS','r1 and r2 coincide.'); end
xn = xn ./ rabs; yn = yn ./ rabs; zn = zn ./ rabs;

% projections
sp1 = x1.*xn + y1.*yn + z1.*zn;
sp2 = x2.*xn + y2.*yn + z2.*zn;

% closest point vector to origin (impact vector)
px = x1 - xn.*sp1;
py = y1 - yn.*sp1;
pz = z1 - zn.*sp1;

% impact parameter
p = hypot(hypot(px,py),pz);

% longitude (robust at poles)
glon = atan2(py, px) * 180/pi;
rho  = hypot(px,py);
r    = hypot(rho, pz);

% ellipsoid axes (km)
sph = opts.ellipsoid;
a0  = sph.SemimajorAxis;  % remax
b0  = sph.SemiminorAxis;  % remin

switch opts.method
    case "radial"
        % original mapping: geocentric -> geodetic via tan(phi) = (a/b)^2 tan(theta)
        % handle rho==0 to avoid 0/0
        glat = zeros(size(p));
        % where rho>0
        idx = rho > 0;
        d2  = (a0/b0)^2;
        glat(idx) = atan( d2 * (pz(idx)./rho(idx)) ) * 180/pi;
        % poles
        glat(~idx & (pz>=0)) = +90;
        glat(~idx & (pz<0))  = -90;

        % radial height above ellipsoid along the same geocentric direction
        h = r .* ( 1 - (a0*b0) ./ sqrt( (a0^2)*(pz.^2) + (b0^2)*(rho.^2) ) );

    case "geodetic"
        % standard geodetic: use MATLAB built-in for accuracy/stability (normal height)
        % note: ecef2geodetic returns deg + height same length unit as ellipsoid
        try
            [glat, glon_bi, h] = ecef2geodetic(px, py, pz, sph, 'degrees');
            % prefer lon from ecef2geodetic for consistency (wrap may differ by 360)
            glon = glon_bi;
        catch ME
            if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
                error(['ecef2geodetic not found. Install Mapping/Aerospace Toolbox, ', ...
                    'or use opts.method="radial".']);
            else
                rethrow(ME);
            end
        end
end

% meta info
meta = struct();
meta.method = char(opts.method);
meta.ellipsoid = sph;
meta.s1 = sp1; meta.s2 = sp2;
meta.betweenSegment = (sp1.*sp2) < 0;  % true => closest point lies between the two satellites
meta.n = cat(3, xn, yn, zn);           % LOS unit vector

end