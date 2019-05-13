function [nLon, nLat] = UTM2LonLat( nE, nN, nZone, bSHemi, sEllipsoid, nFalseE, nFalseN )
% Convert E,N,Zone,Hemi(UTM) to Lon,Lat using WGS84. Geographic not geomagnetic.
% This function encapsulates the entire transformation without having to call
% the m_map library. It is easier to distribute. The guts come from m_map with
% attribution below.
%
% David Myer, Sept 2013
%-------------------------------------------------------------------------------
% See also LonLat2UTM

    % 6/26/2015 - Zones around Norway are special and have different false
    % eastings & northings. I occasionally encounter these done the correct way.
    if ~exist( 'nFalseE', 'var' )
        nFalseE = [];
    end
    if ~exist( 'nFalseN', 'var' )
        nFalseN = [];
    end

    % 12/2014 - support char input of hemisphere too
    if ischar(bSHemi)
        bSHemi = (bSHemi(1) == 'S');
    end
    
    % DGM 2/27/2015 - Support multiple ellipsoids
    if ~exist( 'sEllipsoid', 'var' ) || isempty( sEllipsoid )
        sEllipsoid = 'wgs84';
    end
    switch( lower( sEllipsoid ) )
    case 'wgs84' , nEllipsoid = [6378137.0, 1/298.257];
    case {'intl24', 'ed50'} % EuropeanDatum 1950 uses Intl 1924 ellipsoid params
                   nEllipsoid = [6378388.0, 1/297.000];
    case 'normal', nEllipsoid = [1.0, 0];
    case 'sphere', nEllipsoid = [6370997.0, 0];
    case 'grs80' , nEllipsoid = [6378137.0, 1/298.257];
    case 'grs67' , nEllipsoid = [6378160.0, 1/247.247];
    case 'wgs72' , nEllipsoid = [6378135.0, 1/298.260];
    case 'wgs66' , nEllipsoid = [6378145.0, 1/298.250];
    case 'wgs60' , nEllipsoid = [6378165.0, 1/298.300];
    case 'clrk66', nEllipsoid = [6378206.4, 1/294.980];
    case 'clrk80', nEllipsoid = [6378249.1, 1/293.466];
    case 'intl67', nEllipsoid = [6378157.5, 1/298.250];
    otherwise
        error( 'UTM2LonLat: Unknown / unsupported ellipsoid: "%s"', sEllipsoid );
    end
    
    [nLat,nLon] = mu_utm2ll( nE, nN, nZone, bSHemi, nEllipsoid, nFalseE, nFalseN );
    
    return;
end


%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% BELOW EXCERPTED FROM m_map's mp_utm.m function. Original header follows:
% 6/26/2015 DGM modified to support the weird zones around Norway
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
%-------------------------------------------------------------------------------
% MP_UTM   Universal Transverse Mercator projection
%           This function should not be used directly; instead it is
%           is accessed by various high-level functions named M_*.
%
% mp_utm.m, Peter Lemmond (peter@whoi.edu)
%
% created mp_utm.m 13Aug98 from mp_tmerc.m, v1.2d distribution, by:
%
% Rich Pawlowicz (rich@ocgy.ubc.ca) 2/Apr/1997
%
% This software is provided "as is" without warranty of any kind. But
% it's mine, so you can't sell it.
%
% Mathematical formulas for the projections and their inverses are taken from
%
%      Snyder, John P., Map Projections used by the US Geological Survey, 
%      Geol. Surv. Bull. 1532, 2nd Edition, USGPO, Washington D.C., 1983.
%
%
% 10/Dec/98 - PL added various ellipsoids.
function [lat,lon] = mu_utm2ll( x, y, zone, hemisphere, ellipsoid, FALSE_EAST, FALSE_NORTH )
%mu_utm2ll		Convert X/Y UTM coordinates to geodetic lat,lon 
%
%	[lat,lon] = mu_utm2ll (x,y, zone, hemisphere,ellipsoid)
%
%	input is X/Y vectors, zone number, hemisphere(N=0,S=1),
%		ellipsoid info [eq-rad, flat]
%	output is lat/lon vectors
%
%	see also	mu_ll2utm, utmzone


% some general constants

DEG2RADS    = 0.01745329252;
RADIUS      = ellipsoid(1);
FLAT        = ellipsoid(2);
K_NOT       = 0.9996;
if ~exist( 'FALSE_EAST', 'var' ) || isempty( FALSE_EAST )
    FALSE_EAST  = 500000;
end
if ~exist( 'FALSE_NORTH', 'var' ) || isempty( FALSE_NORTH )
    FALSE_NORTH = 10000000;
end

if ((zone < 1) || (zone > 60))
  error ('utm zones only valid from 1 to 60');
end

% compute some geodetic parameters

e2  = 2*FLAT - FLAT*FLAT;
e4  = e2 * e2;
e6  = e4 * e2;
eps = e2 / (1-e2);
em1 = sqrt(1-e2);
e1  = (1-em1)/(1+em1);
e12 = e1*e1;

lambda_not  = ((-180 + zone*6) - 3) * DEG2RADS;

% remove the false things

x = x - (RADIUS>1)*FALSE_EAST;
if (hemisphere)
  y = y - (RADIUS>1)*FALSE_NORTH;
end

% compute the footpoint latitude

M = y/K_NOT;
mu = M/(RADIUS * (1 - 0.25*e2 - 0.046875*e4 - 0.01953125*e6));
foot = mu + (1.5*e1 - 0.84375*e12*e1)*sin(2*mu) ...
    + (1.3125*e12 - 1.71875*e12*e12)*sin(4*mu) ...
    + (1.57291666667*e12*e1)*sin(6*mu) ...
    + (2.142578125*e12*e12)*sin(8*mu);

% some other terms

sinF = sin(foot);
cosF = cos(foot);
tanF = tan(foot);

N = RADIUS ./ sqrt(1-e2*(sinF.*sinF));
T = tanF.*tanF;
T2 = T.*T;
C = eps * cosF.*cosF;
C2 = C.*C;
denom = sqrt(1-e2*(sinF.*sinF));
R = RADIUS * em1*em1 ./ (denom.*denom.*denom);
D = x./(N*K_NOT);
D2 = D.*D;
D4 = D2.*D2;

% can now compute the lat and lon

lat = foot - (N.*tanF./R) .* (0.5*D2 - (5 + 3*T + 10*C - 4*C2 - 9*eps).*D4/24 ...
    + (61 + 90*T + 298*C + 45*T2 - 252*eps - 3*C2) .* D4 .* D2/720);

lon = lambda_not + (D - (1 + 2*T +C).*D2.*D/6 + ...
    (5 - 2*C + 28*T - 3*C2 + 8*eps + 24*T2).*D4.*D./120)./cosF;


% convert back to degrees;

lat=lat/DEG2RADS;
lon=lon/DEG2RADS;

return
end
