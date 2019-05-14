function Bz = get_PolygonFields_HED_FD_FHT(freqs,xyPolyTx,zTx,xyRx,zRx,sig,mu,z,filterName,GQorder,nSplinePointsPerDecade)

%
% Assumes single values for: zTx,xyRx,zRx
%
% x,y,z is right-handed with z positive down, so visualize as x=up, y=right, z=down
%
% Calculates the frequency domain magnetic field of a polygon
% shaped loop of wire in the horizontal plane. 
% It does so by carrying out a line integral of the
% vertical magnetic field produced by infinitesimal horizontal electric dipoles
% along the wire segments. The integration is performed using Gauss
% quadarature. Use get_HED_FD_FHT to compute Bz field from a point HED 
%

bPlot = false; % set to true to plot polygon and quadrature points

 
% Get quadrature weights:
[xq,wq] = GLegIntP(GQorder);
 
% Initialize the value of the integral to zero
Bz = zeros(1,length(freqs));
 
if bPlot
   figure; 
end

%
% KWK May 2019: noting that Bz at range R has sin(theta) symmetry, 
% efficient calculations are possible by first computing all quadrature
% points, computing Bz at theta = 90 (so sin(90) = 1), interpolating in R
% using splines, then applying the sin(theta) scaling. The old code did
% this but was hard wired for the SkyTEM Tx geometry, so we generalize this
% in the new code so that it works for any Tx geometry.
% 

% Loop over the number of Tx loop segments and collect quadrature points
% and weights in r and theta:

nvertices = size(xyPolyTx,1);
nq        = length(xq);
wquad     = zeros(nq*nvertices,1);
rquad     = zeros(nq*nvertices,1);
thetaquad = zeros(nq*nvertices,1);

for k=1:nvertices 
        
    % get the segments:  
    if (k == nvertices)
        v = [ xyPolyTx(k,:) ; xyPolyTx(1,:)  ];
    else
        v = [ xyPolyTx(k,:) ; xyPolyTx(k+1,:) ];
    end
   
    % Get length of wire segment and its midpoint:
    dv   = diff(v);
    dx   = dv(1);
    dy   = dv(2);
    lenv = norm(dv);
    midpoint = sum(v,1)/2;
    
    wquad((k-1)*nq + (1:nq))     = wq*lenv/2;
    
    
    % change quadarature interval from -1 to +1 to a to b. At the same
    % time, break out x and y location of quadrature points along the wire
    % segment:
    xquad = xq*dx/2 + midpoint(1);
    yquad = xq*dy/2 + midpoint(2);
    
    % Get angle theta between wire segment quad points and receiver azimuth: 
    dxRx = xyRx(1) - xquad;
    dyRx = xyRx(2) - yquad;
    
    rquad((k-1)*nq + (1:nq)) = sqrt(dxRx.^2+dyRx.^2); % distance from each quad point to Rx
    
    dxdp = v(2,1) - xquad;
    dydp = v(2,2) - yquad;    

    drxq = [dxRx dyRx];
    dv2q = [dxdp dydp];

    acrossb = (dydp.*dxRx - dxdp.*dyRx);
    adotb   = dot(dv2q',drxq')';
    
   	thetaquad((k-1)*nq + (1:nq)) = atan2d(acrossb,adotb);    % tan(theta) =  a x b / a dot b   % this gets 4 quadarant theta
 
        
    if bPlot
        plot(v([1:end 1],1),v([1:end 1],2),'k-'); hold on
        plot(xquad,yquad,'*')
        plot(midpoint(1),midpoint(2),'ro')
        axis equal
        plot(xyRx(1),xyRx(2),'k*')
    end
    
end

% Get spline points to use for ranges:
nDecadesRange = log10(max(rquad)/min(rquad));
nSplinePts = ceil(nDecadesRange)*nSplinePointsPerDecade;
% the above has at least nSplinePointsPerDecade points (for example when far field
% solution requested)
rSpline = logspace(log10(min(rquad)),log10(max(rquad)),nSplinePts);

% Now loop over spine points and compute Bz:
BzSpline = zeros(length(freqs),nSplinePts);
for j = 1:length(rSpline)
    BzSpline(:,j) = get_HED_FD_FHT(freqs,zTx,rSpline(j),zRx,90,sig,mu,z,filterName); % computed at theta=90 since Bz ~ sin(theta)
end
 
% Now interpolate to rquad points, scale by sin(thetaquad), mulitply by
% quadrature weight sand sum:

BzS = interp1(rSpline',BzSpline(:,:).',rquad ,'spline'); % note it is fastest to do this for all ranges at once so spline fit coefficents are computed only once
 
for iFreq = 1:length(freqs)
    Bz(iFreq) = sum(BzS(:,iFreq).*wquad.*sind(thetaquad));
end
 
 
% normalize by true polygon area:
area = polyarea(xyPolyTx(:,1),xyPolyTx(:,2));

Bz = Bz/area;

%keyboard
 
end


function [x,w]=GLegIntP(nGL)
% Gauss-Legendre integration points in the interval [-1,1].
%
% Description
%     [#x#,#w#]=GLegIntP(#nGL#)
%     calculates the abscissas (x) and weights (w) for the Gauss-Legendre
%     quadrature. The algorithm presented here is described in J. A.
%     Gubner, Gaussian quadrature and the eigenvalue problem, October 16,
%     2014.
%
% Input arguments
%     #nGL# (scalar) is the number of the integration points to be
%     calculated.
%
% Output arguments
%     #x# ([#nGL# x 1]) contains the coordinates of the integration points.
%     #w# ([#nGL# x 1]) contains the weights of the integration points.
%
% Parents (calling functions)
%     GLegQuad > GLegIntP
%
% Children (called functions)
%     GLegIntP >
%

%__________________________________________________________________________
% Contact author
%
%  (c) 2015 by George Papazafeiropoulos
%  First Lieutenant, Infrastructure Engineer, Hellenic Air Force
%  Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%
% Email: gpapazafeiropoulos@yahoo.gr
%
% Website: http://users.ntua.gr/gpapazaf/
%

beta=(1:nGL-1)./sqrt(4*(1:nGL-1).^2-1);
J=diag(beta,-1)+diag(beta,1);
[V,D]=eig(J);
[x,ix]=sort(diag(D));
w=2*V(1,ix)'.^2;

end