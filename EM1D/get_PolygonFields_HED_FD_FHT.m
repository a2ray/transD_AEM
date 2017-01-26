function Bz = get_PolygonFields_HED_FD_FHT(freqs,xyPolyTx,zTx,xyRx,zRx,sig,mu,z,filterName,GQorder)

% Calculates the frequency domain magnetic field of a polygon
% shaped loop of wire in the horizontal plane. 
% It does so by carrying out a line integral of the
% vertical magnetic field produced by infinitesimal horizontal electric dipoles
% alond the wire segments. The integration is performed using Gauss
% quadarature. Use get_HED_FD_FHT to compute Bz field from a point HED 
%

bPlot = false; % set to true to plot polygon and quadrature points
 
% Get quadrature weights:
[xq,wq] = GLegIntP(GQorder);
            
%number of polygon vertices 
nvertices = size(xyPolyTx,1);

%number of frequencies
nf = length(freqs);

%initialize the value of the integral to zero
Bz = zeros(1,nf);
 
if bPlot
   figure; 
end
 
%loop over the number of segments:
for k=1:nvertices
    
    % get the segments:  
    if (k == nvertices)
        v = [ xyPolyTx(k,:) ; xyPolyTx(1,:)  ];
    else
        v = [ xyPolyTx(k,:) ; xyPolyTx(k+1,:) ];
    end
   
    % Get length of wire segment and its midoint:
    
    dv   = diff(v);
    dx   = dv(1);
    dy   = dv(2);
    lenv = norm(dv);
    
    midpoint = sum(v,1)/2;
    
    % change quadarature interval from -1 to +1 to a to b. At the same
    % time, break out x and y location of quadrature points along the wire
    % segment:
    xquad = xq*dx/2 + midpoint(1);
    yquad = xq*dy/2 + midpoint(2);
    wquad = wq*lenv/2;  
    
    if bPlot
        plot(v([1:end 1],1),v([1:end 1],2),'k-'); hold on
        plot(xquad,yquad,'*')
        plot(midpoint(1),midpoint(2),'ro')
        axis equal
        
        plot(xyRx(1),xyRx(2),'k*')
    end
    % Now loop over quadrature points and compute Bz:
    
    BzQ = 0*Bz;
    
    for i = 1:length(xquad)
        
        % Get angle theta between wire segment and receiver azimuth. Since Bz
        % has cos(theta) dependence, the direction doesn't matter.
        
        dxRx = xyRx(1) - xquad(i);
        dyRx = xyRx(2) - yquad(i);
        
        dxdp = v(2,1) - xquad(i);
        dydp = v(2,2) - yquad(i);
        
        drxq = [dxRx dyRx];
        dv2q = [dxdp dydp];
      
        acrossb = (dydp*dxRx - dxdp*dyRx);
        adotb   = dv2q*drxq';
        theta   = atan2d( acrossb,adotb);    % tan(theta) =  a x b / a dot b   % this gets 4 quadarant theta
 
        rRx     = norm(drxq);
        
        BzP     = get_HED_FD_FHT(freqs,zTx,rRx,zRx,theta,sig,mu,z,filterName);
        
        BzQ     = BzQ + BzP*wquad(i); 
        
    end
    
    Bz = Bz + BzQ;
    
end

% normalize by true polygon area:
area = polyarea(xyPolyTx(:,1),xyPolyTx(:,2));

Bz = Bz/area;
 
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