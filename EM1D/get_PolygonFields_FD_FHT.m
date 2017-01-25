function Bz = get_PolygonFields_FD_FHT(freqs,xyPolyTx,zTx,xyRx,zRx,sig,mu,z,filterName,GQorder)

%This function calculates the frequency domain magnetic field of a polygon
%shaped loop of wire. It does so by approximating the loop as a set of
%infinitesimal vertical dipoles, dm, and integrating over the area of the
%loop by Gauss quadrature. It uses get_VMD_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,filterName)
%to evaluate the secondary magnetic field (the Earth response) of a unit
%dipole at a receiver a distance rRx away from the dipole. Then it divides
%up the domain (the polygon's area) into triangles, finds the Gauss
%quadrature points within each triangle and their associated weights, and
%uses Gauss quadrature using get_VMD_FD_FHT at those points with those
%weights.

 
%compute center of polygon defined by vertices
centerPoly = mean(xyPolyTx);

%now each pair of neighboring vertices on the polygon has a third point
%(the center) - forming a triangle

%number of polygon vertices 
nvertices = size(xyPolyTx,1);

%number of frequencies
nf = length(freqs);

%initialize the value of the integral to zero
I = zeros(1,nf);

%kwk debug: plot triangles:
figure;

%loop over the number of triangles = the number of polygon vertices
for k=1:nvertices
    %find the 3 vertices of the current triangle (2 polygon + center)
    if (k == nvertices)
        v = [ xyPolyTx(k,:) ; xyPolyTx(1,:) ; centerPoly ];
    else
        v = [ xyPolyTx(k,:) ; xyPolyTx(k+1,:) ; centerPoly ];
    end
    
    %compute the Gauss quadrature points and weights for this triangle
    %get Gauss points and weights on the standard triangle first
    xw = TriGaussPoints(GQorder);
    %convert to our coordinates
    for l=1:size(xw(:,1))
        x_(l) = v(1,1)*(1-xw(l,1)-xw(l,2))+v(2,1)*xw(l,1)+v(3,1)*xw(l,2);
        y_(l) = v(1,2)*(1-xw(l,1)-xw(l,2))+v(2,2)*xw(l,1)+v(3,2)*xw(l,2);
    end

    %kwk debug: plot triangles and quad points:    
     plot(v([1:end 1],1),v([1:end 1],2),'k-'); hold on
      plot(x_,y_,'*')
   
    %compute distances (horizontal range only from each integration point to receiver
    distRx_ = sqrt( (x_ - xyRx(1)).^2 + ...
                   (y_ - xyRx(2)).^2 ); 
    
    %perform numerical integration (evaluate vertical magnetic dipole at
    %each integration point and multiply by its integration weight)
    tmp = 0.0*I;
    for l=1:length(distRx_)
        tmp = tmp + get_VMD_FD_FHT(freqs,zTx,distRx_(l),zRx,sig,mu,z,filterName)*xw(l,3);
    end
    tmp = tmp*polyarea(v(:,1),v(:,2));
    I = I + tmp;
    %on to the next triangle!
end

% normalize by true polygon area:
area = polyarea(xyPolyTx(:,1),xyPolyTx(:,2));

I = I/area;

Bz = I;
 


