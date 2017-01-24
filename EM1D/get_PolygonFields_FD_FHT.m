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

nf = length(freqs);
%initialize the value of the integral to zero
%method 1
I = zeros(1,nf);
%method 2
I2 = I;

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
    %method 1
    [X,Y,Wx,Wy]=triquad(GQorder,v);
    %method 2
    [x,w] = simplexquad(GQorder,v);

    %kwk debug: plot triangles and quad points:    
     plot(v([1:end 1],1),v([1:end 1],2),'k-'); hold on
     plot(x(:,1),x(:,2),'o'); hold on
      plot(X,Y,'+'); hold on
  
    
% KWK debug R is just horizontal range to receiver!
% Rz = ones(GQorder)*(zRx - zTx); 
   % R = sqrt((X-xyRx(1)).^2 + (Y-xyRx(2)).^2);  %+ (Rz).^2).^(0.5);
    
    %compute vector of distances (for I2)
    distRx = sqrt( (x(:,1) - xyRx(1)).^2 + ...
                   (x(:,2) - xyRx(2)).^2 ); 
 
    
    %evaluate the vertical magnetic dipole at each integration point
    %method 2
    for l=1:length(distRx)
        I2 = I2 + get_VMD_FD_FHT(freqs,zTx,distRx(l),zRx,sig,mu,z,filterName)*w(l);
    end
    %method 1
%     feval = zeros(size(R,1),size(R,2),nf);
%     for l=1:size(R,1)
%         for m=1:size(R,2)
%             feval(l,m,:) = get_VMD_FD_FHT(freqs,zTx,R(l,m),zRx,sig,mu,z,filterName);
%         end
%     end
%     %compute the integral over this triangle (finally!)
%     for j=1:nf
%         I(j) = I(j) + Wx'*feval(:,:,j)*Wy;
%     end
    %on to the next triangle!
end

% normalize by true polygon area:
area = polyarea(xyPolyTx(:,1),xyPolyTx(:,2));

%I = I/area;
I2 = I2/area;
% BzCircle = get_VMD_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,filterName);
Bz = I2;
 



%keyboard


