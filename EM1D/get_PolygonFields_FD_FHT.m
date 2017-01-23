function Bz = get_PolygonFields_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,filterName,vertices,GQorder,A)

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

%receiver location is the last row of vertices
r = vertices(end,:);
%compute center of polygon defined by vertices
centerPoly = mean(vertices);
%now each pair of neighboring vertices on the polygon has a third point
%(the center) - forming a triangle

%number of polygon vertices + 1 (the receiver location is included last)
nvertices = size(vertices,1);
%initialize the value of the integral to zero
nf = length(freqs);
I = zeros(1,nf);
%loop over the number of triangles = the number of polygon vertices
for k=1:nvertices-1
    %find the 3 vertices of the current triangle (2 polygon + center)
    v = [ vertices(k,:) ; vertices(k+1,:) ; centerPoly ];
    [X,Y,Wx,Wy]=triquad(GQorder,v);
    %compute distance from each integration point to the receiver
    Rx = ones(GQorder)*r(1,1);
    Ry = ones(GQorder)*r(1,2);
    Rz = ones(GQorder)*(zRx - zTx);
    R = ((X-Rx).^2 + (Y-Ry).^2 + (Rz).^2).^(0.5);
    %evaluate the vertical magnetic dipole at each integration point
    feval = zeros(size(R,1),size(R,2),nf);
    %keyboard
    for l=1:size(R,1)
        for m=1:size(R,2)
            feval(l,m,:) = get_VMD_FD_FHT(freqs,zTx,R(l,m),zRx,sig,mu,z,filterName);
            %keyboard
            %feval(l,m,:) = tmp;
        end
    end
    %compute the integral over this triangle (finally!)
    %keyboard
    for j=1:nf
        I(j) = I(j) + Wx'*feval(:,:,j)*Wy;
    end   
end
%normalize by the area of the loop
I = I/A;
keyboard
Bz = get_VMD_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,filterName);
keyboard
h = 70;


