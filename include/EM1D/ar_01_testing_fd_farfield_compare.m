
clear all
clc

% frequencies
nFreqsPerDecade     = 10;
freqLowLimit        = 1d-3; 
freqHighLimit       = 1d6;
freqs = 10.^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit));

% model
z   = [-1d5    0   100  ];   
sig = [1d-12   1/10 1/10 ];
% sig = [1d-12   1d-12 1d-12 ];
% sig = [1 1 1 ];

mu = ones(size(sig));

% geometry
zTx     = -30; 
zRx     =  -32;
xyRx    = [40 0]; 
rTx     = 17;
nangles = 200;
farfieldapprox = true;

angles  = linspace(0,2*pi,nangles);
xTx     = rTx*cos(angles(1:end-1));
yTx     = rTx*sin(angles(1:end-1));
xyPolyTx = [xTx' yTx'];
% xyPolyTx = [ -15.09 -2.00 ; 
%              -8.11 -10.16 ;
%               8.11 -10.16 ; 
%              15.09  -2.00 ;
%              15.09   2.00 ;
%               8.11  10.16 ; 
%              -8.11  10.16 ; 
%              -15.09  2.00 ];

HankelFilterName = 'kk201Hankel.txt';    
%HankelFilterName = 'wer201Hankel.txt';    
GQorder                = 5;   % number of gauss quadrature points per Tx wire segment for line integration
nSplinePointsPerDecade = 10;  % number of spline points per decade of ranges for spline interpolation of Bz(r)


% For testing far field solution, make the transmitter polygon tiny so that
% it looks like a point dipole at the receiver location:
if farfieldapprox
    xyPolyTx = xyPolyTx/1000;
end    
% Then the fields should be nearly the same (<< 1% difference). 

% Polyon transmitter:
tic
BzFDPoly = get_PolygonFields_HED_FD_FHT(freqs,xyPolyTx,zTx,xyRx,zRx,sig,mu,z,HankelFilterName,GQorder,nSplinePointsPerDecade);

t = toc;
fprintf('%32s %.3f s\n','get_PolygonFields_FD_FHT time:',t);

% Perfectly circular loop:  
area     = polyarea(xyPolyTx(:,1),xyPolyTx(:,2));   % circle approximation requires area of polygon
rTxLoop  = sqrt(area/pi);                           % from area = pi*r^2
rRx      = norm(xyRx);                              % receiver is at same range from center of polygon 

tic
BzFDCirc = get_LoopFields_circle_FD_FHT(freqs,rTxLoop,zTx,rRx,zRx,sig,mu,z,HankelFilterName);
t = toc;
fprintf('%32s %.3f s\n','get_LoopFields_FD_FHT time:',t);

% Point dipole:
tic
BzPoint = get_VMD_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,HankelFilterName);
t = toc;
fprintf('%32s %.3f s\n','get_VMD_FD_FHT time:',t);
% Note that when rRx >> rTxLoop, then the point dipole field and loop and
% polygon fields should look the same

figure;
subplot(2,1,1)
loglog(freqs,abs(BzFDPoly),'b-',freqs,abs(BzFDCirc),'ro',freqs,abs(BzPoint),'g+')
legend('Polygon','Circle','Point')
xlabel('Frequency (Hz)')
ylabel('Amplitude')

subplot(2,1,2)
loglog(freqs,abs(BzFDCirc-BzFDPoly)./abs(BzFDCirc)*100,'b-')
hold on;
loglog(freqs,abs(BzFDCirc-BzPoint)./abs(BzPoint)*100,'r-')
loglog(freqs,abs(BzFDPoly-BzPoint)./abs(BzPoint)*100,'g-')
xlabel('Frequency (Hz)')
ylabel('Relative Difference (%)')
legend('(Polygon-Circle)/Circle','(Circle-Point)/Point','(BzFDPoly-Point)/Point')

 
