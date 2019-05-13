
clear all
clc

nFreqsPerDecade     = 10;
freqLowLimit        = 1d-3; 
freqHighLimit       = 1d6;
freqs = 10.^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit));


z   = [-1d5    0   100  ];   
sig = [1d-12   1/10 1/10 ];
mu = ones(size(sig));

zTx     = -30; 
zRx     = -31.90;
xyRx    = [17 0]; 

xyPolyTx = [ -15.09 -2.00 ; 
             -8.11 -10.16 ;
              8.11 -10.16 ; 
             15.09  -2.00 ;
             15.09   2.00 ;
              8.11  10.16 ; 
             -8.11  10.16 ; 
             -15.09  2.00 ];
           
HankelFilterName = 'kk201Hankel.txt';    
GQorder        = 2;


% For testing far field solution, make the transmitter polygon tiny so that
% it looks like a point dipole at the receiver location:
%xyPolyTx    = xyPolyTx/100;    
% Then the fields should be nearly the same (<< 1% difference). 

% Polyon transmitter:
tic
BzFDPoly = get_PolygonFields_HED_FD_FHT(freqs,xyPolyTx,zTx,xyRx,zRx,sig,mu,z,HankelFilterName,GQorder);
t = toc;
fprintf('%32s %.3f s\n','get_PolygonFields_FD_FHT time:',t);

% Perfectly circular loop:  
area     = polyarea(xyPolyTx(:,1),xyPolyTx(:,2));   % circle approximation requires area of polygon
rTxLoop  = sqrt(area/pi);                           % from area = pi*r^2
rRx      = norm(xyRx);                              % receiver is at same range from center of polygon 

tic
BzFDCirc = get_LoopFields_FD_FHT(freqs,rTxLoop,zTx,rRx,zRx,sig,mu,z,HankelFilterName);
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

 
