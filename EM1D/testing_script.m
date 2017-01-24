
clear all
clc

nFreqsPerDecade = 10;
freqLowLimit        = 1d-3; 
freqHighLimit       = 1d6;
freqs = 10.^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit));

z   = [-1d5     0   100  ];   
sig = [1d-20   1/10 1/10 ];
mu = ones(size(sig));

zTx   = -30; 
rRx   = 2+sqrt(488/pi);
zRx = -31.90;
A = 488; %area of polygon in square meters

HankelFilterName = 'kk201Hankel.txt';

%the last vertex is actually the receiver location
vertices = [ -15.09 -2.00 ; -8.11 -10.16 ; 8.11 -10.16 ; 15.09 -2.00 ; 15.09 2.00 ;...
        8.11 10.16 ; -8.11 10.16 ; -15.09 2.00 ; -17.00 0.00 ];
    
GQorder = 2;

%Bz = get_VMD_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,HankelFilterName);

tic
BzFD = get_PolygonFields_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,HankelFilterName,vertices,GQorder,A);
toc

%norm(Bz-BzFD)
