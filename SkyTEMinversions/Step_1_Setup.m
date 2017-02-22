% This script retrieves data from the _dat.xyz SkyTEM data file specified below under
% fileName, and loads all other setup info needed to run PT_RJMCMC

clear all
clc

fileName = '../TaylorGlacierSkyTEMdata/TaylorGlacier_dat.xyz';
Q = SkyTEMDataProc(fileName,5);
%gate times
times = Q.times;
%separate out the high mode and low mode times into separate arrays, get
%rid of NaNs from data and sd arrays, convert from rel error to abs error
k = 1;
l = 1;
for j=1:length(times)
    if( isnan(Q.dataHM(j)) == 0 )
        HighMode.times(k) = times(j);
        HighMode.data(k) = Q.dataHM(j);
        HighMode.sd(k) = Q.dataHM(j)*Q.dataErrHM(j);
        k = k + 1;
    end
    if( isnan(Q.dataLM(j)) == 0 )
        LowMode.times(l) = times(j);
        LowMode.data(l) = Q.dataLM(j);
        LowMode.sd(l) = Q.dataLM(j)*Q.dataErrLM(j);
        l = l + 1;
    end
end
% %data (high and low modes)
% HighMode.data = Q.dataHM;  
% LowMode.data = Q.dataLM;
% %data error (in correct units - not relative error)
% HighMode.sd = Q.dataHM.*Q.dataErrHM;
% LowMode.sd = Q.dataLM.*Q.dataErrLM;
%Height of transmitter and receivers (low and high modes) and
%xy-coordinates of the receiver (same datum for z as the model)
zTx     = -30; 
xyRx    = [17 0]; 
LowMode.zRx = -31.90;
HighMode.zRx = -30.19;         
% coordinates of the vertices of the transmitter loop polygon
xyPolyTx = [ -15.09 -2.00 ; 
             -8.11 -10.16 ;
              8.11 -10.16 ; 
             15.09  -2.00 ;
             15.09   2.00 ;
              8.11  10.16 ; 
             -8.11  10.16 ; 
             -15.09  2.00 ];
                 
%Parameters needed for the Bayesian inversion
numIterations = 5e5;
saveEvery = 1e4;   
log10rho_min = 0.5;  %min resistivity allowed
log10rho_max = 5;    %max resistivity allowed
zMin = 0.0;          %min interface depth allowed
zmax = 400;          %max interface depth allowed
kMax = 20;           %max number of layers allowed
kMin = 1;            %min number of layers allowed
nFreqsPerDecade = 15;   %density of Fourier domain sampling
HankelFilterName = 'kk201Hankel.txt';   %Hankel filters
CosSinFilterName = 'kk201CosSin.txt';   %Cosine-Sine filters
LoopQuadOrder = 3;      %Gauss quadrature order for loop integeration
lowPassFilters = [ 4.500E+05  3.000E+05 ];  %Butterworth low-pass filters

%low-mode transmitter ramp: [ time    fraction-of-full-power ]
LowMode.ramp = [
-1.000E-03   0.00 
-9.800E-04   0.46 
-9.550E-04   0.79 
-9.200E-04   0.94 
-9.100E-04   0.95 
-8.000E-04   1.00 
0.000E+00    1.00 
4.400E-07    0.98 
8.000E-07    0.93 
1.100E-06    0.87 
1.580E-06    0.69 
2.160E-06    0.46 
2.940E-06    0.22 
3.900E-06    0.07 
5.340E-06    0.01 
6.480E-06    0.00 
];

%high-mode transmitter ramp
HighMode.ramp = [ 
 -1.000E-02   0.00 
-8.700E-03   0.36 
-7.100E-03   0.64 
-5.500E-03   0.81 
-3.300E-03   0.93 
-1.200E-03   0.98 
-5.000E-05   1.00 
 0.000E+00   1.00 
 3.700E-06   0.94 
 1.540E-05   0.72 
 3.080E-05   0.43 
 4.220E-05   0.21 
 5.090E-05   0.01 
 5.500E-05   0.00 
];


