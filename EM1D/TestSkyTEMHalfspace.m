  
%%
% Generate synthetic data from a forward model
%
% Note for SkyTEM: This will be replaced with a routine that gets the
% SkyTEM data and writes out a data file prior to PT_RJMCMC inversion.
%

%
% Model setup:
%
z   = [-1d5     0   100  ];   % m, vertical position of layer top boundaries, 1st one ignored since it's the top of air.
sig = [1d-12   1/10 1/10 ];   % (S/m), layer conductvities. First is air (1d-12).

%
% SkyTEM system setup:
%

lowPassFiltersLM = [ 4.500E+05  3.000E+05 ];
lowPassFiltersHM = [ 4.500E+05  3.000E+05 ];


tLM = [ 1.300E-05
 1.650E-05
 2.100E-05
 2.700E-05
 3.450E-05
 4.350E-05
 5.500E-05
 6.950E-05
 8.800E-05
 1.115E-04
 1.405E-04
 1.775E-04
 2.240E-04
 2.820E-04
 3.555E-04
 4.480E-04
 5.645E-04
 7.110E-04
 8.955E-04
];

tHM = [ 8.570E-05
 1.092E-04
 1.382E-04
 1.752E-04
 2.217E-04
 2.797E-04
 3.532E-04
 4.457E-04
 5.622E-04
 7.087E-04
 8.932E-04
 1.126E-03
 1.418E-03
 1.786E-03
 2.250E-03
 2.833E-03
 3.568E-03
 4.493E-03
 5.808E-03
];

 rampLM = [
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

rampHM=[ 
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


% 10 ohm-m
dLM=[ 5.738414E-06
 4.236351E-06
 3.111256E-06
 2.243019E-06
 1.608307E-06
 1.161711E-06
 8.251093E-07
 5.775482E-07
 3.972306E-07
 2.682531E-07
 1.798862E-07
 1.183093E-07
 7.655095E-08
 4.898255E-08
 3.071389E-08
 1.892388E-08
 1.144417E-08
 6.789652E-09
 3.946301E-09
];

dHM = [
8.457780E-07
4.794991E-07
2.861063E-07
1.729966E-07
1.053307E-07
6.465972E-08
3.950176E-08
2.406247E-08
1.459279E-08
8.806743E-09
5.283412E-09
3.147092E-09
1.863403E-09
1.094271E-09
6.370838E-10
3.673768E-10
2.096507E-10
1.181879E-10
6.198442E-11
];




xyPolyTx = [ -15.09 -2.00 ; 
             -8.11 -10.16 ;
              8.11 -10.16 ; 
             15.09  -2.00 ;
             15.09   2.00 ;
              8.11  10.16 ; 
             -8.11  10.16 ; 
             -15.09  2.00 ];
           
zTx     = -30; 
xyRx    = [17 0]; 

zRxLM = -31.90;
zRxHM = -30.19;         % m, reciever coil vertical position (use same datum as model z).
 
% Normalize SkyTEM solutions by loop area:
area = polyarea(xyPolyTx(1:end,1),xyPolyTx(1:end,2) );
dHM = dHM/area;
dLM = dLM/area;

%--------------------------------------------------------------------------
%-------- You shouldn't need to change anything below here ----------------
%--------------------------------------------------------------------------
%
% Some numerial parameters used by the forward code (don't change these):
%
nFreqsPerDecade = 15;
LoopQuadOrder   = 5;   % Gauss quadrature order for polygon loop integration
HankelFilterName = 'kk201Hankel.txt';
CosSinFilterName = 'kk201CosSin.txt';

mu = ones(size(sig)); % relative magnetic permeability of each layer. Set this to 1 always!

% 
% Compute the model response, which is normalized by loop dipole moment, so
% Bz has units T/(Am^2) = Vs/(Am^4)
% Note Bz returned below is actually dBz/dt, so the units are V/(Am^4)
%

tic
BzLM = get_LoopFields_TD_FHT(tLM,xyPolyTx,zTx,xyRx,zRxLM,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade,LoopQuadOrder,rampLM,lowPassFiltersLM);    
toc                     
% BzLM_NoLowPass = get_LoopFields_TD_FHT(tLM,xyPolyTx,zTx,xyRx,zRxLM,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade,LoopQuadOrder,rampLM);  
% BzLM_NoRamp    = get_LoopFields_TD_FHT(tLM,xyPolyTx,zTx,xyRx,zRxLM,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade,LoopQuadOrder,[],lowPassFiltersLM);  

tic
BzHM = get_LoopFields_TD_FHT(tHM,xyPolyTx,zTx,xyRx,zRxHM,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade,LoopQuadOrder,rampHM,lowPassFiltersHM);    
toc     

% 
% Plot SkyTEM responses:                     
%
lw = 2;
figure;
subplot(2,1,1)
loglog(tLM,abs(BzLM),'b-','linewidth',1); 
hold on;
loglog(tLM,abs(dLM),'c--','linewidth',lw); 

% loglog(tLM,abs(BzLM_NoLowPass),'m-','linewidth',1); 
% loglog(tLM,abs(BzLM_NoRamp),'g-','linewidth',1); 

loglog(tHM,abs(BzHM),'r-','linewidth',1); 
hold on;
loglog(tHM,abs(dHM),'k--','linewidth',lw);

 

ylabel('dBz/dt (V/Am^4)')
xlabel('Time (s)')
set(gca,'fontsize',14)
legend('LM Scripps', 'LM SkyTEM','HM Scripps', 'HM SkyTEM')
%legend('LM Scripps', 'LM SkyTEM','LM Scripps No LowPass','LM Scripps No Ramp','HM Scripps', 'HM SkyTEM') 

subplot(2,1,2)
loglog(tLM,100*abs(BzLM+dLM)./abs(dLM),'b-','linewidth',lw)
hold on
loglog(tHM,100*abs(BzHM+dHM)./abs(dHM),'r-','linewidth',lw)
xlabel('Time (s)')
ylabel('Relative Difference %')
legend('LM', 'HM')
 set(gca,'fontsize',14)
 title('|Scripps - SkyTEM | / |SkyTEM|')
 
 
%% Compare with Analytical formula central loop sounding on surface of a halfspace:
clear all

% For testing with halfpsace equation:

% % % formula from Christiansen et el groundwater chapter 
% % %eqn 6.26:


nFreqsPerDecade = 15;
LoopQuadOrder   = 3;   % Gauss quadrature order for polygon loop integration
loopRadius      = 15; % 10 m circular loop for Christiansen formula

nPolygonSides   = 20;  % number of polygon sides. more means better circle approximation. must be >2

HankelFilterName = 'kk201Hankel.txt';
CosSinFilterName = 'kk201CosSin.txt';

z   = [-1d5     0     ];   % m, vertical position of layer top boundaries, 1st one ignored since it's the top of air.
sig = [1d-20   1/10   ]; 


%-------------------------------
% Don't change anything below here:

tlong       = logspace(-6,-1,40);
xyRx        = [0 0];  % has to be 0 since Christiansen formula is only for receiver coil in center of transmitter loop
zRx         = 0; % Rx depth. has to be 0 since Christiansen formula (from Ward and Hohmann actually) is valid only on surface of halspace
zTx         = 0; % Tx depth. has to be 0 since Christiansen formula (from Ward and Hohmann actually) is valid only on surface of halspace
mu          = ones(size(sig)); % relative magnetic permeability of each layer. Set this to 1 always!

% Make a polygon approximation of a circle to test get_PolygonFields versus
% Christiansen et al analytic formula for perfect circle
 
L = -pi/2+linspace(0,2*pi,nPolygonSides+1);
L = L(1:end-1);
xyPolyTx(:,1) = loopRadius*cos(L)';
xyPolyTx(:,2) = loopRadius*sin(L)';
 
% Nothing to change below here:
sig0 = sig(2);
theta = sqrt(4*pi*1d-7*sig0 ./ (4*tlong) );
a     = loopRadius;  
dbzdt = -1/(sig0*a^3)*(3*erf(theta*a) - 2/sqrt(pi)*theta.*a.*(3+2*theta.^2*a^2).*exp(-theta.^2*a^2));
dbzdt = dbzdt/(pi*a^2);


% Compute loop response for the polygon, without any waveform ramp or
% filters:
BzPoly = get_LoopFields_TD_FHT(tlong,xyPolyTx,zTx,xyRx,zRx,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade,...
                             LoopQuadOrder);    


figure;
subplot(2,1,1)
loglog(tlong,-(dbzdt),'k-','linewidth',1);
hold on;
loglog(tlong,-(BzPoly),'r--','linewidth',1);
legend('Analytic','Polygon')
xlabel('Time (s)')
title(sprintf('LoopQuadOrder=%i, nPolygonSides=%i, nFreqsPerDecade=%i',LoopQuadOrder,nPolygonSides,nFreqsPerDecade))
set(gca,'fontsize',14)
 
subplot(2,1,2)
loglog(tlong,100*abs(dbzdt'-BzPoly)./abs(dbzdt)','-','linewidth',1);
hold all
ylabel('Relative Difference (%)')
xlabel('Time (s)')
set(gca,'fontsize',14)


%% Test FD Polygon integration without time domain:
clear all;


xyPolyTx = [ -15.09 -2.00 ; 
             -8.11 -10.16 ;
              8.11 -10.16 ; 
             15.09  -2.00 ;
             15.09   2.00 ;
              8.11  10.16 ; 
             -8.11  10.16 ; 
             -15.09  2.00 ];
           
zTx     = -30; 
zRx     = -31.90;
xyRx    = [17 0]; 

nFreqsPerDecade = 3;
 
HankelFilterName = 'kk201Hankel.txt';
 

z   = [-1d5     0     ];   % m, vertical position of layer top boundaries, 1st one ignored since it's the top of air.
sig = [1d-20   1/10   ]; 
mu  = ones(size(sig)); % relative magnetic permeability of each layer. Set this to 1 always!


freqs = 10.^(-3:1/nFreqsPerDecade:6); 
 
tic
nOrderMax = 8;
for order = 1:nOrderMax
    BzPoly(order,:) = get_PolygonFields_HED_FD_FHT(freqs,xyPolyTx,zTx,xyRx,zRx,sig,mu,z, HankelFilterName, order);
   % BzPoly(order,:) = get_PolygonFields_VMD_FD_FHT(freqs,xyPolyTx,zTx,xyRx,zRx,sig,mu,z, HankelFilterName, order);
end
toc
 
figure;
for order = 1:nOrderMax-1
    loglog(freqs,100*abs(BzPoly(order,:)-BzPoly(nOrderMax,:))./abs(BzPoly(nOrderMax,:)),'-','linewidth',1);
    hold all
end
legend(gca,num2str([1:nOrderMax-1]'))
title(sprintf(' Rx range=%.0f',xyRx(1)))
ylabel('Relative Difference (%)')
xlabel('Frequency (Hz)')
set(gca,'fontsize',14)