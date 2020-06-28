clear
%% SkyTEM model and data
S.z    = [-1d6, 0];
S.zMax = 200;
S.rho = [1d12];
x.z    = [100 ];
x.rhoh = [1 1];
x.rhov = [1 1];
writeData = false; % won't overwrite

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

%% SkyTEM polygon and receiver xy details
S.xyPolyTx = [ -15.09 -2.00 ; 
             -8.11 -10.16 ;
              8.11 -10.16 ; 
             15.09  -2.00 ;
             15.09   2.00 ;
              8.11  10.16 ; 
             -8.11  10.16 ; 
             -15.09  2.00 ];

 S.xyRx    = [17 0]; 

%% SkyTEM HM system setup:
%
S.HM_times = [ 8.570E-05
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

S.HM_ramp = [ 
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

S.HM_zRx = -30.19; % m, reciever coil vertical position (use same datum as model z).

%% SkyTEM LM setup:
S.LM_ramp=[
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
        
S.LM_times=[ 1.300E-05
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

S.LM_zRx = -31.9; % not sure if Rx should be same for HM and LM ...        
%% for both HM and LM in common
area     = polyarea(S.xyPolyTx(:,1),S.xyPolyTx(:,2));   % circle approximation requires area of polygon
S.rTxLoop  = sqrt(area/pi);                           % from area = pi*r^2
S.rRx      = norm(S.xyRx);
S.zTx = -30;   % m, transmitter vertical position (use same datum as model z).
S.lowPassFilters = [300000 450000];
dHM = dHM/area;
dLM = dLM/area;

%% --------------------------------------------------------------------------
%-------- You shouldn't need to change anything below here ----------------
%--------------------------------------------------------------------------
%
% Some numerial parameters used by the forward code (don't change these):
%
S.nFreqsPerDecade = 10;
S.HankelFilterName = 'kk201Hankel.txt';
S.CosSinFilterName = 'kk201CosSin.txt';
S.LoopQuadOrder = 3;

% 
% Compute the model response, which is normalized by loop dipole moment, so
% Bz has units T/(Am^2) = Vs/(Am^4)
% Note Bz returned below is actually dBz/dt, so the units are V/(Am^4)
%
 
% model polygon
S.modelLMpoly = true;
S.modelHMpoly = true;
S.modelLMloop = false;
S.modelHMloop = false;
tic
BzPoly = get_field(S,x);
toc                      

% model loop
S.modelLMpoly = false;
S.modelHMpoly = false;
S.modelLMloop = true;
S.modelHMloop = true;
tic
BzLoop = get_field(S,x);
toc                      

%
% Plot SkyTEM responses:                     
%
figure;
subplot(1,2,1);
plot_model(S,x,2','b')
set(gca,'ydir','rev')
ylabel('Depth m')
set(gca,'fontsize',14)
grid on
subplot(1,2,2);
loglog([S.LM_times;S.HM_times],abs(BzPoly),'b*','linewidth',1); hold all
loglog([S.LM_times;S.HM_times],abs(BzLoop),'r*','linewidth',1); 
loglog(S.LM_times,abs(dLM),'linewidth',1); 
loglog(S.HM_times,abs(dHM),'linewidth',1); 

grid on
ylabel('dBz/dt (V/Am^4)')
xlabel('Time (s)')
set(gca,'fontsize',14)
legend('Polygon', 'Loop','SkyTEM LM', 'SkyTEM HM')
set(gcf, 'Units','pixels', 'Position',[0 0 1200 600])
figure()
semilogx([S.LM_times;S.HM_times], 100*abs(abs(BzPoly)-abs(BzLoop))./(abs(BzPoly)),'*')
ylabel('% diff')
xlabel('Time (s)')
title('polygon - loop response')
set(gca,'fontsize',14)
%--------------------------------------------------------------------------


 
 