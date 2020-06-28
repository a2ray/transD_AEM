clear
addpath('../RJMCMC')
%% model
S.z    = [-1d5, 0];
S.zMax = 100;
S.rho = [1d12];
x.z    = [20 50];
x.rhoh = [1 0 2];
x.rhov = [1 0 2];
writeData = false; % won't overwrite

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
S.HM_times = [
            7.53900E-05 9.60000E-05
            9.63900E-05 1.22000E-04
            1.22390E-04 1.54000E-04
            1.54390E-04 1.96000E-04
            1.96390E-04 2.47000E-04
            2.47390E-04 3.12000E-04
            3.12390E-04 3.94000E-04
            3.94390E-04 4.97000E-04
            4.97390E-04 6.27000E-04
            6.27390E-04 7.90000E-04
            7.90390E-04 9.96000E-04
            9.96390E-04 1.25500E-03
            1.25539E-03 1.58100E-03
            1.58139E-03 1.99100E-03
            1.99139E-03 2.50800E-03
            2.50839E-03 3.15800E-03
            3.15839E-03 3.97700E-03
            3.97739E-03 5.00800E-03
            5.00839E-03 6.30600E-03
            6.30639E-03 7.93900E-03
            7.93939E-03 9.73900E-03
        ];
S.HM_times = exp(mean(log(S.HM_times),2));
S.HM_ramp = [
            -1.000E-02 0.000E+00
            -8.386E-03 4.568E-01
            -6.380E-03 7.526E-01
            -3.783E-03 9.204E-01
            0.000E+00 1.000E+00
            3.960E-07 9.984E-01
            7.782E-07 9.914E-01
            1.212E-06 9.799E-01
            3.440E-06 9.175E-01
            1.981E-05 4.587E-01
            3.619E-05 7.675E-03
            3.664E-05 3.072E-03
            3.719E-05 8.319E-04
            3.798E-05 1.190E-04
            3.997E-05 0.000E+00
            1.000E-02 0.000E+00
            ];

S.HM_zRx = -37.0; % m, reciever coil vertical position (use same datum as model z).
%% for both HM and LM in common
area     = polyarea(S.xyPolyTx(:,1),S.xyPolyTx(:,2));   % circle approximation requires area of polygon
S.rTxLoop  = sqrt(area/pi);                           % from area = pi*r^2
S.rRx      = norm(S.xyRx);
S.zTx = -35;   % m, transmitter vertical position (use same datum as model z).
S.lowPassFilters = [300000 450000];

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
 
% model loop
S.modelLMpoly = false;
S.modelHMpoly = false;
S.modelLMloop = false;
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
loglog(S.HM_times,abs(BzLoop),'r*','linewidth',1); 
grid on
ylabel('dBz/dt (V/Am^4)')
xlabel('Time (s)')
set(gca,'fontsize',14)


 
 