  
%%
% Generate synthetic data from a forward model
%
% Note for SkyTEM: This will be replaced with a routine that gets the
% SkyTEM data and writes out a data file prior to PT_RJMCMC inversion.
%

clear all
clc

%
% Model setup:
%
z   = [-1d5  0  300  350 ];   % m, vertical position of layer top boundaries, 1st one ignored since it's the top of air.
sig = [1d-12  1/5000 1/10 1/10 1/500 ];   % (S/m), layer conductvities. First is air (1d-12).

%
% SkyTEM system setup:
%
times = logspace(-5,-2,18);  % s, gate times.

rTx = 12;    % m, trasmitter loop radius
zTx = -35;   % m, transmitter vertical position (use same datum as model z).
zRx = -35.5; % m, reciever coil vertical position (use same datum as model z).
rRx = 12;    % m, receiver coil distance from center of loop. For SkyTEM this should be same as loop radius rTx.

%
% Inversion setup:
%
noiseToAdd = 0.05;  % relative noise to add to model response when creating synthetic noisy data (0.05 = 5% relative noise).
kMin = 1;           % Minimum number of additional layer interfaces to add
kMax = 35;          % Maximum number of additional layer interfaces to add
nFixedLayers = 1;   % Number of fixed layers in z. 1 means only top layer is fixed, which is the air layer in z,sig above.
minThickness = 0;   % Inversion layers will be fixed to be within these thickness bounds.
maxThickness = 600;   
log10rho_min = -1;  % min log10(resistivity). You shouldn't need to change these.
log10rho_max = 5;   % max log10(resistivity). You shouldn't need to change these.

numIterations = 4e5;  % 100,000 Number of RJ-MCMC iterations to carry out (should be at least 100,000)
saveEvery     = 1e4;   % 10,000 RJ-MCMC models from each parallel tempering chain are save every N iterations (should be at least 10000)


ReferenceModel.z   = z;     % Reference model used in plotting code. Store the "truth" here, or the LCI inversion model if you want to overlay it on the Bayesian result
ReferenceModel.rho = 1./sig;
 
%
% Name of synthetic data file to output:
%
DataFile = 'Test_SynthData.mat'; 

outputFolder = 'SynthTrash';  % Name of folder to store PT_RJMCMC results


%--------------------------------------------------------------------------
%-------- You shouldn't need to change anything below here ----------------
%--------------------------------------------------------------------------
%
% Some numerial parameters used by the forward code (don't change these):
%
nFreqsPerDecade = 10;
HankelFilterName = 'kk101Hankel.txt';
CosSinFilterName = 'kk101CosSin.txt';
mu = ones(size(sig)); % relative magnetic permeability of each layer. Set this to 1 always!


% 
% Compute the model response, which is normalized by loop dipole moment, so
% Bz has units T/(Am^2) = Vs/(Am^4)
% Note Bz returned below is actually dBz/dt, so the units are V/(Am^4)
%

rampTime = [0 1; 5d-6 0.5;  10d-6 0]; 
tic
Bz = get_LoopFields_TD_FHT_Synth(times,rTx,zTx,rRx,zRx,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade,rampTime);    
toc                      
% For comparison, make another model with no conductive layers:                     
sig_resistor = [1d-12  1/1000 1/1000 1/1000 1/1000];      
Bz_resistor  = get_LoopFields_TD_FHT(times,rTx,zTx,rRx,zRx,sig_resistor,mu,z,HankelFilterName,CosSinFilterName,nFreqsPerDecade,rampTime);            

%
% Add noise to create synthetic noisy data:
%
noise = randn(size(Bz))*noiseToAdd.*abs(Bz);
data  = Bz + noise;
sd    = noiseToAdd.*abs(Bz);  
Ramp = rampTime;

%
% Plot SkyTEM responses:                     
%
figure;
subplot(1,2,1);
zz = zeros(2*length(z),1);
zz(1:2:end-1)   = z;
zz(2:2:end-2) = z(2:end);
zz(1) = zz(2)-50;
zz(end) = zz(end-1)+50;

rr = 0*zz;
rr(1:2:end-1)   = 1./sig;
rr(2:2:end-2)   = 1./sig(1:end-1);
rr(end)         = 1./sig(end);

rr2 = 0*zz;
rr2(1:2:end-1)   = 1./sig_resistor;
rr2(2:2:end-2)   = 1./sig_resistor(1:end-1); 
rr2(end)         = 1./sig_resistor(end);

semilogx(rr,zz,'b-','linewidth',2);
hold on
semilogx(rr2,zz,'r-','linewidth',2);
axis ij;
xlabel('Log10(resistivity')
ylabel('Depth (m)')
  set(gca,'fontsize',14)  
legend('Conductive Model','Resistive Model')
set(gca,'xlim',10.^[log10rho_min log10rho_max])

subplot(1,2,2);
loglog(times,abs(Bz),'b-','linewidth',1); 
hold on;
loglog(times,abs(Bz_resistor),'r-','linewidth',1); 
loglog(times,abs(data),'bo'); 
loglog([times; times],[abs(data)+sd abs(data)-sd;]','b-','linewidth',2); % error bar plotting
ylabel('dBz/dt (V/Am^4)')
xlabel('Time (s)')
set(gca,'fontsize',14)
legend('Conductive Model', 'Resistive Model', 'Synthetic Noisy Data')


%--------------------------------------------------------------------------
%% Create input data file needed by RJMCMC code:
%

% Fixed model layers:
rho = 1/sig(1:nFixedLayers);   
z   = z(1:nFixedLayers+1);    % Layer boundaries for fixed layers

zMin = z(nFixedLayers+1)+minThickness; % minimum inversion layer position. z(2)+2 means inversion layers added below z(2) at 2 m minimum thickness.
zMax = maxThickness;

if saveEvery > numIterations
    saveEvery = numIterations;
end

% Save config to .mat file:
save(DataFile,'data','sd','times','rTx','zTx','rRx','zRx', ...
              'HankelFilterName','CosSinFilterName','nFreqsPerDecade',...
              'rho','z','kMin','kMax','zMin','zMax','log10rho_min',...
              'log10rho_max','ReferenceModel','numIterations','saveEvery','Ramp');



 
 
 