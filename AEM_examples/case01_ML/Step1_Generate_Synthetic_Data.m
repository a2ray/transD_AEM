%% well to model
rng(2)
M = jitteryx;
writeData = true; % won't overwrite
S.ReferenceModel = M;
%%
% Generate synthetic data from a forward model
%
% Note for SkyTEM: This will be replaced with a routine that gets the
% SkyTEM data and writes out a data file prior to PT_RJMCMC inversion.
%

%
% SkyTEM system setup:
%
looparea = 337.04;
S.rTxLoop = sqrt(looparea/pi);    % m, trasmitter loop radius
S.zTx = -35;   % m, transmitter vertical position (use same datum as model z).
S.zRx = -37; % m, reciever coil vertical position (use same datum as model z).
S.rRx = 13.37;    % m, receiver coil distance from center of loop. For SkyTEM this should be same as loop radius rTx.
% Ramp settings
remove_LM = 3; remove_HM = 0;
t_HM = load('HM_times.txt');
t_HM = t_HM(remove_HM+1:end);
t_LM = load('LM_times.txt');
t_LM = t_LM(remove_LM+1:end);
ramp_HM = load('HM_ramp.txt');
ramp_LM = load('LM_ramp.txt');
S.modelLMpoly = false;
S.modelHMpoly = false;
S.modelLMloop = true;
S.modelHMloop = true;
S.LM_times = t_LM;
S.HM_times = t_HM;
S.LM_ramp  = ramp_LM;
S.HM_ramp  = ramp_HM;
S.lowPassFilters = [155e3 250e3];
%
% Inversion setup:
%
noiseToAdd = 0.05;  % relative noise to add to model response when creating synthetic noisy data (0.05 = 5% relative noise).
S.kMin = 1;           % Minimum number of additional layer interfaces to add
S.kMax = 50;          % Maximum number of additional layer interfaces to add
S.zMin = 0;   % There is no min thickness
S.zMax = M.z(end); % depth at which to place last interface  
S.log10rho_min = -0.5;  % min log10(resistivity). You shouldn't need to change these.
S.log10rho_max = 2.5;   % max log10(resistivity). You shouldn't need to change these.
S.nTemps = 4;         % number of chains in PT[H,V,intfcCount,meanH,meanV,medianH,medianV,kOut] = plot_RJMCMC(4,s_ll,k_ll ,nBurnIn,nThin,5,50,'isotropic','NOnormalize','logPDF',S);
S.Tmax = 2.0;         % Tmax for PT
S.numIterations = 250000;  % 100,000 Number of RJ-MCMC iterations to carry out (should be at least 100,000)
S.saveEvery     = 1000;   % 10,000 RJ-MCMC models from each parallel tempering chain are save every N iterations (should be at least 10000)
S.jeffereys_prior = true;
S.debug_prior = false;
S.MLerroradjust = true;
% MCMC proposal sizes
%step sizes in model space, decreasing temperature
S.UstepSize = linspace(0.04, 0.005, S.nTemps); %for update in log10 rho
S.BstepSize = linspace(0.6, 0.4, S.nTemps); %for birth / death in log10 rho
S.MstepSize = linspace(10, 2, S.nTemps); %for move interface in m
 
%
% Name of synthetic data file to output:
%
DataFile = 'Test_SynthData.mat'; 

outputFolder = 'Test_SynthData';  % Name of folder to store PT_RJMCMC results
% plot the model
figure;
plot_model(S,M,2','r')
plot_model(S,x,2','b')
plot_model(S,bg,2,'k')
xlim([-S.log10rho_max -S.log10rho_min])
set(gca,'ydir','rev')
xlabel('Log10 S/m')
ylabel('Depth m')
set(gca,'fontsize',14)
grid on

%--------------------------------------------------------------------------
%-------- You shouldn't need to change anything below here ----------------
%--------------------------------------------------------------------------
%
% Some numerial parameters used by the forward code (don't change these):
%
S.nFreqsPerDecade = 10;
S.HankelFilterName = 'kk201Hankel.txt';
S.CosSinFilterName = 'kk201CosSin.txt';

% 
% Compute the model response, which is normalized by loop dipole moment, so
% Bz has units T/(Am^2) = Vs/(Am^4)
% Note Bz returned below is actually dBz/dt, so the units are V/(Am^4)
%
 
tic
S.data = []; %needed for simple format misfit
Bz = get_field(S,M);
Bzclean = get_field(S,x);
BzBg = get_field(S,bg);
toc

%
% Add noise to create synthetic noisy data:
%
noise = randn(size(Bz))*noiseToAdd.*abs(Bz);
data  = abs(Bz) + noise;
nLM = length(t_LM);
S.d_LM = data(1:nLM);
S.d_HM = data(nLM+1:end);
S.sd_LM    = noiseToAdd.*abs(Bz(1:nLM));
S.sd_HM    = noiseToAdd.*abs(Bz(nLM+1:end));
misfit = getMisfit(M, S)
sprintf('chi2/2: %f rms: %4.2f', misfit(1), misfit(2))
%%
% Plot SkyTEM responses:                     
%
figure;
subplot(1,2,1);
plot_model(S,M,2','r')
plot_model(S,x,2','b')
plot_model(S,bg,2,'k')
xlim([-S.log10rho_max -S.log10rho_min])
set(gca,'ydir','rev')
xlabel('Log10 S/m')
ylabel('Depth m')
set(gca,'fontsize',14)
grid on
subplot(1,2,2);
loglog(S.LM_times,abs(Bzclean(1:nLM)),'b--','linewidth',1); 
hold on
loglog(S.HM_times,abs(Bzclean(nLM+1:end)),'b--','linewidth',1, 'HandleVisibility','off'); 
loglog(S.LM_times,abs(BzBg(1:nLM)),'k-','linewidth',1);
loglog(S.HM_times,abs(BzBg(nLM+1:end)),'k-','linewidth',1, 'HandleVisibility','off'); 
loglog(S.LM_times,abs(Bz(1:nLM)),'r-','linewidth',1); 
loglog(S.HM_times,abs(Bz(nLM+1:end)),'r-','linewidth',1, 'HandleVisibility','off'); 
errorbar([S.LM_times; S.HM_times], [S.d_LM; S.d_HM], [S.sd_LM; S.sd_HM], '.'); 
grid on
ylabel('dBz/dt (V/Am^4)')
xlabel('Time (s)')
set(gca,'fontsize',14)
legend('Blocky model', 'Background', 'Realistic model','Synthetic Noisy Data')
set(gcf, 'Units','pixels', 'Position',[0 0 1200 600])
%--------------------------------------------------------------------------
%% Create input data file needed by RJMCMC code:
%

if S.saveEvery > S.numIterations
    S.saveEvery = S.numIterations;
end

if writeData
    save(DataFile, '-struct', 'S')
end    


 
 
 