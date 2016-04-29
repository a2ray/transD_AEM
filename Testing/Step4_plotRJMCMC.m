 
%
% Read in PT_RJMCMC results and plot them up:
%
nBurnIn = 500; % skip first nBurnIn samples. Set this to value determined from Figure made in Step 3.
nThin   = 1;   % thin chain by this factor. makes plotting faster, but dont' skip too many samples or the plot will look jagged. 


%---------------------------------------------------

S = load(DataFile);
[~,FileRoot] = fileparts(DataFile);

load(strcat(outputFolder,'/',FileRoot,'_PT_RJMCMC_12.mat'))
 
 

% Trim to what has actually been computed so far:
iComputed = find(k_ll,1,'last');  
k_ll = k_ll(1:iComputed);
s_ll = s_ll(1:iComputed);


[H,V,intfcCount,meanH,meanV,medianH,medianV,kOut] = ...
plot_RJMCMC(4,s_ll,k_ll ,nBurnIn,nThin,5,50,'isotropic','NOnormalize','logPDF',S);

% function [pdf_matrixH,pdf_matrixV,intfcCount,meanModelH,meanModelV,medianModelH,medianModelV,kOut] =
% plot_rjmcmc_new_parallel(nProcs,samples,kTracker,burnin,thin,binDepthInt,bits,isotropic,normalize,S)
