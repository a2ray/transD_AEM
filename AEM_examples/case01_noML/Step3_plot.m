% convergence stats
plot_convergence_PT_RJMCMC('Test_SynthData/Test_synthData')
subplot(421); ylim([0, 3])
load('Test_SynthData/Test_SynthData_PT_RJMCMC_4.mat')

% posterior ensemble
nBurnin = 50000; nThin = 2;
[H,V,intfcCount,meanH,meanV,medianH,medianV,kOut] = ...
    plot_RJMCMC(4,s_ll,k_ll ,nBurnin,nThin,1,100,'isotropic','normalize','NologPDF',S_ll);

% data fit for last model
firstEmpty = find(k_ll==0,1);
if (isempty(firstEmpty))
    firstEmpty = length(k_ll)+1;
end 
s_ll = s_ll(1:firstEmpty-1);
Bz = get_field(S_ll, S_ll.ReferenceModel);
BzLast = get_field(S_ll, s_ll{end});

figure()
loglog(S.LM_times,abs(BzLast(1:nLM)),'k-','linewidth',1); hold on
loglog(S.HM_times,abs(BzLast(nLM+1:end)),'k-','linewidth',1, 'HandleVisibility','off'); 
loglog(S.LM_times,abs(Bz(1:nLM)),'r-','linewidth',1); 
loglog(S.HM_times,abs(Bz(nLM+1:end)),'r-','linewidth',1, 'HandleVisibility','off'); 
errorbar([S.LM_times; S.HM_times], [S.d_LM; S.d_HM], [S.sd_LM; S.sd_HM], '.'); 