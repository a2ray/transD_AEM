% convergence stats
plot_convergence_PT_RJMCMC('Test_SynthData/Test_synthData')
subplot(421); ylim([0, 3])
load('Test_SynthData/Test_SynthData_PT_RJMCMC_6.mat')

% posterior ensemble
nBurnin = 50000; nThin = 1;
[H,V,intfcCount,meanH,meanV,medianH,medianV,kOut] = ...
    plot_RJMCMC(1,s_ll,k_ll ,nBurnin,nThin,1,100,'isotropic','normalize','NologPDF',S_ll);

% data fit for last model
Bz = get_field(S_ll, S_ll.ReferenceModel);
BzLast = get_field(S_ll, s_ll{end});

figure()
loglog(S_ll.times,abs(S_ll.data),'ro');
hold on
loglog(S_ll.times,abs(Bz),'k-','linewidth',1);
loglog(S_ll.times,abs(BzLast),'b-','linewidth',1);
grid()
ylabel('dBz/dt (V/Am^4)')
xlabel('Time (s)')
set(gca,'fontsize',14)
legend('Data', 'True model', 'Last model')
title('Data fit')
