

% Step 3:
% Read in the _12.mat file and plot up the misfit to see how many samples
% need to be trimmed during the burn-in period at the start:
%
[~,FileRoot] = fileparts(DataFile);

load(strcat(outputFolder,'/',FileRoot,'_PT_RJMCMC_12.mat'))
 
% Trim to what has actually been computed so far:
iComputed = find(k_ll,1,'last'); 
en_ll =  en_ll(1:iComputed,:);

figure
plot(en_ll(:,2))
xlabel('iterations #')
ylabel('RMS misfit')

title('Use this to find burn-in number of samples');


% If you want to plot up some model realizations:
%
% nBurnIn = 1000;
% nEnd    = 10000;
% nThin = 10;

%
% iComputed = find(k_ll,1,'last'); 
% s_ll =  s_ll(1:iComputed,:);
% 
% iRange = nStart:nThin:nEnd;
% 
% 
% figure
%
% x.z = [];
% x.rhoh = [];
% Truth.zMax = zMax;
% hTrue = plot_model(Truth,x,3,'r');
% hold on;
% 
% for i = iRange;% 1:length(s_ll)-1
%     h = plot_model(S_ll,s_ll{i},.25,'k');
%   
% end
% hTrue = plot_model(Truth,x,3,'r');
%  hLast = plot_model(S_ll,s_ll{iRange(end)},3,'b');
% 
% set(gca,'ydir','reverse')
% xlabel('log10 (ohm-m)')
% ylabel('Depth (m)')
%  legend([hTrue;hLast],{'Truth','Last RJMCMC'})