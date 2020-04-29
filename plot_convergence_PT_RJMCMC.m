function plot_convergence_PT_RJMCMC(filePrefix)

nFiles = length(dir([filePrefix,'*_PT_RJMCMC*mat'])) - 1;%we don't want the swaps file

figure
for iFile = 1:nFiles
    load ([filePrefix,'_PT_RJMCMC_',num2str(iFile)])
    numAcceptWindows = length(AR_ll);
    ARwindowLength = size(en_ll,1)/numAcceptWindows;
    uAR = zeros(numAcceptWindows,1);
    bAR = zeros(numAcceptWindows,1);
    dAR = zeros(numAcceptWindows,1);
    mAR = zeros(numAcceptWindows,1);
    swaps = zeros(numAcceptWindows,1);
    
    firstEmpty = find(k_ll==0,1);
    if (isempty(firstEmpty))
        firstEmpty = length(k_ll)+1;
    end 
    s_ll = s_ll(1:firstEmpty-1);
    en_ll = en_ll(1:firstEmpty-1,:);
    k_ll = k_ll(1:firstEmpty-1);
    s1 = subplot(4,2,1);
    plotcol = 2;    
    if S_ll.MLerroradjust
        plotcol=1;
    end    
    plot(en_ll(:,plotcol))
    if iFile == nFiles
        plot(en_ll(:,plotcol),'-k','linewidth',2)
    end
    if S_ll.MLerroradjust
        title('-Log Likelihood')
    else
        title('RMS misfit per chain')
    end    
    hold all
    s2 = subplot(4,2,2);
    plot(k_ll)
    hold all
    title('interfaces per chain')
    if iFile == nFiles
        plot(k_ll,'-k','linewidth',2)
    end
    s3 = subplot(4,2,3);
    for iWindow = 1:numAcceptWindows
        uAR(iWindow) = AR_ll{iWindow}.uAR;
        bAR(iWindow) = AR_ll{iWindow}.bAR;
        dAR(iWindow) = AR_ll{iWindow}.dAR;
        mAR(iWindow) = AR_ll{iWindow}.mAR;
        swaps(iWindow) = AR_ll{iWindow}.swapRate;
    end    
    plot((1:numAcceptWindows)*ARwindowLength, uAR)
    title('MCMC update rate per chain')
    hold all
    if iFile == nFiles
        plot((1:numAcceptWindows)*ARwindowLength, uAR,'-k','linewidth',2)
    end
    s4 = subplot(4,2,4);
    plot((1:numAcceptWindows)*ARwindowLength, bAR)
    title('birth rate per chain')
    hold all
    if iFile == nFiles
        plot((1:numAcceptWindows)*ARwindowLength, bAR,'-k','linewidth',2)
    end
    s5 = subplot(4,2,5);
    plot((1:numAcceptWindows)*ARwindowLength, dAR)
    title('death rate per chain')
    hold all
    if iFile == nFiles
        plot((1:numAcceptWindows)*ARwindowLength, dAR,'-k','linewidth',2)
    end
    s6 = subplot(4,2,6);
    plot((1:numAcceptWindows)*ARwindowLength, mAR)
    title('iface move rate per chain')
    hold all
    if iFile == nFiles
        plot((1:numAcceptWindows)*ARwindowLength, mAR,'-k','linewidth',2)
    end
    s7 = subplot(4,2,7);
    plot((1:numAcceptWindows)*ARwindowLength,swaps)
    title('Swap rate per chain')
    hold all
    if iFile == nFiles
        plot((1:numAcceptWindows)*ARwindowLength, swaps,'-k','linewidth',2)
        linkaxes([s1 s2 s3 s4 s5 s6 s7],'x')
        set (gcf,'position',[0 0  1200 1200])
        xlim([0 firstEmpty])
    end
end