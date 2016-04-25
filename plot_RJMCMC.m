function [pdf_matrixH,pdf_matrixV,intfcCount,meanModelH,meanModelV,medianModelH,medianModelV,kOut] = ...
    plot_rjmcmc_new_parallel(nProcs,samples,kTracker,burnin,thin,binDepthInt,bits,isotropic,normalize,logPDF,S)

%KWK debug:
S.rhMin = S.log10rho_min;
S.rvMin = S.log10rho_min;
S.rhMax = S.log10rho_max;
S.rvMax = S.log10rho_max;
S.isotropic = true;

truth.z = S.ReferenceModel.z;
truth.rhoh = S.ReferenceModel.rho;

%%
zFixed = S.z;
%keep only the thinned samples
samples = samples(burnin+1:thin:end);
kTracker = kTracker(burnin+1:thin:end);


nSamples=size(samples,1);
nBins=ceil((S.zMax-zFixed(end))/binDepthInt);%can't place any interfaces till first zMin,though

s = cell(nProcs,1);%each proc contains a certain no. of samples
k = s; intfVec = s; binCount = s;

%initialize cells that will grow in the second dimension, each proc 
%processes 1/nProc samples
histcellH = cell(nProcs,1); confidH=zeros(nBins,2);
histcellV = cell(nProcs,1); confidV=zeros(nBins,2);

histcellH(:) = {cell(nBins,1)};
histcellV(:) = {cell(nBins,1)};

samplesPerProc = fix(nSamples/nProcs);
for ii=1:nProcs
    if ii~=nProcs
        s{ii} = samples((ii-1)*samplesPerProc +1 : ii*samplesPerProc);
        k{ii} = kTracker((ii-1)*samplesPerProc +1 : ii*samplesPerProc);
    else
        s{ii} = samples((ii-1)*samplesPerProc +1 : nSamples);
        k{ii} = kTracker((ii-1)*samplesPerProc +1 : nSamples);
    end
    
    intfVec{ii}  = zeros (sum(k{ii}),1);
    intfVecW{ii} = zeros (size(intfVec{ii})); %store weights at this temp
    binCount{ii} = zeros(nBins);
    %declared as cells because the length *could* change
    histcellH{ii}(:) = {zeros(1,nSamples)};
    histcellV{ii}(:) = {zeros(1,nSamples)};
end    
%we've split up the samples into diffefent cells now

kOut = kTracker;

medianModelH = zeros(nBins,1);medianModelV = zeros(nBins,1);
%maxH         = zeros(nBins,1);maxV         = zeros(nBins,1);
meanModelH         = zeros(nBins,1);meanModelV         = zeros(nBins,1);

depth_int  = (S.zMax-zFixed(end))/nBins;

parfor procInd = 1:nProcs
%     fid = fopen(['proc_',num2str(procInd),'_status'],'w');
    for ii=(1:length(s{procInd}))
        x=s{procInd}{ii};
        %get the interfaces
        %find first 0 in intfVector, to start from
        startLoc = find (~intfVec{procInd},1,'first');
        intfVec{procInd}(startLoc:startLoc + k{procInd}(ii)-1) ...
                = x.z;
         %attach S.zMax and last resistivity if not part of model
         if x.z(end)<S.zMax
             x.z = [x.z,S.zMax];
         end    
        binIndex=1;
        for jj=1:length(x.z)
            while x.z(jj) >= zFixed(end)+ depth_int*binIndex
                    binCount{procInd}(binIndex) = binCount{procInd}(binIndex) +1;
                    c = binCount{procInd}(binIndex);
                    histcellH{procInd}{binIndex}(c) = x.rhoh(jj);
                    histcellV{procInd}{binIndex}(c) = x.rhov(jj);
                    binIndex = binIndex +1;
            end
            if binIndex <= nBins
                binCount{procInd}(binIndex) = binCount{procInd}(binIndex) +1;
                c = binCount{procInd}(binIndex);
                histcellH{procInd}{binIndex}(c) = x.rhoh(jj);
                histcellV{procInd}{binIndex}(c) = x.rhov(jj);
            end    
        end
        if mod(ii,1000) == 0 & procInd == 1
%             fprintf(fid,'done %d out of %d\n',ii,length((s{procInd})));
	       fprintf('done %d out of %d\n',ii,length((s{procInd})));
        end
    end

    %crop to correct size
    for iiBin=1:nBins
        histcellH{procInd}{iiBin} = histcellH{procInd}{iiBin}(1:binCount{procInd}(iiBin));
        histcellV{procInd}{iiBin} = histcellV{procInd}{iiBin}(1:binCount{procInd}(iiBin));
    end    
end%parfor

tempH = cell(nBins,1); tempV = tempH; 
tempInt = [];

%now we have all the depth bins with the samples of rhoh and rhov
%find their histograms!!
edges = [S.rhMin:(S.rhMax-S.rhMin)/(bits):S.rhMax];
rho_int = edges(2)-edges(1);
pdf_matrixH = zeros(nBins,length(edges));
pdf_matrixV = zeros(nBins,length(edges));

%concatenate each process' histogram counts in the right bin (horizontally,
%dim 2)
for iBin=1:nBins
    for procInd = 1:nProcs
        tempH{iBin} = cat(2,tempH{iBin},histcellH{procInd}{iBin});
        tempV{iBin} = cat(2,tempV{iBin},histcellV{procInd}{iBin});
        if iBin == 1 %we only need do this once for the intfVec
          
           tempInt = [tempInt;intfVec{procInd}];

        end    
        
    end
    %find the histograms
    pdf_matrixH(iBin,:) = histc(tempH{iBin},edges);
    pdf_matrixV(iBin,:) = histc(tempV{iBin},edges);
     
end

%get rid of histc bin counts at the last edge value
pdf_matrixH(:,end) =[]; 
pdf_matrixV(:,end) =[];

intfVec = tempInt;
clear tempInt 
    
    
for i=1:nBins
    %normalize to pdf    
    pdf_matrixH(i,:) = pdf_matrixH(i,:)/sum(pdf_matrixH(i,:));
    pdf_matrixV(i,:) = pdf_matrixV(i,:)/sum(pdf_matrixV(i,:));
    if strcmp('normalize',normalize)
        pdf_matrixH(i,:) = pdf_matrixH(i,:)/max(pdf_matrixH(i,:));
        pdf_matrixV(i,:) = pdf_matrixV(i,:)/max(pdf_matrixV(i,:));
    end
    %find 5 and 95% confidence intervals
    [~,confidH(i,1)]=ismember (0,cumsum(pdf_matrixH(i,:))/sum(pdf_matrixH(i,:))>= 0.05,'legacy');
    [~,confidH(i,2)]=ismember (0,cumsum(pdf_matrixH(i,:))/sum(pdf_matrixH(i,:))>= 0.95,'legacy');
    [~,confidV(i,1)]=ismember (0,cumsum(pdf_matrixV(i,:))/sum(pdf_matrixV(i,:))>= 0.05,'legacy');
    [~,confidV(i,2)]=ismember (0,cumsum(pdf_matrixV(i,:))/sum(pdf_matrixV(i,:))>= 0.95,'legacy');
%     [~,confidH(i,1)]=ismember (0,cumsum(pdf_matrixH(i,:))/sum(pdf_matrixH(i,:))>= 0.01,'legacy');
%     [~,confidH(i,2)]=ismember (0,cumsum(pdf_matrixH(i,:))/sum(pdf_matrixH(i,:))>= 0.99,'legacy');
%     [~,confidV(i,1)]=ismember (0,cumsum(pdf_matrixV(i,:))/sum(pdf_matrixV(i,:))>= 0.01,'legacy');
%     [~,confidV(i,2)]=ismember (0,cumsum(pdf_matrixV(i,:))/sum(pdf_matrixV(i,:))>= 0.99,'legacy');
%find median model
     medianModelH(i) = median(tempH{i});
     medianModelV(i) = median(tempV{i});
%      [~,maxH(i)]     = max(pdf_matrixH(i,:));
%      [~,maxV(i)]     = max(pdf_matrixV(i,:));
     meanModelH(i)   = mean(tempH{i});
     meanModelV(i)   = mean(tempV{i});
%     modeModelH(i)   = mode(histcellH{i});
%     modeModelV(i)   = mode(histcellV{i});
    
end
confidH=confidH+1;confidV=confidV+1;%an ismember thing
if (strcmp(isotropic,'isotropic'))
    pdf_matrixV = pdf_matrixH;
    confidV     = confidH; 
end

%
% Plot the gridded results:
% 

hFig = figure;
if strcmpi('isotropic',isotropic)
   hRhov = subplot(1,3,1);
   hInterfaces = subplot(1,3,2);
   hLayers = subplot(1,3,3);
else
   hRhov = subplot(1,4,1);
   hRhoh = subplot(1,4,2);
   hInterfaces = subplot(1,4,3);
   hLayers = subplot(1,4,4);
end


    axes(hRhov);
    if strcmpi(logPDF,'logpdf')
        dat = log10([pdf_matrixV,pdf_matrixV(:,end);pdf_matrixV(end,:),pdf_matrixV(end,end)]);
        cbStr = 'log10(PDF)';
    else
        dat = ([pdf_matrixV,pdf_matrixV(:,end);pdf_matrixV(end,:),pdf_matrixV(end,end)]);
        cbStr = 'PDF';
    end
    pcolor(edges,zFixed(end)+(0:nBins)*depth_int,dat)
    set (gca,'ydir','reverse','layer','top')
    shading flat
    hold on
      plot(edges(1)+confidV(:,1)*rho_int,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'w','linewidth',1)
      plot(edges(1)+confidV(:,2)*rho_int,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'w','linewidth',1)
    %  plot(edges(1)+rho_int/2+(maxV-1)*rho_int,S.zMin+depth_int/2+(0:nBins-1)*depth_int,'m')
    %  plot(medianModelV,S.zMin+depth_int/2+(0:nBins-1)*depth_int,'.-r')
    % %  plot(modeModelV,S.zMin+depth_int/2+(0:nBins-1)*depth_int,'.-r')
    %   plot(meanModelV,  S.zMin+depth_int/2+(0:nBins-1)*depth_int,'*-y')
    if ~strcmpi(isotropic,'isotropic')
        title ('Vertical Resistivity')
    else
        title ('Resistivity')
    end
    xlabel ('Log_{10} (ohm-m)')
    %caxis ((2*[-sqrt(var(pdf_matrixV(:))),sqrt(var(pdf_matrixV(:)))]+mean(pdf_matrixV(:))));
    h=caxis;
    depthlim=ylim;
    hold on
     ylabel ('Depth (m)','fontsize',11);
        caxis([-3 0])
     
if ~strcmpi(isotropic,'isotropic')
    axes(hRhoh);
    pcolor(edges,zFixed(end)+(0:nBins)*depth_int,[pdf_matrixH,pdf_matrixH(:,end);pdf_matrixH(end,:),pdf_matrixH(end,end)])
    shading flat
    set (gca,'yticklabel',[])
    set (gca,'ydir','reverse','layer','top')
    title ('Horizontal resistivity')
    xlabel ('Log_{10} (ohm-m)')
    ylabel ('Depth (m)','fontsize',11);
    hold on
      plot(edges(1)+confidH(:,1)*rho_int,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'w','linewidth',1)
      plot(edges(1)+confidH(:,2)*rho_int,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'w','linewidth',1)
    %  plot(edges(1)+rho_int/2+(maxH-1)*rho_int,S.zMin+depth_int/2+(0:nBins-1)*depth_int,'m')
    %  plot(medianModelH,S.zMin+depth_int/2+(0:nBins-1)*depth_int,'.-r')
    % %  plot(modeModelH,S.zMin+depth_int/2+(0:nBins-1)*depth_int,'.-r')
    %   plot(meanModelH,  S.zMin+depth_int/2+(0:nBins-1)*depth_int,'*-y')

    hold on
 
    
end
%caxis(h);
h1=colorbar;
set(h1, 'Position', [.05 .11 .01 .8150])
set(get(h1,'ylabel'),'String', cbStr,'fontsize',11,'interpreter','tex');
 

%freezeColors
%histogram of interfaces
nBins=ceil((S.zMax-S.zMin)/binDepthInt);%this time with S.zMin as there are no interfaces prior to this
depth_int  = (S.zMax-S.zMin)/nBins;

axes(hInterfaces);
kBins=[S.zMin+(0:nBins)*depth_int];
[intfcCount]=histc(intfVec,kBins);
intfcCount(end)=[];
intfcCount=intfcCount/sum(intfcCount)/(kBins(2)-kBins(1));
colormap ('jet');
area([kBins(1):depth_int:kBins(end-1)]+0.5*(depth_int),intfcCount);
view (90,90)
xlim(depthlim);
ylabel ({'Interface probability'},'fontsize',11);
%yl=[0 4*(S.zMax-S.zMin)^-1]; ylim(yl);
%set(gca,'Ytick',[0:max(yl)/2:max(yl)],'YTickLabel',sprintf('%0.5f|',[0:max(yl)/2:max(yl)]))
hold all
plot (kBins,(S.zMax-S.zMin)^-1*ones(length(kBins),1),'--k')
set(gca, 'fontsize',11)
title ('Interface depth')
 

%plot truth
%plot rhov

if ~isempty(truth)
    
    if ~strcmpi(isotropic,'isotropic')
        subplot (hRhoh);
        plot_Truth(S,truth,2,'m--');
        hold on
    end

    %plot rhoh
    subplot (hRhov);
    plot_Truth(S,truth,2,'m--') ;
    hold on

end

%set(gcf, 'Units','inches', 'Position',[0 0 6.6 2.8])
%colormap(flipud(colormap ('gray')))
colormap('parula')

%plot hist of number of layers
axes(hLayers)
[a,b]=hist(kOut,S.kMin:S.kMax);
bar(b,a/sum(a)/(b(2)-b(1)))
hold on 
plot ([0 max(b)],(S.kMax-S.kMin+1)^-1*ones(1,2),'--k')
xlabel ('Number of interfaces','fontsize',11)
xlim([0,S.kMax])
%ylim([0,10*(S.kMax-S.kMin+1)^-1])
set(gca, 'fontsize',11)
ylabel ('Probability of interfaces','fontsize',11)
%set(gcf, 'Units','inches', 'Position',[0 0 3.3 2])
end

% function plot_model(S,x,which,lw)
%     if nargin<4
%         lw = 1;
%     end
%     if ~S.isotropic && nargin<3
%         beep;
%         disp ('specify H or V')
%     elseif S.isotropic
%         which = 'H';
%     end     
%     numInt = length(x.z);
%     earthmodel  = [S.z,x.z,S.rho(1,:),10.^x.rhoh,S.rho(2,:),10.^x.rhov];
%     S.numlayers = length(S.z)+numInt;
%     z=earthmodel(1:S.numlayers);
%     rho = [earthmodel(S.numlayers+1:2*S.numlayers);earthmodel(2*S.numlayers+1:3*S.numlayers)];
%     plotz=[];
%     plotz(1)=z(2);
%     plotrhoh=[];plotrhov=[];
%     for k=2:length(z)-1
%         plotz(2*k-2)=z(k+1);plotz(2*k-1)=z(k+1);
%         plotrhoh(2*k-3)=rho(1,k);plotrhoh(2*k-2)=rho(1,k);
%         plotrhov(2*k-3)=rho(2,k);plotrhov(2*k-2)=rho(2,k);
%     end
%     plotrhoh(2*k-1)=rho(1,end);plotrhov(2*k-1)=rho(2,end);hold all
%     if (which=='V')
%         plot(log10([plotrhov';plotrhov(end)]),([plotz,S.zMax]),'-k','linewidth',lw)
%     else
%         plot(log10([plotrhoh';plotrhoh(end)]),([plotz,S.zMax]),'-k','linewidth',lw)
%     end
% end 


function h = plot_Truth(S,x,lw,lc)
    if nargin<4
        lw = 1;
    end
 
    numInt = length(x.z);
     
    z   = x.z;
    rho = x.rhoh; 
    %[earthmodel(S.numlayers+1:2*S.numlayers)]; %;earthmodel(2*S.numlayers+1:3*S.numlayers)];
    plotz=[];
    plotz(1)=z(2);
    plotrhoh=[];plotrhov=[];
    for k=2:length(z)-1
        plotz(2*k-2)=z(k+1);plotz(2*k-1)=z(k+1);
        plotrhoh(2*k-3)=rho(1,k);plotrhoh(2*k-2)=rho(1,k);
        %plotrhov(2*k-3)=rho(2,k);plotrhov(2*k-2)=rho(2,k);
    end
    plotrhoh(2*k-1)=rho(1,end);%plotrhov(2*k-1)=rho(2,end);hold all
 
       h =  plot(log10([plotrhoh';plotrhoh(end)]),([plotz,S.zMax]),lc,'linewidth',lw);
%     end
end    

