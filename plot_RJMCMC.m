function [pdf_matrixH,pdf_matrixV,intfcCount,meanModelH,meanModelV,medianModelH,medianModelV,kOut] = ...
    plot_rjmcmc_new_parallel(nProcs,samples,kTracker,burnin,thin,binDepthInt,bits,isotropic,normalize,logPDF,S,z0)

% Trim to what has actually been computed so far:
iComputed = find(kTracker>0,1,'last');
kTracker = kTracker(1:iComputed);
samples = samples(1:iComputed);

%%
if ~exist('z0')
    z0 = 0.0;
end    
zFixed = S.z-z0;
%keep only the thinned samples
samples = samples(burnin+1:thin:end);
kTracker = kTracker(burnin+1:thin:end);

nSamples=size(samples,1);
nBins=ceil((S.zMax-zFixed(end))/binDepthInt);%can't place any interfaces till first zMin,though

s = cell(nProcs,1);%each proc contains a certain no. of samples
k = s; intfVec = s; binCount = s;

%initialize cells that will grow in the second dimension, each proc 
%processes 1/nProc samples
histcellH = cell(nProcs,1); 

histcellH(:) = {cell(nBins,1)};

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
    binCount{ii} = zeros(nBins);
    %declared as cells because the length *could* change
    histcellH{ii}(:) = {zeros(1,nSamples)};
end    
%we've split up the samples into diffefent cells now

kOut = kTracker;
depth_int  = (S.zMax-zFixed(end))/nBins;

for procInd = 1:nProcs
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
                    binIndex = binIndex +1;
            end
            if binIndex <= nBins
                binCount{procInd}(binIndex) = binCount{procInd}(binIndex) +1;
                c = binCount{procInd}(binIndex);
                histcellH{procInd}{binIndex}(c) = x.rhoh(jj);
            end    
        end
        if mod(ii,1000) == 0 & procInd == 1
	       fprintf('done %d out of %d\n',ii,length((s{procInd})));
        end
    end

    %crop to correct size
    for iiBin=1:nBins
        histcellH{procInd}{iiBin} = histcellH{procInd}{iiBin}(1:binCount{procInd}(iiBin));
    end    
end%parfor

tempH = cell(nBins,1);
tempInt = [];

%now we have all the depth bins with the samples of rhoh and rhov
%find their histograms!!
edges = [S.rhMin:(S.rhMax-S.rhMin)/(bits):S.rhMax];
pdf_matrixH = zeros(nBins,length(edges));

%concatenate each process' histogram counts in the right bin (horizontally,
%dim 2)
for iBin=1:nBins
    for procInd = 1:nProcs
        tempH{iBin} = cat(2,tempH{iBin},histcellH{procInd}{iBin});
        if iBin == 1 %we only need do this once for the intfVec
           tempInt = [tempInt;intfVec{procInd}];
        end    
    end
    %find the histograms
    pdf_matrixH(iBin,:) = histc(tempH{iBin},edges);
end

%get rid of histc bin counts at the last edge value
pdf_matrixH(:,end) =[]; 

intfVec = tempInt;
clear tempInt 

confidH_New = zeros(nBins,2);
for i=1:nBins
    %normalize to pdf    
    pdf_matrixH(i,:) = pdf_matrixH(i,:)/sum(pdf_matrixH(i,:));
    if strcmp('normalize',normalize)
        pdf_matrixH(i,:) = pdf_matrixH(i,:)/max(pdf_matrixH(i,:));
    end
    
    cdfL = cumsum(pdf_matrixH(i,:))/sum(pdf_matrixH(i,:));
    ind = find(cdfL<0.05,1,'last');
    if isempty(ind) || ind == 1
        confidH_New(i,1) = edges(1);
    else
        confidH_New(i,1) = edges(ind); %interp1(cdf(ind-1:ind),edges(ind)rho_int/2,0.05);
    end
    cdfR = cumsum(fliplr(pdf_matrixH(i,:)))/sum(pdf_matrixH(i,:));
    ind = find(cdfR<0.05,1,'last');
    edgesR = fliplr(edges);
    if  isempty(ind) || ind == 1
        confidH_New(i,2) = edgesR(1);
    else
        confidH_New(i,2) =  edgesR(ind); %interp1(cdf(ind-1:ind),edges(ind-1:ind)+rho_int/2,0.95);
    end
      
end

%
% Plot the gridded results:
% 

hFig = figure;

hRhov = subplot(1,3,1);
hInterfaces = subplot(1,3,2);
hLayers = subplot(1,3,3);


axes(hRhov);
if strcmpi(logPDF,'logpdf')
    dat = log10([pdf_matrixH,pdf_matrixH(:,end);pdf_matrixH(end,:),pdf_matrixH(end,end)]);
    cbStr = 'log10(PDF)';
else
    dat = ([pdf_matrixH,pdf_matrixH(:,end);pdf_matrixH(end,:),pdf_matrixH(end,end)]);
    cbStr = 'PDF';
end

%-1 for sigma from rho
pcolor(-edges,zFixed(end)+(0:nBins)*depth_int,dat)
set (gca,'ydir','reverse','layer','top')
shading flat
hold on

title ('Conductivity')

xlabel ('Log10 S/m')
depthlim=ylim;
hold on
ylabel ('Depth (m)');
set(gca,'fontsize',14)


plot(-confidH_New(:,1),zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'r','linewidth',2)
plot(-confidH_New(:,2),zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'r','linewidth',2)


h1=colorbar;
set(h1, 'Position', [.05 .11 .01 .8150])
set(get(h1,'ylabel'),'String', cbStr,'fontsize',14,'interpreter','tex');

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
xlabel('Depth (m)')
ylabel ('Interface probability');
%yl=[0 4*(S.zMax-S.zMin)^-1]; ylim(yl);
%set(gca,'Ytick',[0:max(yl)/2:max(yl)],'YTickLabel',sprintf('%0.5f|',[0:max(yl)/2:max(yl)]))
hold all
plot (kBins,(S.zMax-S.zMin)^-1*ones(length(kBins),1),'--k')
set(gca, 'fontsize',14)
title ('Interface depth')

colormap('parula')

%plot hist of number of layers
axes(hLayers)
[a,b]=hist(kOut,S.kMin:S.kMax);
bar(b,a/sum(a)/(b(2)-b(1)))
hold on 
plot ([0 max(b)],(S.kMax-S.kMin+1)^-1*ones(1,2),'--k')
xlabel ('Number of interfaces','fontsize',11)
xlim([0,S.kMax])
ylabel ('Probability of interfaces','fontsize',14)
set(gca, 'fontsize',14)
set(gcf, 'Units','pixels', 'Position',[0 0 1700 600])

%legacy 
pdf_matrixV = 0;
meanModelH = 0; meanModelV=0;medianModelH=0;medianModelV=0;
end


