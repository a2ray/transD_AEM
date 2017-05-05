function [pdf_matrixH,pdf_matrixV,intfcCount,meanModelH,meanModelV,medianModelH,medianModelV,kOut] = ...
    plot_rjmcmc_new_parallel(nProcs,samples,kTracker,burnin,thin,binDepthInt,bits,isotropic,normalize,logPDF,S)

%KWK debug:
S.rhMin = S.log10rho_min;
S.rvMin = S.log10rho_min;
S.rhMax = S.log10rho_max;
S.rvMax = S.log10rho_max;
S.isotropic = true;

%truth.z = S.ReferenceModel.z;
%truth.rhoh = S.ReferenceModel.rho;

%%
zFixed = S.z;
%keep only the thinned samples
samples = samples(burnin+1:thin:end);
kTracker = kTracker(burnin+1:thin:end);


nSamples=size(samples,1);
nBins=ceil((S.zMax-zFixed(end))/binDepthInt);%can't place any interfaces till first zMin,though

%% porosity vs salinity pdfs

tic
m = 2; %Archie's Law exponent (sig = sig_w*phi^m)
nPhiSampling = 20; %how many values of porosity and pore fluid resistivity you want to evaluate
phi_ = linspace(-2,-0.69,nPhiSampling); %log10(porosity) values to calculate sig_w over
rhoW_ = linspace(-1.5,0.5,nPhiSampling); %log10(pore fluid resisitivity) values to calc porosity over
zBrine = 195; %depth (m) at which we want to evaluate brine resistivity
rhoBrine = zeros(length(samples),1);
Porosity_PFR = zeros(nSamples,nPhiSampling,2); %holds the porosity and pore fluid resistivity values
for j=1:nSamples
    %determine the resistivity of this model at the brine depth (zBrine)
    for l=1:length(samples{j}.z)
        if( samples{j}.z(l) > zBrine )
            rhoBrine(j) = samples{j}.rhov(l);
            break
        elseif( samples{j}.z(l) < zBrine && l == length(samples{j}.z) )
            rhoBrine(j) = samples{j}.rhov(l+1);
        end
    end
    %compute, for each pore fluid resistivity, and for this model, the log10 porosity
    Porosity_PFR(j,:,1) = (1/m) * ( rhoW_ - rhoBrine(j) );
    %compute, for each porosity, and for this model, the log10 pore fluid resistivity
    Porosity_PFR(j,:,2) = m*phi_ + rhoBrine(j);
end
toc
%now build the histograms for each of the porosity and pore fluid resistivity values
nbinsPoroPFR = 50;
a = min(min(Porosity_PFR(:,:,2)));
b = max(max(Porosity_PFR(:,:,2)));
binEdgesPFR = a:(b-a)/nbinsPoroPFR:b;
a = min(min(Porosity_PFR(:,:,1)));
b = max(max(Porosity_PFR(:,:,1)));
binEdgesPoro = a:(b-a)/nbinsPoroPFR:b;
PoroPDF = zeros(nPhiSampling,nbinsPoroPFR+1);
PFRPDF = zeros(nPhiSampling,nbinsPoroPFR+1);
%keyboard
for j=1:nPhiSampling
    PoroPDF(j,:) = histc(Porosity_PFR(:,j,1),binEdgesPoro);
    PFRPDF(j,:) = histc(Porosity_PFR(:,j,2),binEdgesPFR);
end
%create the meshgrids for plotting
[RHOW_,PORO_] = meshgrid(binEdgesPFR,phi_);
[PORO,RHOW] = meshgrid(binEdgesPoro,rhoW_);
figure(1)
pcolor(PORO,RHOW,PoroPDF)
xlabel('log_{10} porosity (%)')
ylabel('log_{10} pore fluid resistivity (ohm-m) (chosen values)')
title('Porosity likelihood for chosen values of pore fluid resistivity')
figure(2)
pcolor(RHOW_,PORO_,PFRPDF)
xlabel('log_{10} pore fluid resistivity (ohm-m) (predicted values)')
ylabel('log_{10} porosity (%) (chosen values)')
title('Pore fluid resistivity likelihood for chosen values of porosity')
%keyboard

%%

s = cell(nProcs,1);%each proc contains a certain no. of samples
k = s; intfVec = s; binCount = s;

%initialize cells that will grow in the second dimension, each proc 
%processes 1/nProc samples
histcellH = cell(nProcs,1); confidH=zeros(nBins,2);confidH_New=zeros(nBins,2);
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


% parameters for statistics of conductive layer
conductor_top = 290.0;
conductor_bottom = 360.0;
conductor_samples = [];


%% messing around with covariance plots

%specify a certain subset of models to look at (e.g. those with n interfaces)
specialN = 5; %look at models with exactly specialN interfaces
specFlag = false; %use this when a model has specialN interfaces
specHist = histcellH; %histogram for only these models
specBinCount = binCount; %bin counts for only these models
specZ = zeros(specialN,length(samples));
specRho = zeros(specialN+1,length(samples));
j = 1;
for l=1:length(samples)   
    %look at models with a specific number of interfaces
    if( length(samples{l}.z) == specialN )
        specZ(:,j) = samples{l}.z;
        specRho(:,j) = samples{l}.rhov;
        j = j+1;
    end
end

% specZ = specZ(:,1:j);
% specZthick = diff([ zeros(1,size(specZ,2)) ; specZ ]);
% specRho = specRho(:,1:j);
% Z = 0:10:600;
% Thick = 0:5:max(specZthick);
% RHO = 0:0.1:5;
% [RHO_,Z_] = meshgrid(RHO,Z);
% [RHOT_,Thick_] = meshgrid(RHO,Thick);
% pdf_zRho = 0*Z_;
% pdf_ThickRho = 0*Thick_;
% for l=1:length(RHO)-1
%     for m=1:length(Thick)-1
%         for j=1:size(specZthick,2)
%             if( specZthick(2,j) > Thick(m) && specZthick(2,j) <= Thick(m+1) && ...
%                     specRho(2,j) > RHO(l) && specRho(2,j) <= RHO(l+1) )
%                 pdf_ThickRho(m,l) = pdf_ThickRho(m,l) + 1;
%             end
%         end
%     end
%     for k=1:length(Z)-1
%         for j=1:size(specZ,2)
%             if( specZ(1,j) > Z(k) && specZ(1,j) <= Z(k+1) && ...
%                     specRho(2,j) > RHO(l) && specRho(2,j) <= RHO(l+1) )
%                 pdf_zRho(k,l) = pdf_zRho(k,l) + 1;
%             end
%         end
%     end
% end
% figure(11)
% pcolor(RHO_,Z_,pdf_zRho)
% h = gca;
% h.YDir = 'reverse';
% shading interp
% figure(12)
% pcolor(RHOT_,Thick_,pdf_ThickRho)
% shading interp

%keyboard

%%

medianModelH = zeros(nBins,1);medianModelV = zeros(nBins,1);
%maxH         = zeros(nBins,1);maxV         = zeros(nBins,1);
meanModelH         = zeros(nBins,1);meanModelV         = zeros(nBins,1);

depth_int  = (S.zMax-zFixed(end))/nBins;


parfor procInd = 1:nProcs
    fprintf('we are here\n')
%     fid = fopen(['proc_',num2str(procInd),'_status'],'w');
    specFlag = false;
    for ii=(1:length(s{procInd}))
        x=s{procInd}{ii};
        if( length(x.z) == 5)
            specFlag = true;
        end
        %get the interfaces
        %find first 0 in intfVector, to start from
        startLoc = find (~intfVec{procInd},1,'first');
        intfVec{procInd}(startLoc:startLoc + k{procInd}(ii)-1) ...
                = x.z;
         %attach S.zMax and last resistivity if not part of model
         if x.z(end)<S.zMax
             x.z = [x.z,S.zMax];
         end
        if( specFlag == true )
            binIndex=1;
            for jj=1:length(x.z)
                while x.z(jj) >= zFixed(end)+ depth_int*binIndex
                        binCount{procInd}(binIndex) = binCount{procInd}(binIndex) +1;
                        specBinCount{procInd}(binIndex) = specBinCount{procInd}(binIndex) +1;
                        c = binCount{procInd}(binIndex);
                        c_ = specBinCount{procInd}(binIndex);
                        histcellH{procInd}{binIndex}(c) = x.rhoh(jj);
                        specHist{procInd}{binIndex}(c_) = x.rhoh(jj);
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
        else
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
        end
        if mod(ii,1000) == 0 & procInd == 1
%             fprintf(fid,'done %d out of %d\n',ii,length((s{procInd})));
	       fprintf('done %d out of %d\n',ii,length((s{procInd})));
        end
    end

    %crop to correct size
    for iiBin=1:nBins
        histcellH{procInd}{iiBin} = histcellH{procInd}{iiBin}(1:binCount{procInd}(iiBin));
        specHist{procInd}{iiBin} = specHist{procInd}{iiBin}(1:specBinCount{procInd}(iiBin));
        histcellV{procInd}{iiBin} = histcellV{procInd}{iiBin}(1:binCount{procInd}(iiBin));
    end    
end%parfor

tempH = cell(nBins,1); tempV = tempH; tempSpec = tempH; 
tempInt = [];

%now we have all the depth bins with the samples of rhoh and rhov
%find their histograms!!
edges = [S.rhMin:(S.rhMax-S.rhMin)/(bits):S.rhMax];
rho_int = edges(2)-edges(1);
pdf_matrixH = zeros(nBins,length(edges));
pdf_matrixSpec = zeros(nBins,length(edges));
pdf_matrixV = zeros(nBins,length(edges));

%concatenate each process' histogram counts in the right bin (horizontally,
%dim 2)
for iBin=1:nBins
    for procInd = 1:nProcs
        tempH{iBin} = cat(2,tempH{iBin},histcellH{procInd}{iBin});
        tempSpec{iBin} = cat(2,tempSpec{iBin},specHist{procInd}{iBin});
        tempV{iBin} = cat(2,tempV{iBin},histcellV{procInd}{iBin});
        if iBin == 1 %we only need do this once for the intfVec
          
           tempInt = [tempInt;intfVec{procInd}];

        end    
        
    end
    %find the histograms
    pdf_matrixH(iBin,:) = histc(tempH{iBin},edges);
    pdf_matrixSpec(iBin,:) = histc(tempSpec{iBin},edges);
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
     %compute the mode
     pp = histogram(tempH{i},edges);
     [~,jj] = max(pp.Values);
     modeModelH(i)   = (pp.BinEdges(jj)+pp.BinEdges(jj+1))/2;
     %compute 25th and 75th percentiles
     confid75H(i)    = prctile(tempH{i},75);
     confid25H(i)    = prctile(tempH{i},25);
      
     %for generating the marginal PDF over the conductive region
     if( depth_int*i > conductor_top && depth_int*i < conductor_bottom )
         conductor_samples = [ conductor_samples tempH{i} ];
         %Resistance = Resistance + tempH{i}*depth_int;
     end
          
end

Resistance = zeros(1,length(samples));
Conductance = Resistance;
dz = 1;%depth_int;
%compute resistance = rho*thickness over the depth of the conductor
tic
for l=1:length(samples)
    Z = [ 0 samples{l}.z S.zMax ];
    z = conductor_top;
    for j=1:length(Z)-1
        if( z > conductor_bottom )
            break
        end
        while( z > Z(j) && z <= Z(j+1) )           
            Resistance(l) = Resistance(l) + 10.^(samples{l}.rhoh(j));
            Conductance(l) = Conductance(l) + 10.^(-samples{l}.rhoh(j));
            z = z + dz;
        end
    end   
end
Resistance = Resistance*dz;
Conductance = Conductance*dz;
toc

%plot the marginal PDF over the depth of the conductor
conductorMarginalPDF = histc(conductor_samples,edges);
conductorMarginalPDF = conductorMarginalPDF/((edges(3)-edges(2))*sum(conductorMarginalPDF));
figure(99)
plot(edges,conductorMarginalPDF,'LineWidth',2)
xlabel('log_{10} resistivity (Ohm-m)')
ylabel('Bin count')
title('Marginal PDF over the depth of the conductor')

%compute the median, 25th, and 75th percentiles over a "conductor" portion
%of marginal PDF
a = 0.32; %minimum rho consistent with conductor portion of marginal PDF
b = 2.24; %maximum rho consistent with conductor portion of marginal PDF
ConductiveSignal = conductor_samples(conductor_samples>=a & conductor_samples<=b);
fprintf('\nmedian conductivity: %f\n',median(ConductiveSignal))
fprintf('25th percentile of conductivity: %f,\n', prctile(ConductiveSignal,25))
fprintf('75th percentile of conductivity: %f\n\n', prctile(ConductiveSignal,75))

%conductance statistics and plotting
fprintf('\nMedian Conductance = %f\n',median(Conductance))
fprintf('25th percentile of Conductance = %f\n',prctile(Conductance,25))
fprintf('75th percentile of Conductance = %f\n\n',prctile(Conductance,75))
x1 = 210.0;
x2 = 237.5;
x3 = 262.5;
maxT = 60;%x3-x1;
minT = 20;%x2-x1;
likelyT = 42.5;%(maxT+minT)/2;
fprintf('\nAn upper estimate on log10 resistivity: %f\n',-log10(prctile(Conductance,25)/maxT))
fprintf('A lower estimate on log10 resistivity: %f\n',-log10(prctile(Conductance,75)/minT))
fprintf('A likely estimate on log10 resistivity: %f\n\n',-log10(median(Conductance)/likelyT))

figure(97)
tmp = Conductance( Conductance < prctile(Conductance,100));
histogram(tmp)
xlabel('Conductance  ( sigma x thickness )')
ylabel('Bin count')
title('Conductivity-thickness product histogram')


%retrieve the smooth inversion result
RegSol = SkyTEMInvProc(S.fileName,S.FID);
%0 to z_max, gridded finely
binnedZ = zFixed(end)+depth_int/2+(0:nBins-1)*depth_int;
RegSolInterp = spline(RegSol.z,log10(RegSol.inv),binnedZ);


confidH=confidH+1;confidV=confidV+1;%an ismember thing
if (strcmp(isotropic,'isotropic'))
    pdf_matrixV = pdf_matrixH;
    confidV     = confidH; 
end


%estimate DOI ( x in [0,1] )
DOI = (confidH_New(:,2) - confidH_New(:,1))/(S.rhMax-S.rhMin);
figure(101)
plot(DOI,binnedZ)



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
%       plot(edges(1)+confidV(:,1)*rho_int-rho_int/2,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'r','linewidth',2)
%       plot(edges(1)+confidV(:,2)*rho_int+rho_int/2,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'r','linewidth',2)
       plot(confidH_New(:,1),zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'r','linewidth',2)
       plot(confidH_New(:,2),zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'r','linewidth',2)
       %plot(confid25H,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'y','linewidth',2)
       %plot(confid75H,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'y','linewidth',2)
       %plot(meanModelH,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'b','linewidth',2)
       plot(medianModelH,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'k','linewidth',2)
       %plot(modeModelH,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'b','linewidth',2)
       plot(RegSolInterp,zFixed(end)+depth_int/2+(0:nBins-1)*depth_int,'g','linewidth',2)
       plot(linspace(-1,5,length(binnedZ)),RegSol.DOI(1)*ones(length(binnedZ)),'m','linewidth',2)
       plot(linspace(S.rhMin,S.rhMax,length(binnedZ)),RegSol.DOI(2)*ones(length(binnedZ)),'m','linewidth',2)


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

% if ~isempty(truth) && ~isempty(truth.z)
%     
%     if ~strcmpi(isotropic,'isotropic')
%         subplot (hRhoh);
%         plot_Truth(S,truth,2,'m--');
%         hold on
%     end
% 
%     %plot rhoh
%     subplot (hRhov);
%     plot_Truth(S,truth,2,'m--') ;
%     hold on
% 
% end

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


% function h = plot_Truth(S,x,lw,lc)
%     if nargin<4
%         lw = 1;
%     end
%  
%     numInt = length(x.z);
%      
%     z   = x.z;
%     rho = x.rhoh; 
%     %[earthmodel(S.numlayers+1:2*S.numlayers)]; %;earthmodel(2*S.numlayers+1:3*S.numlayers)];
%     plotz=[];
%     plotz(1)=z(2);
%     plotrhoh=[];plotrhov=[];
%     for k=2:length(z)-1
%         plotz(2*k-2)=z(k+1);plotz(2*k-1)=z(k+1);
%         plotrhoh(2*k-3)=rho(1,k);plotrhoh(2*k-2)=rho(1,k);
%         %plotrhov(2*k-3)=rho(2,k);plotrhov(2*k-2)=rho(2,k);
%     end
%     plotrhoh(2*k-1)=rho(1,end);%plotrhov(2*k-1)=rho(2,end);hold all
%  
%        h =  plot(log10([plotrhoh';plotrhoh(end)]),([plotz,S.zMax]),lc,'linewidth',lw);
% %     end
% end    
% 
