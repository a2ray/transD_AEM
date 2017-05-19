% This script retrieves data from the _dat.xyz SkyTEM data file specified below under
% fileName, and loads all other setup info needed to run PT_RJMCMC

clear all
clc

fileName = 'TaylorGlacierSkyTEMdata/TaylorGlacier';
%fileName = 'TaylorGlacierSkyTEMdata/TaylorValley_dat.xyz';
lineNum = 2590;
Q = SkyTEMDataProc(fileName,lineNum);
%load('SynthDataNoise.txt');
%FID number
FID = Q.FID;
%gate times
times = Q.times;
%separate out the high mode and low mode times into separate arrays, get
%rid of NaNs from data and sd arrays, convert from rel error to abs error
k = 1;
l = 1;
for j=1:length(times)
    if( length(Q.dataHM) > 0 && isnan(Q.dataHM(j)) == 0 )
        HighMode.times(k) = times(j);
        HighMode.data(k) = Q.dataHM(j);
        HighMode.sd(k) = Q.dataHM(j)*Q.dataErrHM(j);
        k = k + 1;
    end
    if( length(Q.dataLM) > 0 && isnan(Q.dataLM(j)) == 0 )
        LowMode.times(l) = times(j);
        LowMode.data(l) = Q.dataLM(j);
        LowMode.sd(l) = Q.dataLM(j)*Q.dataErrLM(j);
        l = l + 1;
    end
end
%keyboard
% %data (high and low modes)
%HighMode.data = SynthDataNoise(1:19)';   
%LowMode.data = SynthDataNoise(20:end)'; 
%keyboard
% %data error (in correct units - not relative error)
% HighMode.sd = Q.dataHM.*Q.dataErrHM;
% LowMode.sd = Q.dataLM.*Q.dataErrLM;
%Height of transmitter and receivers (low and high modes) and
%xy-coordinates of the receiver (same datum for z as the model)
zTx     = -Q.alt; 
xyRx    = [17 0]; 
LowMode.zRx = (zTx - 1.90);
HighMode.zRx = (zTx - 0.19);         
% coordinates of the vertices of the transmitter loop polygon
xyPolyTx = [ -15.09  2.00 ;
             -15.09 -2.00 ; 
             -8.11 -10.16 ;
              8.11 -10.16 ; 
             15.09  -2.00 ;
             15.09   2.00 ;
              8.11  10.16 ; 
             -8.11  10.16 ];

%normalize the data by the area of the loop
area = polyarea(xyPolyTx(:,1),xyPolyTx(:,2));
if( length(Q.dataHM) > 0 )
    HighMode.data = HighMode.data/area;
    HighMode.sd = HighMode.sd/area;
end
if( length(Q.dataLM) > 0 )
    LowMode.data = LowMode.data/area;
    LowMode.sd = LowMode.sd/area;
end
                 
%Parameters needed for the Bayesian inversion
numIterations = 5e5;
saveEvery = 1e4;   
log10rho_min = -1;  %min resistivity allowed
log10rho_max = 5;    %max resistivity allowed
zMin = 0.0;          %min interface depth allowed
zMax = 600;          %max interface depth allowed
kMax = 35;           %max number of layers allowed
kMin = 1;            %min number of layers allowed
nFreqsPerDecade = 10;   %density of Fourier domain sampling
HankelFilterName = 'kk101Hankel.txt';   %Hankel filters
CosSinFilterName = 'kk101CosSin.txt';   %Cosine-Sine filters
LoopQuadOrder = 2;      %Gauss quadrature order for loop integeration
lowPassFilters = [ 4.500E+05  3.000E+05 ];  %Butterworth low-pass filters

%low-mode transmitter ramp: [ time    fraction-of-full-power ]
LowMode.ramp = [
-1.000E-03   0.00 
-9.800E-04   0.46 
-9.550E-04   0.79 
-9.200E-04   0.94 
-9.100E-04   0.95 
-8.000E-04   1.00 
0.000E+00    1.00 
4.400E-07    0.98 
8.000E-07    0.93 
1.100E-06    0.87 
1.580E-06    0.69 
2.160E-06    0.46 
2.940E-06    0.22 
3.900E-06    0.07 
5.340E-06    0.01 
6.480E-06    0.00 
];

%high-mode transmitter ramp
HighMode.ramp = [ 
 -1.000E-02   0.00 
-8.700E-03   0.36 
-7.100E-03   0.64 
-5.500E-03   0.81 
-3.300E-03   0.93 
-1.200E-03   0.98 
-5.000E-05   1.00 
 0.000E+00   1.00 
 3.700E-06   0.94 
 1.540E-05   0.72 
 3.080E-05   0.43 
 4.220E-05   0.21 
 5.090E-05   0.01 
 5.500E-05   0.00 
];

%number of layers in the model that are not inverted for (only air here)
nFixedLayers = 1; 
%none of this synthetic model matters but the first two: top of air and the ground (z = 0)
z   = [-1d5  0      50    100  150  ];   
%only the first conductivity (that of air) matters here
sig = [1d-12  1/1000 1/100 1/10 1/1000 ];   
% Fixed model layers:
rho = 1/sig(1:nFixedLayers);   
z   = z(1:nFixedLayers+1);    % Layer boundaries for fixed layers

%name of folder to save output from Bayesian inversion
outputFolder = 'SkyTEMinversions/Trash';
%name of file to save the structure that contains all of these parameters,
%data, etc
DataFile = ['SkyTEM-' num2str(lineNum) '-FID' num2str(FID) '.mat']; 
%save this setup to the .mat file
save(DataFile,'HighMode','LowMode','xyPolyTx','zTx','xyRx','LoopQuadOrder',...
              'HankelFilterName','CosSinFilterName','nFreqsPerDecade',...
              'kMin','kMax','zMin','zMax','log10rho_min','log10rho_max',...
              'numIterations','saveEvery','lowPassFilters','z','rho','FID','fileName');
          
%% Data plotting - modeled data validation

% U = load(['Trash/SkyTEM-' num2str(lineNum) '-FID' num2str(FID) '_PT_RJMCMC_4.mat']);
% iComputed = find(U.k_ll,1,'last');
% NtoPlot = 50;
% burnIn = 5000;  %if you don't want to use this, set it to 1
% BzHigh = zeros(NtoPlot,length(HighMode.times));
% for k=1:500
%     rand;
% end
% for l=1:NtoPlot
%     indexes(l) = ceil(iComputed*rand);
%     while( indexes(l) < burnIn )
%         indexes(l) = ceil(iComputed*rand);
%     end
%     ModelSig = 1./[rho 10.^(U.s_ll{indexes(l)}.rhoh)];
%     ModelZ = [z U.s_ll{indexes(l)}.z];
%     mu  = ones(size(ModelSig));
%     BzHigh(l,:) = get_LoopFields_TD_FHT(HighMode.times,xyPolyTx,zTx,xyRx,HighMode.zRx,...
%         ModelSig,mu,ModelZ,HankelFilterName,CosSinFilterName,nFreqsPerDecade,...
%         LoopQuadOrder,HighMode.ramp,lowPassFilters);
% end

% % ModelSig = 1./[rho 10.^(U.s_ll{modnum}.rhoh)];
% % ModelZ = [z U.s_ll{modnum}.z];
% % mu  = ones(size(ModelSig));
% % BzHigh = get_LoopFields_TD_FHT(HighMode.times,xyPolyTx,zTx,xyRx,HighMode.zRx,...
% %     ModelSig,mu,ModelZ,HankelFilterName,CosSinFilterName,nFreqsPerDecade,...
% %     LoopQuadOrder,HighMode.ramp,lowPassFilters);
% % BzHigh = BzHigh';

% eHM = HighMode.data+HighMode.sd;
% eHM = [ eHM ; HighMode.data-HighMode.sd];

% %plotting
% figure(7)
% loglog(HighMode.times,HighMode.data,'-ro','LineWidth',2)
% hold on
% loglog([HighMode.times;HighMode.times],eHM,'-k','LineWidth',2)
% axis([ 0.8*min(HighMode.times) 1.2*max(HighMode.times) 0.8*min(HighMode.data) 1.2*max(HighMode.data) ])
% for l=1:NtoPlot
%     loglog(HighMode.times,abs(BzHigh(l,:)),'-*','Color',[0.8 0.8 0.8])
% end
% xlabel('time (s)')
% ylabel('dBz/dt')
% title(['Datafit for ' num2str(NtoPlot) ' models'])
% hold off
% 
% v = ceil(NtoPlot*rand);
% r = (abs(HighMode.data)-abs(BzHigh(v,:)))./(HighMode.sd);
% figure(8)
% semilogx(HighMode.times,r,'o','LineWidth',2)
% axis([ 0.8*min(HighMode.times), 1.2*max(HighMode.times), 0.8*min(r), 1.2*max(r)])
% xlabel('time (s)')
% ylabel('Error-normalized residual')
% title(['Normalized residuals for randomly selected model'])

% %retrieve the smooth inversion result
% RegSol = SkyTEMInvProc(fileName,FID);
% RegSol.z = [ z(1) RegSol.z ];
% RegSol.sig = 1./[ rho RegSol.inv ];
% RegSol.mu  = ones(size(RegSol.sig));
% BzReg = get_LoopFields_TD_FHT(HighMode.times,xyPolyTx,zTx,xyRx,HighMode.zRx,...
%     RegSol.sig,RegSol.mu,RegSol.z,HankelFilterName,CosSinFilterName,nFreqsPerDecade,...
%     LoopQuadOrder,HighMode.ramp,lowPassFilters);
% BzReg = BzReg';
% 
% %plotting
% figure(9)
% loglog(HighMode.times,HighMode.data,'-ro')
% hold on
% loglog([HighMode.times;HighMode.times],eHM,'-','LineWidth',2)
% axis([ 0.8*min(HighMode.times) 1.2*max(HighMode.times) 0.8*min(HighMode.data) 1.2*max(HighMode.data) ])
% loglog(HighMode.times,abs(BzReg),'-*')
% xlabel('time (s)')
% ylabel('dBz/dt')
% title(['Datafit for model number ' num2str(modnum)])
% hold off
% 
% figure(10)
% semilogx(HighMode.times,(abs(HighMode.data)-abs(BzReg))./(HighMode.sd),'o','LineWidth',2)
% xlabel('time (s)')
% ylabel('Error-normalized residual')
% title(['Normalized residuals for model number ' num2str(modnum)])


%%  Exporting coordinates of data sounding
          
%output the lon-lat of the datapoints I have inverted
%Taylor Valley
% InvertedXY = [ Q.XYcoord(18650,:) ; Q.XYcoord(18800,:) ;...
%     Q.XYcoord(280,:) ; Q.XYcoord(6680,:) ];
% [Lon, Lat] = UTM2LonLat(InvertedXY(:,1),InvertedXY(:,2),58,'S');
% Qout = [Lon Lat];
% dlmwrite('InvertedLonLatData.txt',Qout,'Delimiter',',','Precision',8)

% %Taylor Glacier
% InvertedXY = [ Q.XYcoord(2425,:) ; Q.XYcoord(5,:) ;...
%     Q.XYcoord(3470,:) ; Q.XYcoord(3000,:) ; Q.XYcoord(500,:) ];
% [Lon, Lat] = UTM2LonLat(InvertedXY(:,1),InvertedXY(:,2),58,'S');
% Qout = [Lon Lat];
% dlmwrite('InvertedLonLatData.txt',Qout,'Delimiter',',','Precision',8)


PropXY = [ Q.XYcoord(lineNum,:) ];
[Lon, Lat] = UTM2LonLat(PropXY(:,1),PropXY(:,2),58,'S');
Qout = [Lon Lat];
dlmwrite(['LonLat/LonLat' num2str(lineNum) '.txt'],Qout,'Delimiter',',','Precision',8)
dlmwrite(['xy/xy' num2str(lineNum) '.txt'],PropXY,'Delimiter',',','Precision',8)

%plots the convergence of these chains
% tmp = ['SkyTEMinversions/Trash/SkyTEM-' num2str(lineNum) '-FID' num2str(FID) '_PT_RJMCMC'];
% plot_convergence_PT_RJMCMC(tmp)

% %% Correlation and Covariance
% 
% U = load(['Trash/SkyTEM-' num2str(lineNum) '-FID' num2str(FID) '_PT_RJMCMC_4.mat']);
% iComputed = find(U.k_ll,1,'last');
% NtoPlot = 1000;
% burnIn = 5000;  %if you don't want to use this, set it to 1
% 
% BzHigh = zeros(NtoPlot,length(HighMode.times));
% r = BzHigh;
% 
% for k=1:500
%     rand;
% end
% 
% for l=1:NtoPlot
%     indexes(l) = ceil(iComputed*rand);
%     while( indexes(l) < burnIn )
%         indexes(l) = ceil(iComputed*rand);
%     end
%     ModelSig = 1./[rho 10.^(U.s_ll{indexes(l)}.rhoh)];
%     ModelZ = [z U.s_ll{indexes(l)}.z];
%     mu  = ones(size(ModelSig));
%     BzHigh(l,:) = get_LoopFields_TD_FHT(HighMode.times,xyPolyTx,zTx,xyRx,HighMode.zRx,...
%         ModelSig,mu,ModelZ,HankelFilterName,CosSinFilterName,nFreqsPerDecade,...
%         LoopQuadOrder,HighMode.ramp,lowPassFilters);
%     r(l,:) = (abs(HighMode.data)-abs(BzHigh(l,:)))./HighMode.sd;
%     if( mod(l,100) == 0 )
%         fprintf('%d\n',l)
%     end
% end
% sd = HighMode.sd;
% save('Residuals.mat','r','sd')

