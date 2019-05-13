function [] = GlacierGradient(Datafile)

%This function takes as input the path to a .mat file containing a
%structure with the following fields:
%
%XY:         an array of the (x,y) coordinates of the soundings that have been inverted (in UTM)
%
%lineNum:    an array of the line numbers (the inputs to the function SkyTEMDataProc)
%            of the inverted soundings
%
%CTP:        likely log10 resistivity (from median of conductance and average of thinnest and
%            thickest layers), upper estimate of log10 resistivity (from 25th percentile of
%            conductance and thinnest layer estimate), lower estimate of
%            log10 resisitivity (from 75th percentile of conductance and
%            thickest layer estimate)
%
%Marg:       median, 25th and 75th percentiles of marginal PDF over likely
%            depth of the conductive layer, computed from only the portion of the 
%            marginal PDF that is consistent with conductive sediment ( <~ 150 Ohm-m )
%
%This function does a bunch of plotting with this data

load(Datafile);

xy = TG.xy - TG.xy(1,:);  %choose a new xy origin
d = sqrt( xy(:,1).^2 + xy(:,2).^2 );  %distance of each sounding from origin
[d,qq] = sort(d); %sort by distance from the origin

% figure(1)
% plot(d,MargMedian(qq),'*','Linewidth',2)
% hold on
% plot(d,MargP75(qq),'*','Linewidth',2)
% plot(d,MargP25(qq),'*','Linewidth',2)
% xlabel('Distance from mouth of glacier (m)')
% ylabel('log_{10} resistivity (Ohm-m)')
% title('Resistivity of conductor vs distance from glacier mouth (marginal PDF method)')
% 
% figure(2)
% plot(d,CTPlikely,'*','Linewidth',2)
% hold on
% plot(d,CTPupper,'*','Linewidth',2)
% plot(d,CTPlower,'*','Linewidth',2)
% xlabel('Distance from mouth of glacier (m)')
% ylabel('log_{10} resistivity (Ohm-m)')
% title('Resistivity of conductor vs distance from glacier mouth (Conductance method)')
% 
% ErrorUpper = MargP75(qq) - MargMedian(qq); %for shadedErrorBar
% ErrorLower = -MargP25(qq) + MargMedian(qq);
% 
% figure(3)
% shadedErrorBar(d,MargMedian(qq),[ErrorUpper' ; ErrorLower'])
% xlabel('Distance from mouth of glacier (m)')
% ylabel('log_{10} resistivity (Ohm-m)')
% title('Resistivity of conductor vs distance from glacier mouth (marginal PDF method)')

%% Surf plots of conductivity/conductance/PFR vs r(glacier mouth)

load('C.mat');

for j=1:length(C)
    h = histogram(C{j});
    binEnd(j) = h.BinEdges(end);
end
[~,l] = max(binEnd);
h = histogram(C{l});
edges = h.BinEdges;

for j=1:length(C)
    P(:,j) = histc(C{j},edges);
    P(:,j) = P(:,j)/( sum(P(:,j)) * (edges(2)-edges(1)) );
end

n = length(C);
a = [ d(1:n)' (max(d(1:n))+100) ];   %this is all a silly pcolor thing...
[X_,Y_] = meshgrid(a,edges);
q = zeros(length(P),1);
P = [ P q ];

figure
pcolor(X_,Y_,P)
%surf(X_,Y_,P)
xlabel('distance from glacier mouth (m)')
ylabel('conductance')
zlabel('probability density')
title('PDF of conductance vs distance from glacier mouth')








