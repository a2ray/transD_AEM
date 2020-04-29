%% Compare with Analytical formula central loop sounding on surface of a halfspace:
clear all

% For testing with halfpsace equation:

% % % formula from Christiansen et el groundwater chapter 
% % %eqn 6.26:


nFreqsPerDecade = 15;
LoopQuadOrder   = 3;   % Gauss quadrature order for polygon loop integration
loopRadius      = 15; % 10 m circular loop for Christiansen formula

nPolygonSides   = 500;  % number of polygon sides. more means better circle approximation. must be >2

HankelFilterName = 'kk201Hankel.txt';
CosSinFilterName = 'kk201CosSin.txt';

z   = [-1d5     0     ];   % m, vertical position of layer top boundaries, 1st one ignored since it's the top of air.
sig = [1d-20   1/1   ]; 


%-------------------------------
% Don't change anything below here:

tlong       = logspace(-6,1,40);
xyRx        = [0 0];  % has to be 0 since Christiansen formula is only for receiver coil in center of transmitter loop
zRx         = 0; % Rx depth. has to be 0 since Christiansen formula (from Ward and Hohmann actually) is valid only on surface of halspace
zTx         = 0; % Tx depth. has to be 0 since Christiansen formula (from Ward and Hohmann actually) is valid only on surface of halspace
mu          = ones(size(sig)); % relative magnetic permeability of each layer. Set this to 1 always!

% Make a polygon approximation of a circle to test get_PolygonFields versus
% Christiansen et al analytic formula for perfect circle
 
L = -pi/2+linspace(0,2*pi,nPolygonSides+1);
L = L(1:end-1);
xyPolyTx(:,1) = loopRadius*cos(L)';
xyPolyTx(:,2) = loopRadius*sin(L)';
rRx = norm(xyRx);
rTxLoop  = loopRadius;                           % from area = pi*r^2

% Nothing to change below here:
sig0 = sig(2);
theta = sqrt(4*pi*1d-7*sig0 ./ (4*tlong) );
a     = loopRadius;  
dbzdt = -1/(sig0*a^3)*(3*erf(theta*a) - 2/sqrt(pi)*theta.*a.*(3+2*theta.^2*a^2).*exp(-theta.^2*a^2));
dbzdt = dbzdt/(pi*a^2);


% Compute loop response for the polygon, without any waveform ramp or
% filters:
BzPoly = get_LoopFields_TD_FHT(tlong,xyPolyTx,zTx,xyRx,zRx,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade,LoopQuadOrder);    

BzCircle = get_LoopFields_circle_TD_FHT(tlong,rTxLoop,zTx,rRx,zRx,sig,mu,z,...
                         HankelFilterName,CosSinFilterName,nFreqsPerDecade);
figure;
subplot(2,1,1)
loglog(tlong,-(dbzdt),'k-','linewidth',1);
hold all;
loglog(tlong,-(BzPoly),'r--','linewidth',1);
loglog(tlong,-(BzCircle),'linewidth',1);
grid()
legend('Analytic','Polygon','Circle')
xlabel('Time (s)')
title(sprintf('LoopQuadOrder=%i, nPolygonSides=%i, nFreqsPerDecade=%i',LoopQuadOrder,nPolygonSides,nFreqsPerDecade))
set(gca,'fontsize',14)
 
subplot(2,1,2)
loglog(tlong,100*abs(dbzdt'-BzPoly)./abs(dbzdt)','-','linewidth',1);
hold all
loglog(tlong,100*abs(dbzdt'-BzCircle)./abs(dbzdt)','-','linewidth',1);
grid()
legend('Polygon','Circle')
ylabel('Relative Difference (%)')
xlabel('Time (s)')
set(gca,'fontsize',14)