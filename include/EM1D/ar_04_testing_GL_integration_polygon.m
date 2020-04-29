%% Test convergence of Polygon integration in the time domain:
% if satisifed with analytical comparisons to circular loop
% run this test for the order of integration

clear all;

xyPolyTx = [ -15.09 -2.00 ; 
             -8.11 -10.16 ;
              8.11 -10.16 ; 
             15.09  -2.00 ;
             15.09   2.00 ;
              8.11  10.16 ; 
             -8.11  10.16 ; 
             -15.09  2.00 ];
           
zTx     = -30; 
zRx     = -31.90;
xyRx    = [17 0]; 
 
nFreqsPerDecade = 10;

HankelFilterName = 'kk201Hankel.txt';
CosSinFilterName = 'kk201CosSin.txt';
lowPassFilters = [ 4.500E+05  3.000E+05 ]; 

z   = [-1d5     0     ];   % m, vertical position of layer top boundaries, 1st one ignored since it's the top of air.
sig = [1d-20   1/10   ]; 
mu  = ones(size(sig)); % relative magnetic permeability of each layer. Set this to 1 always!

tlong = logspace(-5,-1, 40); 
 
nOrderMax = 8;
for order = 1:nOrderMax
   BzPoly(order,:) = get_LoopFields_TD_FHT(tlong,xyPolyTx,zTx,xyRx,zRx,sig,mu,z, HankelFilterName,CosSinFilterName,nFreqsPerDecade, order);        
end
 
 
figure;
subplot(2,1,1)
loglog(tlong,abs(BzPoly(nOrderMax,:)),'-','linewidth',1);
title(sprintf('nFreqsPerDecade=%i, Rx range=%.0f',nFreqsPerDecade,xyRx(1)))
ylabel('dBz/dt (V/Am^4)')
xlabel('Time (s)')
set(gca,'fontsize',14)

subplot(2,1,2)
for order = 1:nOrderMax-1
    loglog(tlong,100*abs(BzPoly(order,:)-BzPoly(nOrderMax,:))./abs(BzPoly(nOrderMax,:)),'-','linewidth',1);
    hold all
end
title('Accuracy for various quadrature orders')
legend(gca,num2str([1:nOrderMax-1]'))
ylabel('Relative Difference (%)')
xlabel('Time (s)')
set(gca,'fontsize',14)