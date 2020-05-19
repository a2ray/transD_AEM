clear
rng(1); rhosd = 0.05;
% delz  : the z discretization
% rho0  : log10 resistivity of the surface
% rho1  : log10 resistivity of the anomaly top of the anomaly
% z1    : depth to top of anomaly
% delzA : thickness of anomaly
% z2    : depth of terminating halfspace
% rho2  : log10 resistivity at terminating halfspace 
[S,x,bg] = makemodel('delz',0.25,'rho0',log10(50),'rho1',log10(10),'z1',40,'rhoA',log10(50),'delzA',20,'rho2',log10(3),'z2',250);
jitteryx = x;
jitteryx.rhoh = x.rhoh + rhosd*randn(size(x.rhoh));

