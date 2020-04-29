clear
rng(1); rhosd = 0.05;

[S,x,bg] = makemodel('rho0',log10(50),'del',0.25,'rhoA',log10(100),'delzA',20,'rho2',log10(3),'z2',250);
jitteryx = x;
jitteryx.rhoh = x.rhoh + rhosd*randn(size(x.rhoh));

