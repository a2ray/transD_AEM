
clear all
clc

lineNum = 3275;
FID = 2340;

FileName = ['SkyTEMinversions/Trash/SkyTEM-' num2str(lineNum) '-FID' num2str(FID) '_PT_RJMCMC'];

U1 = load([FileName '_6']);
U2 = load([FileName '_7']);
U3 = load([FileName '_8']);

j1 = find(U1.k_ll,1,'last');
j2 = find(U2.k_ll,1,'last');
j3 = find(U3.k_ll,1,'last');

burnIn = 4000;


s_ll = [U1.s_ll(burnIn:j1); U2.s_ll(burnIn:j2); U3.s_ll(burnIn:j3)];
k_ll = [ U1.k_ll(burnIn:j1) ; U2.k_ll(burnIn:j2) ; U3.k_ll(burnIn:j3) ];

save([FileName,'_Combined'],'s_ll','k_ll')





