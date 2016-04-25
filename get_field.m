
%get field for 1D model x 
function Bz = get_field(S,x)
    
    %  [Er,Eb,Hr,Hb,Ez,Hz] = get_CSEM1D_FD_FHT_aniso_hed_mex(S.f,S.r,S.zRx,S.zTx,[S.z,x.z],1./ModelRho,'kk201Hankel.mat',1);
    
    % Isotropic only, ignoring rhov feld
    ModelSig = 1./[S.rho(1,:),10.^x.rhoh ]; 
    ModelZ = [S.z,x.z];
    
    mu  = ones(size(ModelSig));
    
    Bz = get_LoopFields_TD_FHT(S.times,S.rTx,S.zTx,S.rRx,S.zRx,ModelSig,mu,ModelZ,...
                         S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade);        
    
end   