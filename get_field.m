

%get field for 1D model x 
function Bz = get_field(S,x)
    
% Isotropic only, ignoring rhov feld
ModelSig = 1./[S.rho(1,:),10.^x.rhoh ]; 
ModelZ = [S.z,x.z];

mu  = ones(size(ModelSig));

    
if isfield(S,'data')  % Simple format
    
    Bz = get_LoopFields_TD_FHT(S.times,S.rTx,S.zTx,S.rRx,S.zRx,ModelSig,mu,ModelZ,...
                         S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.Ramp);        

else % SkyTEM High and Low Mode data:
    
    
    
        BzHigh = get_LoopFields_TD_FHT(S.HighMode.times,S.HighMode.rTx,S.HighMode.zTx,S.HighMode.rRx,S.HighMode.zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.HighMode.ramp);    

        BzLow = get_LoopFields_TD_FHT(S.LowMode.times,S.LowMode.rTx,S.LowMode.zTx,S.LowMode.rRx,S.LowMode.zRx,ModelSig,mu,ModelZ,...
                S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LowMode.ramp);   
                     
     Bz = [BzHigh; BzLow];
   
end


    
    
end   