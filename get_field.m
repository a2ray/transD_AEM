

%get field for 1D model x 
function Bz = get_field(S,x)
    
% Isotropic only, ignoring rhov feld
ModelSig = 1./[S.rho(1,:),10.^x.rhoh ]; %hahaha blah blah blah
%ModelSig = 1./( 10.^x.rhoh );
ModelZ = [S.z,x.z];
%ModelZ = x.z;

mu  = ones(size(ModelSig));

    
if isfield(S,'data')  % Simple format
    
    Bz = get_LoopFields_TD_FHT(S.times,S.rTx,S.zTx,S.rRx,S.zRx,ModelSig,mu,ModelZ,...
                         S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.Ramp);        

else % SkyTEM High and Low Mode data:
    
    
    
        BzHigh = get_LoopFields_TD_FHT(S.HighMode.times,S.xyPolyTx,S.zTx,S.xyRx,S.HighMode.zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.HighMode.ramp,S.lowPassFilters);    

        BzLow = get_LoopFields_TD_FHT(S.LowMode.times,S.xyPolyTx,S.zTx,S.xyRx,S.LowMode.zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.LowMode.ramp,S.lowPassFilters);   
                     
     Bz = [BzHigh; BzLow];
   
end


    
    
end   