

%get field for 1D model x 
function Bz = get_field(S,x)
    
% Isotropic only, ignoring rhov feld
ModelSig = 1./[S.rho(1,:),10.^x.rhoh ]; 
%ModelSig = 1./( 10.^x.rhoh );
ModelZ = [S.z,x.z];
%ModelZ = x.z;

% %for testing purposes only
%ModelZ = [ -1e5 0 50 100 ];
%ModelSig = [ 1e-12 1/100 1/10 1/1000 ];
%keyboard
mu  = ones(size(ModelSig));

    
if isfield(S,'data')  % Simple format
    
    Bz = get_LoopFields_TD_FHT(S.times,S.rTx,S.zTx,S.rRx,S.zRx,ModelSig,mu,ModelZ,...
                         S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.Ramp);        

else % SkyTEM High and Low Mode data:
    
    if( S.DataType == 1)
        
        BzLow = get_LoopFields_TD_FHT(S.LowMode.times,S.xyPolyTx,S.zTx,S.xyRx,S.LowMode.zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.LowMode.ramp,S.lowPassFilters);
             
        Bz = BzLow;
        %fprintf('We are in the low mode get_field section!\n')
             
    elseif( S.DataType == 2 )
    
        BzHigh = get_LoopFields_TD_FHT(S.HighMode.times,S.xyPolyTx,S.zTx,S.xyRx,S.HighMode.zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.HighMode.ramp,S.lowPassFilters); 
             
        Bz = BzHigh;
        %fprintf('We are in the high mode get_field section!\n')
        
    else
        
        BzLow = get_LoopFields_TD_FHT(S.LowMode.times,S.xyPolyTx,S.zTx,S.xyRx,S.LowMode.zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.LowMode.ramp,S.lowPassFilters);  
             
        BzHigh = get_LoopFields_TD_FHT(S.HighMode.times,S.xyPolyTx,S.zTx,S.xyRx,S.HighMode.zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.HighMode.ramp,S.lowPassFilters);
             
        Bz = [BzHigh; BzLow];
        %fprintf('We are in the high-low mode get_field section!\n')

    end
         
                         
end


    
    
end   