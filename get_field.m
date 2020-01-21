

%get field for 1D model x 
function Bz = get_field(S,x)
    
% Isotropic only, ignoring rhov feld
ModelSig = 1./[S.rho(1,:),10.^x.rhoh ]; 
ModelZ = [S.z,x.z];
mu  = ones(size(ModelSig));

Bz = [];
if S.modelLMpoly
    BzLow = get_LoopFields_TD_FHT(S.LM_times,S.xyPolyTx,S.zTx,S.xyRx,S.LM_zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.LM_ramp,S.lowPassFilters);
    Bz = [Bz;BzLow];         
end

if S.modelHMpoly
    BzHigh = get_LoopFields_TD_FHT(S.HM_times,S.xyPolyTx,S.zTx,S.xyRx,S.HM_zRx,ModelSig,mu,ModelZ,...
                 S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LoopQuadOrder,S.HM_ramp,S.lowPassFilters);
    Bz = [Bz;BzHigh];         
end
    
if S.modelLMloop
    if isfield(S,'LM_zTx') == 1
        S.zTx = S.LM_zTx;
    end
    if ~isfield(S,'LM_zRx')
        S.LM_zRx = S.zRx;
    end
    BzLow = get_LoopFields_circle_TD_FHT(S.LM_times,S.rTxLoop,S.zTx,S.rRx,S.LM_zRx,ModelSig,mu,ModelZ,...
                         S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.LM_ramp, S.lowPassFilters);
    Bz = [Bz;BzLow];
end
                        
if S.modelHMloop
    if isfield(S,'HM_zTx') == 1
        S.zTx = S.HM_zTx;
    end
    if ~isfield(S,'HM_zRx')
        S.HM_zRx = S.zRx;
    end
    BzHigh = get_LoopFields_circle_TD_FHT(S.HM_times,S.rTxLoop,S.zTx,S.rRx,S.HM_zRx,ModelSig,mu,ModelZ,...
                         S.HankelFilterName,S.CosSinFilterName,S.nFreqsPerDecade,S.HM_ramp, S.lowPassFilters);
    Bz = [Bz;BzHigh];         
end

end   