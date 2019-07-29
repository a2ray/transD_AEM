
%Large Loop misfit function
function misfit = getMisfit(x,S)

if S.debug_prior ~= true
    
    % Get Bz for loop source:
    Bz = get_field(S,x);
    
    idx_LM = ~isnan(S.d_LM);
    idx_HM = ~isnan(S.d_HM);
    
    r = (abs([S.d_LM(idx_LM); S.d_HM(idx_HM)])...
       - abs(Bz([idx_LM; idx_HM])))./[S.sd_LM(idx_LM); S.sd_HM(idx_HM)];
    
    n = sum(idx_LM) + sum(idx_HM);
    chi2by2 = r'*r/2;
    misfit = [chi2by2 sqrt(2*chi2by2/n)];
else
    misfit = [0 0];
end    
        
end