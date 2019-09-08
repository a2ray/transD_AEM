
%Large Loop misfit function
function misfit = getMisfit(x,S)

if S.debug_prior ~= true
    
    % Get Bz for loop source:
    Bz = get_field(S,x);
    
    idx_LM = ~isnan(S.d_LM);
    idx_HM = ~isnan(S.d_HM);
    n_LM = sum(idx_LM);
    n_HM = sum(idx_HM);
    r = (abs([S.d_LM(idx_LM); S.d_HM(idx_HM)])...
       - abs(Bz([idx_LM; idx_HM])))./abs([S.d_LM(idx_LM); S.d_HM(idx_HM)]);
    
    n = n_LM + n_HM;
    chi2by2 = r'*r/2;
    rms = 2*chi2by2/n;
    
    if S.MLerroradjust
        r_LM = r(1:n_LM);
        r_HM = r(n_LM+1:end);
        chi2by2 = 0.5*(n_LM*log(r_LM'*r_LM) +...
                      n_HM*log(r_HM'*r_HM));  
    end
    
    misfit = [chi2by2 rms];
else
    misfit = [0 0];
end    
        
end