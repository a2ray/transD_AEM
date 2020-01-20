
%Large Loop misfit function
function misfit = getMisfit(x,S)

if ~S.debug_prior
    
    % Get Bz for loop source:
    Bz = get_field(S,x);
    
    idx_LM = ~isnan(S.d_LM);
    idx_HM = ~isnan(S.d_HM);
    n_LM = sum(idx_LM);
    n_HM = sum(idx_HM);
    r_LM = abs(S.d_LM(idx_LM)) - abs(Bz(idx_LM));
    r_HM = abs(S.d_HM(idx_HM)) - abs(Bz(idx_HM));
    
    rms = 1.0;
    if (isfield(S, 'sd_LM') && isfield(S, 'sd_HM'))
        n = n_LM + n_HM;
        chi2by2 = 0.5*(norm(r_LM./S.sd_LM(idx_LM))^2 + ... 
                       norm(r_HM./S.sd_HM(idx_HM))^2);
        rms = 2*chi2by2/n;
    end
           
    if S.MLerroradjust
        r_LM = r_LM./abs(S.d_LM(idx_LM));
        r_HM = r_HM./abs(S.d_HM(idx_HM));
        chi2by2 = 0.5*(n_LM*log(r_LM'*r_LM) +...
                       n_HM*log(r_HM'*r_HM));
    end              
                  
    misfit = [chi2by2 rms];
else
    misfit = [0 0];
end    
        
end