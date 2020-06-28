function h = plot_model(S,x,lw,lc)
    if nargin<4
        lc = 'k';
    end    
    if nargin<3
        lw = 2;
    end    
    if ~isfield(S, 'log10rho_max')
        S.log10rho_max = max(x.rhoh)+eps();
        S.log10rho_min = min(x.rhoh)-eps();
    end
    stairs(-[log10(S.rho(end)) x.rhoh], [S.z(end) x.z S.zMax], 'color', lc, 'linewidth', lw)
    xlim([-S.log10rho_max -S.log10rho_min])
    xlabel('Log10 S/m')
    hold on
end    

