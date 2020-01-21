function h = plot_model(S,x,lw,lc)
    if nargin<4
        lc = 'k';
    end    
    if nargin<3
        lw = 2;
    end    
    
    stairs([-log10(S.rho(end)) x.rhoh], [S.z(end) x.z S.zMax], 'color', lc, 'linewidth', lw)
    hold on
end    

