function h = plot_model(S,x,lw,lc)
    if nargin<4
        lw = 1;
    end
%     if ~S.isotropic && nargin<3
%         beep;
%         disp ('specify H or V')
%     elseif S.isotropic
%         which = 'H';
%     end     
    numInt = length(x.z);
    earthmodel  = [S.z,x.z,S.rho(1,:),10.^x.rhoh]; %,S.rho(2,:),10.^x.rhov];
    S.numlayers = length(S.z)+numInt;
    z=earthmodel(1:S.numlayers);
    rho = [earthmodel(S.numlayers+1:2*S.numlayers)]; %;earthmodel(2*S.numlayers+1:3*S.numlayers)];
    plotz=[];
    plotz(1)=z(2);
    plotrhoh=[];plotrhov=[];
    for k=2:length(z)-1
        plotz(2*k-2)=z(k+1);plotz(2*k-1)=z(k+1);
        plotrhoh(2*k-3)=rho(1,k);plotrhoh(2*k-2)=rho(1,k);
        %plotrhov(2*k-3)=rho(2,k);plotrhov(2*k-2)=rho(2,k);
    end
    plotrhoh(2*k-1)=rho(1,end);%plotrhov(2*k-1)=rho(2,end);hold all
%         if (which=='V')
%             plot(log10([plotrhov';plotrhov(end)]),([plotz,S.zMax]),'-k','linewidth',lw)
%         else
       h =  plot(log10([plotrhoh';plotrhoh(end)]),([plotz,S.zMax]),lc,'linewidth',lw);
%     end
end    

