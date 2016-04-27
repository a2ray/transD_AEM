function [Bz] = get_LoopFields_TD_FHT(times,rTx,zTx,rRx,zRx,sig,mu,z,...
                            HankelFilterName,CosSinFilterName,nFreqsPerDecade,rampTime)
%
% Notes:
%           rampTime - (s) single number of ramp-off time or piecewise linear model
%                       of [time, normalized_current] where the first row should be
%                       [0,1] and last row should be [tlast, 0]. For
%                       example for a ramp down over 0.00001 s:
%                      rampTime  = [0,1;  0.00001,0.6;  0.00002,0.2; 0.00001,0.0];
%           *** Assumes time 0 is start of ramp. ****
%
% KWK debug: Help text needs to be updated for new code...
% 
%
% Written by:
%
% Kerry Key
% Scripps Institution of Oceanography
%
%
%--------------------------------------------------------------------------
 
          
%
% Step 1: Get frequency domain fields for a broad sweep:
%

    freqs = 10.^[-3:1/nFreqsPerDecade:6]; 

    [BzFD] = get_LoopFields_FD_FHT(freqs,rTx,zTx,rRx,zRx,sig,mu,z,HankelFilterName);

 
%
% Step 2: Tranform to time domain:
%
%
% Load the time domain transform digital filter weights here:
%
    A = load(CosSinFilterName);
    Filter.base    = A(:,1);
    Filter.CosFilt = A(:,2);
    Filter.SinFilt = A(:,3);
           
%
% Initialize output Bz:
%
   
    Bz = zeros(length(times),length(rRx));
    
    log10omega = log10(2*pi*freqs); % compute splines in log10(omega) domain for stability
 
%
% Loop over receivers
%
    for iRx = 1:length(rRx)
        
        rRxEval = rRx(iRx);
        zRxEval = zRx(iRx);
            
        %
        % Compute the spline coefficients of the frequency domain solution
        % for this receiver:
        %
        PP = spline(log10omega,BzFD(iRx,:).');
        
        %
        % Loop over each requested time offset:
        %
         
        
        timesIn = times;

        nTimesPerDecade = 5;
        times = 10.^[min(log10(timesIn))-1:1/nTimesPerDecade:max(log10(timesIn))+.1]; 

        BzP = zeros(length(times),1);
        
        for itime = 1:length(times)
            
            t      = times(itime);
            w      = log10(Filter.base/t);
                         
            BzFpp  = ppval(PP,w); % get Bz at angular log10 freqs in w
            
            BzFpp  = -imag(BzFpp)*2/pi; % scale for impulse response
            
            %
            % Multiply kernels by filter weights and sum them all up,
            % remembering to divide by the particular range:
            %
            BzP(itime) = sum(BzFpp.*Filter.SinFilt)/t;
         
                
        end
        
        PP = spline(log10(times),BzP);  
        
        nRamps = size(rampTime,1);
        

        if nRamps > 1 || rampTime(1) > 0

                
            if nRamps == 1
                rT = rampTime;
                rampTime(1,1:2) = [0 1];
                rampTime(2,1:2) = [rT 0];
            end
            
            % quadrature weights:
            N = 5;
            [x,w]=GLegIntP(N);
        
        
            % Loop over time then over ramp time steps:
            for itime = 1:length(timesIn)
                
                
                Bz(itime,iRx) = 0;
               
                for iRamp = 1:size(rampTime,1)-1   
                    
                    % Convolve ramp with v(t):
                    
                    rta = rampTime(iRamp  ,1);
                    rtb = rampTime(iRamp+1,1);
                    dt = rtb - rta;
                    dI = rampTime(iRamp+1,2) - rampTime(iRamp,2);
                    dIdt = dI/dt;

                    if rta > timesIn(itime)
                        break
                    end
                    if rtb > timesIn(itime)  % end in this interval
                        rtb = timesIn(itime);
                    end
                    
                    
                    ta = timesIn(itime)-rta;
                    tb = timesIn(itime)-rtb;
                    
                    tb = max(tb,1d-8);  % rtb > rta, so make sure this is not zero...
                    
                    a = log10(ta);
                    b = log10(tb);
                    f = @getRampResp;
                    %Bz(itime,iRx) = Bz(itime,iRx) + quadgk(f,a,b,'AbsTol',1e-20,'RelTol',1e-4)*dIdt;  
                    
                    Bz(itime,iRx) = Bz(itime,iRx) + (b-a)/2*sum( f( (b-a)/2*x+(a+b)/2).*w)*dIdt;
                end

            end

        else  % No ramp time needed
               
            Bz(:,iRx) = ppval(PP,log10(timesIn));
        end

    end
 

%
%  All done, goodbye!
%
    return;
    
    % nested function for evaluating Bz(log10) integration kernel:
    function resp = getRampResp(t)
          
            resp = ppval(PP,t).*10.^t*log(10);
           
    end
 
end

function [x,w]=GLegIntP(nGL)
% Gauss-Legendre integration points in the interval [-1,1].
%
% Description
%     [#x#,#w#]=GLegIntP(#nGL#)
%     calculates the abscissas (x) and weights (w) for the Gauss-Legendre
%     quadrature. The algorithm presented here is described in J. A.
%     Gubner, Gaussian quadrature and the eigenvalue problem, October 16,
%     2014.
%
% Input arguments
%     #nGL# (scalar) is the number of the integration points to be
%     calculated.
%
% Output arguments
%     #x# ([#nGL# x 1]) contains the coordinates of the integration points.
%     #w# ([#nGL# x 1]) contains the weights of the integration points.
%
% Parents (calling functions)
%     GLegQuad > GLegIntP
%
% Children (called functions)
%     GLegIntP >
%

%__________________________________________________________________________
% Contact author
%
%  (c) 2015 by George Papazafeiropoulos
%  First Lieutenant, Infrastructure Engineer, Hellenic Air Force
%  Civil Engineer, M.Sc., Ph.D. candidate, NTUA
%
% Email: gpapazafeiropoulos@yahoo.gr
%
% Website: http://users.ntua.gr/gpapazaf/
%

beta=(1:nGL-1)./sqrt(4*(1:nGL-1).^2-1);
J=diag(beta,-1)+diag(beta,1);
[V,D]=eig(J);
[x,ix]=sort(diag(D));
w=2*V(1,ix)'.^2;

end
 
