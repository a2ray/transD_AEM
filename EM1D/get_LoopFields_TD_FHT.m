function [Bz] = get_LoopFields_TD_FHT(times,rTx,zTx,rRx,zRx,sig,mu,z,...
                            HankelFilterName,CosSinFilterName,nFreqsPerDecade)
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

    freqs = 10.^[-6:1/nFreqsPerDecade:10]; 

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
        for itime = 1:length(times)
            
            tstart = tic;
            fct    = 0;
            t      = times(itime);
            w      = log10(Filter.base/t);
                         
            BzFpp  = ppval(PP,w); % get Bz at angular log10 freqs in w
            
            BzFpp  = -imag(BzFpp)*2/pi; % scale for impulse response
                    
            %
            % Multiply kernels by filter weights and sum them all up,
            % remembering to divide by the particular range:
            %
            Bz(itime,iRx) = sum(BzFpp.*Filter.SinFilt)/t;
         
                
        end
            
    end
 

%
%  All done, goodbye!
%
    return;
    
 
end



