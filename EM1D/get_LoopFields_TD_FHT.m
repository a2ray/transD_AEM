function dBzdt = get_LoopFields_TD_FHT(times,rTxLoop,zTx,rRx,zRx,sig,mu,z,...
                            HankelFilterName,CosSinFilterName,nFreqsPerDecade,rampTime,lowPassFilters)
%
% Computes the vertical magnetic field time-domain response for a
% large loop source. 
%
% Usage:
%
% dBzdt = get_LoopFields_TD_FHT(times,rTxLoop,zTx,rRx,zRx,sig,mu,z,...
%                            HankelFilterName,CosSinFilterName,nFreqsPerDecade,
%                            rampTime,lowPassFilters)
%
%
% Inputs:
%
% times      - time offsets for TDEM responses [s]. Can be array of values or single value.
% rTxLoop    - radius of transmitter loop. [m]. single value
% zTx        - vertical position of transmitter loop (positive down). [m]
% rRx        - horizontal range(s) to the receiver(s). [m]. Can be array of values or single value.
% zRx        - vertical position of receiver(s) (positive down). [m]. Can be array of values or single value.
% sig        - array of conductivities for each layer. [S/m]
% mu         - array of relative magnetic permeabilities for each layer.
%              Normally this is just an array of ones.
% z          - vertical position of the top of each layer. [m]. Note that z
%              is positive down. Use a dummy value for the topmost layer. 
% HankelFilterName - Digital filter coefficients to use for the Hankel Transform.
%                    Options: 'kk201Hankel.txt', 'kk101Hankel.txt', 'kk51Hankel.txt'
%                    Use 'kk201Hankel.txt' to play it safe. Only use the shorter
%                    101 or 51 point filters if you know what you are doing and
%                    have proven that they are accurate for your setup and model
%                    parameters.
% CosSinFilterName - Digital filter coefficients to use for the Fourier (Co)Sine Transforms.
%                    Options: 'kk201CosSin.txt', 'kk101CosSin.txt'
%                    Use 'kk201CosSin.txt' to play it safe. Only use the shorter
%                    101 point filters if you know what you are doing and
%                    have proven that they are accurate for your setup and model
%                    parameters.
% nFreqsPerDecade  - Number of frequencies per decade for sampling the
%                    frequency domain response prior to its time domain
%                    transformation. 10 is usually a safe value to use.
% rampTime         - Time and amplitude sequence of the transmitter waveform:
%                    [t0, a0; t1,a1; t2,a2;....tlast,0]
%                    example: [0 1; 1d-4 0] is a 1d-4s ramp off. More
%                    complicated waveforms with ramp-up at negative times
%                    are possible (for example SkyTEM waveform).
% lowPassFilters   - Array of corner frequencies[Hz] for any low-pass filters
%                    to apply. e.g. [ 4.500E+05  3.000E+05 ]
%
%
%
% Output:
%
% dBzdt         - time derivative of the vertical magnetic field V/(Am^4). Dimensions: (length(times),length(rRx))
%                 Note that dBzdt (V/m^2) is normalized by the dipole moment (Am^2) to give V/(Am^4) 
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
% Magic numbers:
%
    nTimesPerDecade     = 10;   % number of time samples per log 10 decade for the 
                                % TDEM spline interpolation when doing the waveform 
                                % ramp integration. 5 is usually sufficient.
    freqLowLimit        = 1d-3; % lower and upper limits for the sampling the FDEM response in log10(frequency)
    freqHighLimit       = 1d6;
    nQuadOrderTimeInteg = 5;    % Order of Guass quadrature used for time integration. 
                                 % This is the number of quadrature points per waveform time segment.
% 
% Parse inputs:
%
    if ~exist('nFreqsPerDecade','var')
        nFreqsPerDecade = 10;
    end
    if ~exist('lowPassFilters','var')
        lowPassFilters = [];
    end
    if ~exist('rampTime','var')
        rampTime = [];
    end
    
%
% Step 1: Get frequency domain fields for a broad sweep:
%

    freqs = 10.^(log10(freqLowLimit):1/nFreqsPerDecade:log10(freqHighLimit)); 

  % Perfectly circular loop:   
   BzFD = get_LoopFields_FD_FHT(freqs,rTxLoop,zTx,rRx,zRx,sig,mu,z,HankelFilterName);
    
  % Kernel for a point dipole - this still needs to have a Gauss quadrature wrapper for 
  % integrating over polygon of arbitrary loop source.
  % BzFD = get_VMD_FD_FHT(freqs,zTx,rRx,zRx,sig,mu,z,HankelFilterName);

  
%
% Step 1b: Apply front end filters, if input:
%
% Assumes lowPassFilters is arry of corner frequencies of first order
% Butterworth filters:

    if ~isempty(lowPassFilters)
       Hsc = 1;
       s = 1i*2*pi*freqs;
       for i = 1:length(lowPassFilters)
            fc = lowPassFilters(i);
            Hs = 1./( 1+s./(2*pi*fc));
            Hsc = Hsc.*Hs;
       end
        BzFD = BzFD.*conj(Hsc);
        
    end
 

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
% Initialize output dBzdt:
%
   
    dBzdt = zeros(length(times),length(rRx));
    
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
        % Now compute TD responses:
        %
        if isempty(rampTime) || size(rampTime,1) < 2
            
            % No waveform ramp, compute TDEM responses at specific requested times:
            
            % KWK note: if length(times) is really large, this could be
            % sped up by precomputing solution using
            % timesCompute = 10.^[min(log10(times))-1:1/nTimesPerDecade:max(log10(times))+1]; 
            % then interpolating to requested times. 
            % Will be more efficient when length(timesCompute) < length(times)
               
             for itime = 1:length(times)
            
                t      = times(itime);
                w      = log10(Filter.base/t);

                BzFpp  = ppval(PP,w); % get Bz at angular log10 freqs in w

                %kwk debug:BzFpp = get_LoopFields_FD_FHT(Filter.base/t/(2*pi),rTxLoop,zTx,rRx,zRx,sig,mu,z,HankelFilterName).';
               
                BzFpp  = -imag(BzFpp)*2/pi; % scale for impulse response

                %
                % Multiply kernels by filter weights and sum them all up,
                % remembering to divide by the particular time offset:
                %
                dBzdt(itime,iRx)  = sum(BzFpp.*Filter.SinFilt)/t;
                
             end
            
        else    
             % Waveform has a ramp, so compute TDEM response over even time
             % range, than do waveform integration with spline interpolation: 

            timesIn = times;
 
            times = 10.^[min(log10(timesIn))-1:1/nTimesPerDecade:max(log10(timesIn))+1]; 

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
          
            if nRamps == 1
                rT = rampTime;
                rampTime(1,1:2) = [0 1];
                rampTime(2,1:2) = [rT 0];
            end
            
            % quadrature weights:
            [x,w]=GLegIntP(nQuadOrderTimeInteg);
        
        
            % Loop over time then over ramp time steps:
            for itime = 1:length(timesIn)
                
                
                dBzdt(itime,iRx) = 0;
               
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
               
                    dBzdt(itime,iRx) = dBzdt(itime,iRx) + (b-a)/2*sum( f( (b-a)/2*x+(a+b)/2).*w)*dIdt;
                end

            end
             
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
 
