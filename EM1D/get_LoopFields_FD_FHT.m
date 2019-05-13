
function Bz = get_LoopFields_FD_FHT(freqs,rTxLoop,zTx,rRx,zRx,sig,mu,z,filterName)
%
% Computes the vertical magnetic field frequency-domain response for a
% large loop source. 
%
% Usage:
%
% Bz = get_LoopFields_FD_FHT(freqs,rTxLoop,zTx,rRx,zRx,sig,mu,z,filterName)
%
%
% Inputs:
%
% freqs      - frequency(ies) [Hz]. Can be array of values or single value.
% rTxLoop    - radius of transmitter loop. [m]. single value
% zTx        - vertical position of transmitter loop (positive down). [m]
% rRx        - horizontal range(s) to the receiver(s). [m]. Can be array of values or single value.
% zRx        - vertical position of receiver(s) (positive down). [m]. Can be array of values or single value.
% sig        - array of conductivities for each layer. [S/m]
% mu         - array of relative magnetic permeabilities for each layer.
%              Normally this is just an array of ones.
% z          - vertical position of the top of each layer. [m]. Note that z
%              is positive down. Use a dummy value for the topmost layer. 
% filterName - Digital filter coefficients to use for the Hankel Transform.
%              Options: 'kk201Hankel.txt', 'kk101Hankel.txt', 'kk51Hankel.txt'
%              Use 'kk201Hankel.txt' to play it safe. Only use the shorter
%              101 or 51 point filters if you know what you are doing and
%              have proven that they are accurate for your setup and model
%              parameters.
% 
% Output:
%
% Bz         -  vertical magnetic field (T/Am^2). Dimensions: (length(freqs),length(rRx))
%
%
% Written by:
%
% Kerry Key
% Scripps Institution of Oceanography
%
%
%--------------------------------------------------------------------------

    mu0 = 4*pi*10^-7;

%
% Load the digital filter weights:
%
    A = load(filterName);
    Filter.base = A(:,1);
    Filter.J0   = A(:,2);
    Filter.J1   = A(:,3);
  
%
% Now compute the CSEM responses using the FHT method:
% 

    FJ0 = repmat(Filter.J0,1,length(freqs));
    FJ1 = repmat(Filter.J1,1,length(freqs));

    %
    % Loop over receivers
    %
    for iRx = 1:length(rRx)

        rRxEval = rRx(iRx);
        zRxEval = zRx(iRx);
        
        % FHT:
        if rRxEval < rTxLoop
            
            lTxClose    = true;
            lambda      = Filter.base/rTxLoop;    
            BzK         = getBzKernel(freqs,z,sig,mu,lambda,mu0,zTx,zRxEval,lTxClose,rRxEval,rTxLoop);     
            Bz(iRx,:)   = sum(BzK.*FJ1,1)/rTxLoop;   
            
        else
            lTxClose    = false;  
            
            % non-spline:
             lambda      = Filter.base/rRxEval; 
             BzK         = getBzKernel(freqs,z,sig,mu,lambda,mu0,zTx,zRxEval,lTxClose,rRxEval,rTxLoop);             
            % spline: % Testing only, could result in dicey responses...
%             lambdaSp = 10.^[min(log10(lambda))-.5:1/8:max(log10(lambda))+.5]'; 
%             BzK    = getBzKernel(freqs,z,sig,mu,lambdaSp,mu0,zTx,zRx,lTxClose,rRxEval,rTxLoop);
%             PP = spline(log10(lambdaSp),BzK.'); 
%             BzK  = ppval(PP,log10(lambda)).'; 
            Bz(iRx,:)   = sum(BzK.*FJ0,1)/rRxEval;  
            
        end


        % Normalize by dipole moment, apply pre-coefficients too          
        Bz(iRx,:)  =  rTxLoop*mu(1)*mu0*Bz(iRx,:)/ (pi*rTxLoop^2);      


    end % loop over receivers
        
 
%
%  All done, goodbye!
%
    return;
end
 
%- 
%--------------------------------------------------------------------------
function Bz = getBzKernel(freqs,z,sig,mu,lambda,mu0,zTx,zRx,lTxClose,rRxEval,rTxLoop) 
    
    epsilon = 8.8541878176d-12;
        
    dz        = diff(z);
    nz        = length(z);
    hh        = [1d60 dz(2:end) 1d60];
    
    iTxlayer = find(z < zTx,1,'last');
    iRxlayer = find(z < zRx,1,'last');
    
 
    if (iTxlayer == 1) 
        depthTx = 1d150; % large value allows for zeroing source term if in top layer
    else
        depthTx = zTx - z(iTxlayer);  
    end
 
    if (iTxlayer == length(z)) 
        heightTx = 1d150; % large value allows for zeroing source term if in bottom layer
    else
        heightTx = z(iTxlayer+1) - zTx;   
    end
    
    Rp    = complex(zeros(length(z),length(lambda),length(freqs)));
    Rm    = complex(zeros(length(z),length(lambda),length(freqs)));
 
    a     = complex(zeros(length(z),length(lambda),length(freqs)));
    b     = complex(zeros(length(z),length(lambda),length(freqs)));
 
    
    [LAM,SIG,FREQ] = meshgrid(lambda,sig,freqs);
    [~,MU,~]  = meshgrid(lambda,mu,freqs);
    [~,HH,~]  = meshgrid(lambda,hh,freqs);     
    
 
    gamma = sqrt(LAM.^2 - 1i*2*pi*FREQ*mu0.*MU.*(SIG - 1i*2*pi*FREQ*epsilon));
    
    %
    % Compute layer decay coefficients:
    %
    expgh = exp(-gamma.*HH);  
    
    %
    % Compute r+ and r- layer coefficients:
    %
    rp = ( MU(2:end,:,:).*gamma(1:end-1,:,:) - gamma(2:end,:,:).*MU(1:end-1,:,:) ) ./ ...
         ( MU(2:end,:,:).*gamma(1:end-1,:,:) + gamma(2:end,:,:).*MU(1:end-1,:,:) );
    rm = complex(zeros(size(rp,1)+1,size(rp,2),size(rp,3)));
    
    rm(2:length(sig),:,:) = -rp;
    
    %
    % Compute Rp recursions:  (working from bottom up to top)
    %        
    for i = length(sig)-1:-1:iTxlayer
        Rp(i+1,:,:) = Rp(i+1,:,:).*expgh(i+1,:,:);  % post applied to avoid positive exp in a,b formulation 
        Rp(i,:,:)   = ( rp(i,:,:) +        Rp(i+1,:,:).*expgh(i+1,:,:)) ./ ...
                      ( 1     + rp(i,:,:).*Rp(i+1,:,:).*expgh(i+1,:,:));
    end

    %
    % Compute Rm recursions:  (working from top down to bottom)
    %
    for i = 2:iTxlayer
        Rm(i-1,:,:) = Rm(i-1,:,:).*expgh(i-1,:,:);  % post applied to avoid positive exp in a,b formulation 
        Rm(i,:,:)   = ( rm(i,:,:) +       Rm(i-1,:,:).*expgh(i-1,:,:)) ./ ...
                      ( 1    + rm(i,:,:).*Rm(i-1,:,:).*expgh(i-1,:,:));
    end

    %
    % Compute potential coefficients in the src layer:
    %
    rmrp     = Rm(iTxlayer,:,:).*Rp(iTxlayer,:,:);
    onemrmrp = 1 - rmrp.*expgh(iTxlayer,:,:).*expgh(iTxlayer,:,:);  % added exp term since I omitted it earlier
    rhs      = 1./(2.d0*gamma(iTxlayer,:,:));
    srcp     = exp(-gamma(iTxlayer,:,:)*heightTx).*rhs;   
    srcm     = exp(-gamma(iTxlayer,:,:)*depthTx ).*rhs ;

    a(iTxlayer,:,:) = ( rmrp.*srcm.*expgh(iTxlayer,:,:)  + Rp(iTxlayer,:,:).*srcp  ) ./ onemrmrp;  
    b(iTxlayer,:,:) = ( rmrp.*srcp.*expgh(iTxlayer,:,:)  + Rm(iTxlayer,:,:).*srcm  ) ./ onemrmrp;    

    if iRxlayer < iTxlayer

        for i = iTxlayer-1:-1:iRxlayer

            aytop =  (a(i+1,:,:).*expgh(i+1,:,:) + b(i+1,:,:) + srcm ); 
            a(i,:,:)  = aytop ./ ( 1.d0 + Rm(i,:,:).*expgh(i,:,:) );  
            b(i,:,:)  = a(i,:,:).*Rm(i,:,:);
            srcm  = 0d0; % primary source term only for first layer above source layer

        end

    elseif iRxlayer > iTxlayer

        for i = iTxlayer+1:iRxlayer

            aytop = (a(i-1,:,:) + b(i-1,:,:).*expgh(i-1,:,:) + srcp );
            b(i,:,:)  = aytop./ (1.d0 + Rp(i,:,:).*expgh(i,:,:) );
            a(i,:,:)  = b(i,:,:).*Rp(i,:,:);
            srcp  = 0d0;  % primary source term only for first layer below source layer

        end
    end

    %
    % Finally, compute Fz kernel functions:
    %
    if iRxlayer < length(z)
        expp = exp(+gamma(iRxlayer,:,:)*(zRx - z(iRxlayer+1)));   
    else
        expp = complex(zeros(1,length(lambda),length(freqs)));
    end
    if iRxlayer > 1
        expm = exp(-gamma(iRxlayer,:,:)*(zRx - z(iRxlayer)));  
    else
        expm = complex(zeros(1,length(lambda),length(freqs)));
    end
 
    Fz = a(iRxlayer,:,:).*expp + b(iRxlayer,:,:).*expm;

    if iRxlayer == iTxlayer
        Fz = Fz + (exp(-gamma(iRxlayer,:,:).*abs(zRx - zTx))./(2*gamma(iRxlayer,:,:)));  
    end

    BzKernel   = squeeze(Fz.*LAM(1,:,:).^2);
    
    if size(BzKernel,1) == 1
        BzKernel = BzKernel.';
    end
    if lTxClose
        bess = besselj(0,lambda*rRxEval);
        Bz   =  BzKernel.*repmat( bess,1,length(freqs));
    else
        bess = besselj(1,lambda*rTxLoop);
        Bz   =  BzKernel.*repmat(bess,1,length(freqs));
    end         
 
    
end % 


 