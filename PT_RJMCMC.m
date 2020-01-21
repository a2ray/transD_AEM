function PT_RJMCMC(DataFile,outputFolder,restart)

    %*** You need to specify these ***

    %data, frequencies used, depth of Rx, Tx, Azimuth, etc.
    S_0 = load (DataFile);
    
    %defaults which may not exist
    if isfield(S_0,'debug_prior') == 0
        S_0.debug_prior = false;
    end
    if isfield(S_0,'jeffereys_prior') == 0
        S_0.jeffereys_prior = false;
    end
    
    [~,FileRoot] = fileparts(DataFile);

%     %number of iterations
    N = S_0.numIterations;

    %Acceptance ratios in MCMC chains calculated every so many steps
    if isfield(S_0,'ARwindow') == 0
        ARwindow = 500;
    end
    
    %save every so many steps
    saveWindow = S_0.saveEvery;

    nDisplayEvery = 100;  % print message to screen

    %number of parallel chains (and temperatures)
    nTemps = S_0.nTemps;

    %inverse temperature B ladder
    B = logspace(-log10(S_0.Tmax),0,nTemps);

    %probability of swapping every count of the MCMC chain
    pSwap = 1;
        
    %step sizes in model space, decreasing temperature
    UstepSize = S_0.UstepSize; % update in log10 rho
    BstepSize = S_0.BstepSize; % birth/death in log10 rho
    MstepSize = S_0.MstepSize; %for move interface in m

    % Other parameters for the RJMCMC algorithm:

    S_0.isotropic = true;
    S_0.rhMin = S_0.log10rho_min;
    S_0.rvMin = S_0.log10rho_min;
    S_0.rhMax = S_0.log10rho_max;
    S_0.rvMax = S_0.log10rho_max;

    %*** Shouldn't need to modify below this ***

    if length(UstepSize) ~= nTemps || ...
        length(BstepSize) ~= nTemps || ...
        length(MstepSize) ~= nTemps
        beep
        disp('less steps sizes than nTemps')
        return
    end

    if exist(outputFolder) == 0
        mkdir(outputFolder);
    end

    AR = cell(nTemps,1);
    AR(:) = {cell(fix(N/ARwindow),1)};%Acceptance ratios
    DummyAR.uAR = 0; DummyAR.bAR = 0; DummyAR.dAR = 0; DummyAR.mAR = 0;
    DummyAR.TotalAR = 0;
    DummyAR.evalCount = 0;
    DummyAR.swapRate = 0;

    samples = cell(nTemps,1);
    samples(:) = {cell(N,1)};

    kTracker = cell(nTemps,1);
    kTracker(:) = {zeros(N,1)};

    en = zeros(N,2,nTemps);% Chi^2 and RMS errors

    swapCount = cell(nTemps,1);
    swapCount(:) = {zeros(N,1)};

    count = 0;
    S = cell(nTemps,1);
    S(:) = {S_0};
    
    if nargin<3
        rng('default')
        restart = false;
    else
        loadstruct = load([outputFolder,'/',FileRoot,'_PT_RJMCMC','_',num2str(nTemps)']);
        loadState = loadstruct.loadState;
        RandStream.setGlobalStream(loadState.Stream)
    end
    
    %initialize stuff and start point
    for ii=1:nTemps
        AR{ii}(:) = {DummyAR};
        ConvStat{ii}.uA = 0; ConvStat{ii}.bA = 0; ConvStat{ii}.dA = 0; ConvStat{ii}.mA = 0;
        ConvStat{ii}.uc = 0; ConvStat{ii}.bc = 0; ConvStat{ii}.dc = 0; ConvStat{ii}.mc = 0;
        ConvStat{ii}.evalCount = 0;
        if nargin<3
            k{ii} = 1;
            x{ii}.z = S_0.zMin + (S_0.zMax-S_0.zMin)*rand;
            x{ii}.rhov = S_0.rvMin + (S_0.rvMax-S_0.rvMin)*rand(1,k{ii}+1);
            if S{ii}.isotropic
                x{ii}.rhoh = x{ii}.rhov;
            else
                x{ii}.rhoh = S_0.rhMin + (S_0.rhMax-S_0.rhMin)*rand(1,k{ii}+1);
            end
        else
            % continue sampling with restart
            % Trim to what has actually been computed so far:
            loadstruct = load([outputFolder,'/',FileRoot,'_PT_RJMCMC','_',num2str(ii)']);
            iComputed = find(loadstruct.k_ll>0,1,'last');
            kTracker{ii}(1:iComputed) = loadstruct.k_ll(1:iComputed);
            samples{ii}(1:iComputed) = loadstruct.s_ll(1:iComputed);
            en(1:iComputed,:,ii) = loadstruct.en_ll(1:iComputed,:);
            x{ii} = samples{ii}{iComputed};
            k{ii} = kTracker{ii}(iComputed);
            assert(k{ii} == length(x{ii}.z));
            iComputed_AR = floor((iComputed-1)/ARwindow) + 1;
            AR{ii}(1:iComputed_AR) = loadstruct.AR_ll(1:iComputed_AR);
            if ii==1
                swapstruct = load ([outputFolder,'/',FileRoot,'_PT_RJMCMC','_swaps'],'swapCount');
            end
            swapCount{ii}(1:iComputed) = swapstruct.swapCount{ii}(1:iComputed);
            if ii == nTemps
                count = iComputed;
                clear loadstruct swapstruct
            end    
        end

        oldMisfit{ii} = getMisfit(x{ii},S_0);
        S{ii}.rSD1 = UstepSize(ii);
        S{ii}.rSD2 = BstepSize(ii);
        S{ii}.MoveSd = MstepSize(ii);
    end

    %start MCMC
    tic;
    tStart = tic;

    while count<N
        count = count +1;

        if mod(count,nDisplayEvery) == 0 % display text to user
           tLength      = toc(tStart);
           if restart
               aveIterRate = tLength/(count-iComputed);
           else
               aveIterRate     = tLength/count;
           end    
           predictedEnd = (N-count)*aveIterRate/86400 + now;
           fprintf('Iteration %i out of %i. Mean time per iteration: %4.2f s. Predicted completion time: %s\n',count,N,aveIterRate,datestr(predictedEnd))

        end
        %see if swap
        if rand<pSwap

            %then swap among ALL chains
            for iTemp = nTemps:-1:2
                for jTemp = rand(1:jTemp)
                    if iTemp ~= jTemp
                        %now find swap probability according to likelihoods
                        logAlphaSwap = (oldMisfit{iTemp}(1) - oldMisfit{jTemp}(1))*(B(iTemp) - B(jTemp));
                        if log(rand)<logAlphaSwap
                            temp_x      = x{iTemp};
                            temp_k      = k{iTemp};
                            temp_misfit = oldMisfit{iTemp};

                            x{iTemp}          = x{jTemp};
                            k{iTemp}          = k{jTemp};
                            oldMisfit{iTemp}  = oldMisfit{jTemp};

                            x{jTemp}         = temp_x;
                            k{jTemp}         = temp_k;
                            oldMisfit{jTemp} = temp_misfit;

                            swapCount{iTemp}(count) = 1;
                            swapCount{jTemp}(count) = 1;
                        end
                    end
                end % jTemp
            end % iTemp

        end % if pswap

        %start parallel tempering
        %one step
        for jj=1:nTemps

            [x{jj},k{jj},oldMisfit{jj},ConvStat{jj}] = RJ_MCMC_step(x{jj},k{jj},oldMisfit{jj},ConvStat{jj},S{jj},B(jj));
            samples{jj}{count} = x{jj};
            kTracker{jj}(count)= k{jj};
            en(count,:,jj)      = oldMisfit{jj};

            if mod(count,ARwindow) == 0
                 idx = count/ARwindow;
                 AR{jj}{idx}.uAR = ConvStat{jj}.uA/ConvStat{jj}.uc*100; AR{jj}{idx}.bAR = ConvStat{jj}.bA/ConvStat{jj}.bc*100;
                 AR{jj}{idx}.dAR = ConvStat{jj}.dA/ConvStat{jj}.dc*100; AR{jj}{idx}.mAR = ConvStat{jj}.mA/ConvStat{jj}.mc*100;
                 AR{jj}{idx}.TotalAR = (ConvStat{jj}.uA + ConvStat{jj}.bA + ConvStat{jj}.dA + ConvStat{jj}.mA)/ARwindow*100;
                 AR{jj}{idx}.evalCount = ConvStat{jj}.evalCount;
                 AR{jj}{idx}.swapRate  = sum(swapCount{jj}(count-ARwindow+1:count))/ARwindow*100;
                 ConvStat{jj}.uA = 0; ConvStat{jj}.bA = 0; ConvStat{jj}.dA = 0; ConvStat{jj}.mA = 0;
                 ConvStat{jj}.uc = 0; ConvStat{jj}.bc = 0; ConvStat{jj}.dc = 0; ConvStat{jj}.mc = 0;
            end
        end%one step

        %see if time to save
        if mod(count,saveWindow)==0
           for ll = 1:nTemps
            loadState.x{ll} = x{ll};
            loadState.Stream = RandStream.getGlobalStream;
            s_ll = samples{ll}; k_ll = kTracker{ll}; en_ll = en(:,:,ll); AR_ll = AR{ll}; S_ll=S{ll};
            save ([outputFolder,'/',FileRoot,'_PT_RJMCMC','_',num2str(ll)],'s_ll','k_ll','en_ll', 'AR_ll','S_ll','loadState')
           end
          save ([outputFolder,'/',FileRoot,'_PT_RJMCMC','_swaps'],'swapCount')
        end

    end
end


%This is for fixed dimension MCMC, which is commented out in the main
%program - either use this or RJ_MCMC_step in main program, not both
function [x,oldMisfit,accepted] = MCMCstep(x,oldMisfit,accepted,S,B)
    priorViolate=0;
    xNew = x;
    k = length(x.z);
    xNew.rhov = x.rhov + S.rSD1*randn(1,k+1);
%      layer = randi(k+1);
%      xNew.rhov(layer) = x.rhov(layer) + S.rSD1*randn;
    %is isotropic
    xNew.rhoh = xNew.rhov;

    if  any(xNew.rhov < S.rvMin) || any(xNew.rhov > S.rvMax)
        priorViolate=1;
    end

    if ~priorViolate
     newMisfit = getMisfit(xNew,S);
     logalpha = -B*(newMisfit(1) - oldMisfit(1))/2;
    else
        logalpha = -Inf;
    end

    if log(rand)<logalpha
       x=xNew;
       oldMisfit = newMisfit;
       accepted  = accepted +1;
    end
end
%RJ MCMC moves

function [pertNorm,xNew,priorViolate]=birth(k,x,S)
    xNew.z   = zeros(1,k+1) +NaN;%fills a new model with NaNs
    xNew.rhoh = zeros(1,k+1+1) +NaN;
    xNew.rhov = zeros(1,k+1+1) +NaN;
    pertNorm=0;
    priorViolate=0;

    %propose a new layer interface, at a non-existent interface location
    zProp = S.zMin + rand*(S.zMax-S.zMin);
    [~,pos]=ismember(0,(x.z(1:k)>=zProp),'legacy');

    mu       = [x.rhoh(pos+1),x.rhov(pos+1)];
    normDraw = randn(1,2);
    if S.isotropic%just to make sure no priorviolation takes place
	normDraw(1) = normDraw(2);%and pertnorm is isotropic
	mu(1) = mu(2);
    end
    newRho   = mu + normDraw*[S.rSD2 0;0 S.rSD2];
    if  newRho(1) < S.rhMin || newRho(1)>S.rhMax || ...
        newRho(2) < S.rvMin || newRho(2)>S.rvMax
        priorViolate=1;
        return %exit from this function with priorViolate=1
    end
    %positions for copying old values into new model vector
    xNew.z(pos+1)=zProp;
    xNew.z(isnan(xNew.z))=x.z;
    if rand<0.5
        pos=pos+1;
    end
    xNew.rhoh(pos+1)=newRho(1);
    xNew.rhoh(isnan(xNew.rhoh))=x.rhoh;
    xNew.rhov(pos+1)=newRho(2);
    xNew.rhov(isnan(xNew.rhov))=x.rhov;
    %now actually get the model pert norm sq
    pertNorm=sum((normDraw).^2);
end

function [xNew,priorViolate] = move (x,k,S)

    %pick an interface between 1 and current k to move
    l = randi(k);
    zProp = x.z(l) + S.MoveSd*randn;
    if zProp<S.zMin || zProp>S.zMax
        xNew=x;
        priorViolate=1;
        return; %outside bounds of z
    end

    priorViolate = 0;%checks passed!!
    %See where we are
    [~,pos]=ismember(0,(x.z(1:k)>=zProp),'legacy');
    if pos==l || pos == l-1
        %things are fine, nothing needs to be shifted
        xNew=x;
        xNew.z(l)=zProp;
        return;
    else %move the interface depth and properties
        x.z(l)  = [];%remove the interface
        k=k-1;
        if rand<0.5
            l=l+1;
        end
        %delete rho either above or below this interface
        %keep a copy of deleted value to move
        temp_rhoh    = x.rhoh(l);
        x.rhoh(l) = [];
        temp_rhov    = x.rhov(l);
        x.rhov(l) = [];
        %now birth the interface at zProp
        xNew.z   = zeros(1,k+1) +NaN;%fills a new model with NaNs
        xNew.rhoh = zeros(1,k+1+1) +NaN;
        xNew.rhov = zeros(1,k+1+1) +NaN;

        [~,pos]=ismember(0,(x.z(1:k)>=zProp),'legacy');
        %find where we are again just to be safe
        xNew.z(pos+1)=zProp;
        xNew.z(isnan(xNew.z))=x.z;

        if rand<0.5
            pos=pos+1;
        end
        xNew.rhoh(pos+1)=temp_rhoh;
        xNew.rhoh(isnan(xNew.rhoh))=x.rhoh;
        xNew.rhov(pos+1)=temp_rhov;
        xNew.rhov(isnan(xNew.rhov))=x.rhov;

        k=k+1;%not that this matters, we're not returning it.

    end

end

function [pertNorm,xNew] = death (x,k,S)
    xNew=x;
    %pick an interface to delete
    l = ceil(rand*(k));
    if S.isotropic%just to make sure H V perturbations are same
	x.rhoh(l)   = x.rhov(l);
	x.rhoh(l+1) = x.rhoh(l+1);
    end
    %these are the perturbations
    pertNorm = sum(([x.rhoh(l)-x.rhoh(l+1),x.rhov(l)-x.rhov(l+1)]/S.rSD2).^2);
    xNew.z(l)  = [];
    if rand<0.5
        l=l+1;
    end
    xNew.rhoh(l)= [];
    xNew.rhov(l)= [];
    end

function [xNew,priorViolate]=rhoUpdate(k,x,S,smallFlag)
    if strcmp(smallFlag,'small')
        S.rSD1 = S.rSD1/S.localFac;
    end
    xNew=x;
    priorViolate = 0;
    muH       = [x.rhoh];
    muV       = [x.rhov];
    normDrawH = randn(1,k+1);
    normDrawV = randn(1,k+1);
    newRhoH   = muH + normDrawH*(S.rSD1*eye(k+1));
    newRhoV   = muV + normDrawV*(S.rSD1*eye(k+1));
    if S.isotropic%just to make sure no priorviolation takes place
	newRhoH = newRhoV;
    end
    if  any(newRhoH < S.rhMin) || any(newRhoH>S.rhMax) || ...
        any(newRhoV < S.rvMin) || any(newRhoV>S.rvMax)
        priorViolate=1;
        return %exit from this function with priorViolate=1
    end
    xNew.rhoh=newRhoH;
    xNew.rhov=newRhoV;
end

%end RJMCMC moves

function [p,q] = determinPerm(n)
%returns index pairs p,q form strictly upper triangular matrix
    k = randperm(n/2*(n-1));
    q = floor(sqrt(8*(k-1) + 1)/2 + 3/2);
    p = k - (q-1).*(q-2)/2;
%[k;p;q]'
end

%RJ MCMC main step
function [x,k,oldMisfit,ConvStat] = RJ_MCMC_step (x,k,oldMisfit,ConvStat,S,B)
        del=S.rhMax-S.rhMin;
        %roll the dice
        dart=rand();
        CDF=[1/4,2/4,3/4,4/4];
        [~,pos]=ismember(0,(CDF>=dart),'legacy');
        pos=pos+1;%we know where it fell now
        priorViolate=0;
        switch pos
            case 1 %birth
                ConvStat.bc=ConvStat.bc+1;
                if k==S.kMax
                    priorViolate=1; %no birth allowed
                end
                if(~priorViolate)
                    [pertNorm,xNew,priorViolate]=birth(k,x,S);
                end
                %calculate acceptance probability
                if (~priorViolate)
                    newMisfit = getMisfit(xNew,S);
                    ConvStat.evalCount = ConvStat.evalCount+1;
                    %just priorR times propR
                    logalpha  = 2*log(S.rSD2)+log(2*pi)-2*log(del)+(pertNorm/2);
                    if S.isotropic%then priorR times propR is sqrtd
                        logalpha = 0.5*logalpha;
                    end
                    %insert likelihood to get alpha
                    if S.jeffereys_prior
                        logalpha  = logalpha + log(k) - log(k+1);
                    end    
                    logalpha  = logalpha -(newMisfit(1) - oldMisfit(1))*B;
                else
                    logalpha=-Inf;%prior has been violated so reject model
                end
                if log(rand)<logalpha %accept model
                    k=k+1;
                    x=xNew;
                    oldMisfit = newMisfit;
                    ConvStat.bA = ConvStat.bA +1;
                end

            case 2 %death
                ConvStat.dc=ConvStat.dc+1;
                if (k==S.kMin)
                    priorViolate=1; %no death at kMin
                end
                if (~priorViolate)
                    [pertNorm,xNew] = death (x,k,S);
                    %compute acceptance
                    newMisfit = getMisfit(xNew,S);
                    ConvStat.evalCount = ConvStat.evalCount +1;
                    %just priorR times propR
                    logalpha = -2*log(S.rSD2)-log(2*pi)+2*log(del)-(pertNorm/2);
                    if S.isotropic%then priorR times propR is sqrtd
                        logalpha = 0.5*logalpha;
                    end
                    %insert likelihood to get alpha
                    if S.jeffereys_prior
                        logalpha  = logalpha + log(k) - log(k-1);
                    end    
                    logalpha  = logalpha -(newMisfit(1) - oldMisfit(1))*B;
                else
                    logalpha=-Inf;
                end
                if log(rand)<logalpha %accept model
                    k=k-1;
                    x=xNew;
                    oldMisfit = newMisfit;
                    ConvStat.dA = ConvStat.dA +1;

                end
            case 3%move
                ConvStat.mc=ConvStat.mc+1;
                [xNew,priorViolate] = move (x,k,S);
                if ~priorViolate
                    newMisfit = getMisfit(xNew,S);
                    ConvStat.evalCount = ConvStat.evalCount +1;
                    logalpha = -(newMisfit(1) - oldMisfit(1))*B;
                else
                    logalpha=-Inf;
                end
                if log(rand)<logalpha %accept model
                    %k remains the same
                    x=xNew;
                    oldMisfit = newMisfit;
                    ConvStat.mA = ConvStat.mA +1;

                end
            case 4%update
                ConvStat.uc=ConvStat.uc+1;
                [xNew,priorViolate]=rhoUpdate(k,x,S,'large');%always large in PT
                if ~priorViolate
                    newMisfit = getMisfit(xNew,S);
                    ConvStat.evalCount = ConvStat.evalCount +1;
                    logalpha = -(newMisfit(1) - oldMisfit(1))*B;
                else
                    logalpha=-Inf;
                end
                if log(rand)<logalpha %accept model
                    %k remains the same
                    x=xNew;
                    oldMisfit = newMisfit;
                    ConvStat.uA  = ConvStat.uA +1;
                end
        end% all moves switch
        
end
%end RJ MCMC step
