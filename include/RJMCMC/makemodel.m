function [S,x,bg] = makemodel(varargin)
    % make a three region stairstep model which is easy to model

    %# define defaults at the beginning of the code so that you do not need to
    opt = struct('delz',1,'rho0',log10(100),'rho1',log10(10),'z1',40,'rhoA',log10(50),'delzA',20,'rho2',log10(1),'z2',100);

    %# read the acceptable names
    optionNames = fieldnames(opt);

    %# count arguments
    nArgs = length(varargin);
    if round(nArgs/2)~=nArgs/2
       error('this function needs propertyName/propertyValue pairs')
    end

    for pair = reshape(varargin,2,[]) %# pair is {propName;propValue}
       inpName = pair{1};

       if any(strcmp(inpName,optionNames))
          %# overwrite options. If you want you can test for the right class here
          %# Also, if you find out that there is an option you keep getting wrong,
          %# you can use "if strcmp(inpName,'problemOption'),testMore,end"-statements
          opt.(inpName) = pair{2};
       else
          error('%s is not a recognized parameter name',inpName)
       end
    end
    
    % function starts here
    
    % set up z, rho
    z = 0:opt.delz:opt.z2;
    rho = zeros(size(z));
    % overburden
    idx1 = z<=opt.z1;
    rho(idx1) = opt.rho0 + (opt.rho1-opt.rho0)*z(idx1)/opt.z1;
    % anomaly
    idx2 = z>opt.z1 & z<=opt.z1 + opt.delzA;
    rho(idx2) = opt.rhoA;
    % underburden?
    idx3 = z>opt.z1 + opt.delzA & z<=opt.z2;
    % resistivity down to bottom of anolmaly if we had continued linear trend
    rdown =  opt.rho0 + (opt.rho1-opt.rho0)*(opt.z1+opt.delzA)/opt.z1;
    rho(idx3) = rdown + (opt.rho2-rdown)*(z(idx3)-opt.z1-opt.delzA)/(opt.z2-opt.z1-opt.delzA);
    
    % a background model
    rhobg = zeros(size(z));
    idx1 = z<=opt.z1 + opt.delzA;
    rhobg(idx1) = opt.rho0 + (opt.rho1-opt.rho0)*z(idx1)/opt.z1;
    rhobg(~idx1) = rho(~idx1);
    
    S = [];
    S.z = [-1d6, z(1)];
    znew = 0.5*(z(1:end-1) + z(2:end));
    S.zMax = znew(end) + 1;
    S.rho = [1d12];
    x.rhoh = rho;
    x.z = znew;
    
    bg.rhoh = rhobg;
    bg.z = znew;
end