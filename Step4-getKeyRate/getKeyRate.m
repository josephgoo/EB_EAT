%% FUNCTION NAME: getKeyRate
% The core backend function for calculating key rate. This function takes
% protocol description, channel model (expectations) and error-correction leakage, and returns upper/lower bound and gap
% It calls solver1 and solver2 and returns key rate.
%
% Input Data Structure:
% protocol description: [keyMap,krausOperators,observables,obsMask(optional)]
% channel model: [expectations,expMask(optional)probDist/errorRate,pSift]
% EC description: [leakageEC]
% solverOptions: [globalSetting, solver1, solver2]
% (optional) parameter and name list: can be retrieved using findParameter("PARAMETER_NAME",solverOptions)
%
% Output Data Structure:
% upperBound, success flag
%%

function [upperBound,debugInfo]=getKeyRate(protocolDescription,channelModel,leakageEC,solverOptions,p,names)

    tstart_iteration = tic;
    try
        cvx_solver(solverOptions.globalSetting.cvxSolver);
        cvx_precision(solverOptions.globalSetting.cvxPrecision);
        if(solverOptions.globalSetting.verboseLevel==3)
            cvx_quiet(false);
        else
            cvx_quiet(true);
        end
    catch Error
        if(contains(Error.message,'Undefined function'))
           fprintf('cvx installation not found\n')
        end
        if(contains(Error.message,'missing solver'))
           fprintf('%s\n',Error.message)
        end
        upperBound = 0;
        debugInfo.exitStatus = 'cvx solver error';
        debugInfo.errorMessage = getReport(Error);
        fprintf("**** cvx solver setup error! ****\n")
        return
    end
    
    if(solverOptions.globalSetting.verboseLevel>=2)
       solverOptions.solver1.verbose = 'yes';
    else
       solverOptions.solver1.verbose = 'no';
    end
    
    %%%%%%%%%%%%%% call step 1 solver %%%%%%%%%%%%%%
    try
        solver1Status = [];
        rho0 = eye(prod(protocolDescription.dimensions));
        % ketPlus = 1/sqrt(2)*[1;1];
        % ketMinus = 1/sqrt(2)*[1;-1];
        % pz=findParameter("pz",names,p);
        % rho0=(pz)*kron(zket(4,1)*zket(4,1)',zket(2,1)*zket(2,1)')+...
        %     (pz)*kron(zket(4,2)*zket(4,2)',zket(2,2)*zket(2,2)')+...
        %     (1-pz)*kron(zket(4,3)*zket(4,3)',ketPlus*ketPlus')+...
        %     (1-pz)*kron(zket(4,4)*zket(4,4)',ketMinus*ketMinus');
        % rho0=rho0/2;
        [rho,fval,gap,solver1Status]=...
            step1Solver( ...
            rho0, ...
            protocolDescription.krausOp_sp, ...
            protocolDescription.observables, ...
            channelModel.expectations, ...
            protocolDescription.krausOp_p, ...
            solverOptions.solver1);
        
        %check validity of rho and perform perturbation if not valid
        [rho,~]=perturbation_channel(rho);
        
        debugInfo.solver1Status = solver1Status;
    catch Error
        debugInfo.exitStatus = 'step 1 solver error';
        debugInfo.errorMessage = getReport(Error);
        debugInfo.solver1Status = solver1Status;
        fprintf("**** step 1 solver error ****\n")
        fprintf("**** %s ****\n",Error.message)
        gval = 0;
        return
    end 
    
    %%%%%%%%%%%%%% call step 2 solver %%%%%%%%%%%%%%
    solver2Status = [];
    try
        N=numel(protocolDescription.observables);
        cons = zeros(1, N);
        for i =1:N
            cons(i) = abs(real(trace(rho * protocolDescription.observables{i}) - channelModel.expectations(i)));
        end
        solverOptions.solver2.epsilonprime = max(cons);
        
        [gval,zetaEp,solver2Status] = ...
            step2Solver( ...
            rho, ...
            protocolDescription.observables, ...
            channelModel.expectations, ...
            protocolDescription.krausOp_sp, ...
            protocolDescription.krausOp_p, ...
            solverOptions.solver2);
        
        debugInfo.solver2Status = solver2Status;
    catch Error
        debugInfo.exitStatus = 'step 2 solver error';
        debugInfo.errorMessage = getReport(Error);
        debugInfo.solver2Status = solver2Status;
        fprintf("**** step 2 solver error ****\n") 
        fprintf("**** %s ****\n",Error.message)
        gval=0;
        return
    end

    %%%%%%%%%%%%%% call step 2 - Min_tradeoff %%%%%%%%%%%%%%
    mintradeoff=Step2(gval,zetaEp,names,p);

    %%%%%%%%%%%%%% key length upper bound %%%%%%%%%%%%%%
    n=findParameter("n",names,p);
    pz=findParameter("pz",names,p);
    alphabet=findParameter("alphabet",names,p);
    beta=findParameter("beta",names,p);
    delta=findParameter("delta",names,p);
    dS=(alphabet+1)^2;
    alpha=(-beta+delta)/(-1+2*delta-beta*delta);
    gamma=(1-pz)^2;
    eps=findParameter("eps",names,p);

    %Q=findParameter("Q",names,p);
    %h=mintradeoff.fZero*gamma*Q+mintradeoff.fOne*gamma*(1-Q)+mintradeoff.fAbort*(1-gamma);

    V=sqrt(mintradeoff.fVar+2)+log2(2*dS^2+1);
    K=1/(6*(2-beta)^3*log(2))*2^((beta-1)*(log2(dS)+mintradeoff.fMax-mintradeoff.fMin_Sigma))*...
       (log(2^(log2(dS)+mintradeoff.fMax-mintradeoff.fMin_Sigma)+exp(1)^2))^3

    % V=log2(2*dS^2+1)
    % K=1/(6*(2-beta)^3*log(2))*2^((beta-1)*(log2(dS)+mintradeoff.fMax-mintradeoff.fMin))*...
    %     (log(2^(log2(dS)+mintradeoff.fMax-mintradeoff.fMin)+exp(1)^2))^3
    
    upperBound=n*mintradeoff.fMin-leakageEC-n*(((beta-1)*log(2))/2)*V^2-...
        n*(beta-1)^2*K-n*gamma*log2(abs(alphabet^2))+...
        ((beta-delta)/((beta+1)*(1-delta)))*log2(eps.acc)+...
        (alpha/(alpha-1))*(log2(eps.sec))+1;

    upperBound=upperBound/n;
end

%%%%%%%%%%% helper functions %%%%%%%%%%%%%%%%%

%select part of a cell/numerical array using another same-length mask array
%can be used to parse incoming observable (expectation) into groups with obsMask (expMask)
function a_m=applyMask(a,msk,label)
    if(length(a)~=length(msk))
        fprintf('**** mask length mismatch! ****\n')
    end
    a_m=[];
    for i=1:length(msk)
        if(msk(i)==label)
            a_m=[a_m;a(i)];
        end
    end
    
end

%retrieve a single parameter with a given "varName" label from the input parameter list p of this iteration
function varValue=findParameter(varName,names,p)
    
    varValues = findVariables(varName,names,p); %returns a cell array
    
    if (length(varValues)==0)
       fprintf('**** parameter %s not found! ****\n',varName); 
       varValue = [];
    else
       varValue = varValues{1}; %retrieve the parameter value (can be a numerical value or an array/matrix)
    end
end

%checks whether a single parameter with a given "varName" label exists in the input parameter list p of this iteration
function found=hasParameter(varName,names)
    
    found = false;

    for j=1:length(names)
        if(strcmp(varName,names(j)))
            found = true;
            break; %repeated entries will be ignored
        end
    end
    
end






