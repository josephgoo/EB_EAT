%% FUNCTION NAME: Step1
% Main entry point function for crossover min-tradeoff function calculation.

function [gval,debugInfo] = Step1(protocolDescription,channelModel,solverOptions)

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
        lowerBound = 0;
        FWBound = 0;
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
        
        [gval,solver2Status] = ...
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

end