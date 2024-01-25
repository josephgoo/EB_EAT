%% FUNCTION NAME: Step1
% Main entry point function for crossover min-tradeoff function calculation.

function [gval] = Step1(protocolDescription,channelModel,leakageEC,parameters,solverOptions)

    names=parameters.names;
    p={};
    for i=names
        if isfield(parameters.scan,i)
            p=[p,parameters.scan.(i)];
        elseif isfield(parameters.fixed,i)
            p=[p,parameters.fixed.(i)];
        end
    end
    
    
    protocolDescription=protocolDescription(names,p);
    channelModel=channelModel(protocolDescription,names,p);
    
    %call step 1 solver
    rho0 = eye(prod(protocolDescription.dimensions));
    [rho,fval,gap,status]=...
        step1Solver( ...
        rho0, ...
        protocolDescription.krausOp_sp, ...
        protocolDescription.observables, ...
        channelModel.expectations, ...
        protocolDescription.krausOp_p, ...
        solverOptions.solver1);
    
    %check validity of rho and perform perturbation if not valid
    [rho,~]=perturbation_channel(rho);
    
    %%%%%%Apply step2SolverAsymptotic to BB84
    
    
    
    %call step 2 solver
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

end