%% FUNCTION NAME: Main
% Main entry point function for key rate calculation.

% format long 
% clear all;
% close all;
% clc;

%%%%%%%%%%%%%%%%%%%%% Setting MATLAB Library Path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
restoredefaultpath
cd /Users/joseph/Documents
addpath(genpath('/Users/joseph/Documents/GitHub/EB_EAT'))
addpath(genpath('/Users/joseph/Documents/MATLAB/cvx'))
addpath(genpath('/Users/joseph/Documents/MATLAB/YALMIP-master'))

cvx_solver mosek

%%%%%%%%%%%%%%%%%%%%% Setting User Input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
preset='BB84';
[protocolDescription,channelModel,leakageEC,parameters,solverOptions]=feval(preset);

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
    protocolDescription.keyMap, ...
    protocolDescription.observables, ...
    channelModel.expectations, ...
    protocolDescription.krausOp, ...
    solverOptions.solver1)

%check validity of rho and perform perturbation if not valid
[rho,~]=perturbation_channel(rho);

%%%%%%Apply step2SolverAsymptotic to BB84Lossy



% %call step 2 solver
% N=numel(protocolDescription.observables);
% cons = zeros(1, N);
% for i =1:N
%     cons(i) = abs(real(trace(rho * protocolDescription.observables{i}) - channelModel.expectations(i)));
% end
% solverOptions.solver2.epsilonprime = max(cons);
% 
% [val,solver2Status] = ...
% step2SolverAsymptotic( ...
% rho, ...
% protocolDescription.observables, ...
% channelModel.expectations, ...
% [],[], ...
% protocolDescription.keyMap, ...
% protocolDescription.krausOp, ...
% solverOptions.solver2)