%Preset input for entanglement-based BB84

function [protocolDescription,channelModel,leakageEC,parameters,solverOptions] = BB84()
    [protocolDescription,channelModel,leakageEC]=setDescription();
    parameters=setParameters();
    solverOptions=setOptions();
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% user input %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set description files
%returns function handles of description/channel/error-correction files (functions)
%string name must match file/function name; case-sensitive
function [protocolDescription,channelModel,leakageEC]=setDescription()
    
    protocolDescription=str2func('BB84Description');
    channelModel=str2func('BB84Channel');
    leakageEC=str2func('generalEC');

end

%set the input parameters
%three types of parameters are considered: scan, fixed, optimize
function parameters=setParameters()

    parameters.names = ["n","pz","fullstat","f","alphabet","eps","Q","beta","delta"]; %BB84 

    %%%%%%%%%%%%%%%% 1.parameter to scan over %%%%%%%%%%%%%%%%
    %must name at least one parameter to scan (can be a single-point array if only interested in a fixed value)
    
    %parameters.scan.eta = 10.^(-0.2*(0:5:0)/10); %channel transmittance
    % eta is changed from a single value to an array storing 5 values of
    % transmittances from 0km to 20km
    parameters.scan.n = 1e2:(1e10-1e2)/15:1e10; %(finite size) sent data
    

    %%%%%%%%%%%%%%%% 2.fixed parameters %%%%%%%%%%%%%%%%
    %optional; the constant values can be either numerical or array/matrices

    parameters.fixed.Q=0.005;%QBER
    %parameters.fixed.ed = 0.01; %misalignment, defined as single-photon error, i.e. sin^2(theta)
    parameters.fixed.pz = 0.99; %basis choice probability (for Z basis)
    %parameters.fixed.pd = 0; %1e-6; %dark count probability
    %parameters.fixed.etad = 1; %0.045; %detector efficiency
    parameters.fixed.fullstat = 1; %using full statistics or using QBER/Gain observables only
    parameters.fixed.f = 1.16; %error correction efficiency
    parameters.fixed.alphabet = 2; %(finite size) the encoding alphabet size - for qubits the size is 2

    parameters.fixed.beta=1.9;

    %security parameters for finite-size analysis
    eps.sec=(2/3)*1e-8;
    eps.EC=(1/3)*1e-8;
    eps.acc=eps.sec+eps.EC;
    parameters.fixed.eps = eps;

    %%%%%%%%%%%%%%%% 3.optimizable parameters %%%%%%%%%%%%%%%%
    %optional; declaring optimizable parameters automatically invokes local search optimizers
    %must be in the format of [lowerBound, initialValue, upperBound]
    %parameters.optimize.beta=[1,1+1/99999,2];
    parameters.optimize.delta=[0.9,0.9,0.99];
end

%set the running options for the solvers
function solverOptions=setOptions()

    %%%%%%%%%%%%%%%%% global setting %%%%%%%%%%%%%%%%%
    solverOptions.globalSetting.cvxSolver = 'mosek';
    solverOptions.globalSetting.cvxPrecision = 'high';
    
    %output level:
    %0.output nothing (except errors)
    %1.output at each main iteration
    %2.output at each solver 1 FW iteration
    %3.show all SDP solver output
    solverOptions.globalSetting.verboseLevel = 2; 
    
    %%%%%%%%%%%%%%%%% parameter optimizer setting %%%%%%%%%%%%%%%%%
    solverOptions.optimizer.name = 'coordinateDescent'; %choose between 'coordinateDescent' and 'bruteForce'
%     solverOptions.optimizer.name = 'localSearch_Adam';
    solverOptions.optimizer.linearResolution = 3; %resolution in each dimension (for brute force search and coordinate descent)
    solverOptions.optimizer.maxIterations = 1; %max number of iterations (only for coordinate descent)
    solverOptions.optimizer.linearSearchAlgorithm = 'iterative'; %choose between built-in 'fminbnd' and custom 'iterative' algorithms for linear search (only for coordinate descent)
    solverOptions.optimizer.iterativeDepth = 2; %choose depth of iteration levels; function count = depth * linearResolution (only for coordinate descent and if using 'iterative')
    solverOptions.optimizer.maxSteps = 5; %max number of steps (only for gradient descent and ADAM)
    solverOptions.optimizer.optimizerVerboseLevel = 2; %0:only output optimized result; 1:output a progress bar; 2:output at each function call

    %%%%%%%%%%%%%%%%% step1Solver setting %%%%%%%%%%%%%%%%%
    
    %options mainly affecting performance
    solverOptions.solver1.maxgap = 1e-6; %1e-6 for asymptotic, 2.5e-3 for finite;
    solverOptions.solver1.maxiter = 10;
    solverOptions.solver1.initmethod = 1; %minimizes norm(rho0-rho) or -lambda_min(rho), use method 1 for finite size, 2 for asymptotic v1
    
    %default options
    solverOptions.solver1.linearconstrainttolerance = 1e-10;
    solverOptions.solver1.linesearchprecision = 1e-20;
    solverOptions.solver1.linesearchminstep = 1e-3;
    solverOptions.solver1.maxgap_criteria = true; %true for testing gap, false for testing f1-f0
    solverOptions.solver1.removeLinearDependence = 'rref'; %'qr' for QR decomposition, 'rref' for row echelon form, empty '' for not removing linear dependence
    

    %%%%%%%%%%%%%%%%% step2Solver setting %%%%%%%%%%%%%%%%%
    solverOptions.solver2.epsilon = 0;
    solverOptions.solver2.epsilonprime = 1e-12;
    
    
end