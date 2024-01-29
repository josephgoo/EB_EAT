%% FUNCTION NAME: BB84Channel
% Simple channel model for entanglement-based BB84. Only an errorRate
% (misalignment) is considered. The expectations correspond to
% non-squashing model with four POVM outcomes.
%%

function channelModel = BB84Channel(protocolDescription,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["Q","pz"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %allocate the outputs of this channel model file
    %can be automatically filled in by calling addExpectations(x) or addExpectations(x,'mask',maskValue)
    expectations = [];
    expMask = [];
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model begin %%%%%%%%%%%%%%%%%%%%%%%%%
    
    dimA = protocolDescription.dimensions(1);
    dimB = protocolDescription.dimensions(2);
    
    ketPlus = 1/sqrt(2)*[1;1];
    ketMinus = 1/sqrt(2)*[1;-1];
    
    %Bell states
    phiPlus = 1/sqrt(2)*[1;0;0;1];
    phiMinus = 1/sqrt(2)*[1;0;0;-1];
    psiPlus = 1/sqrt(2)*[0;1;1;0];
    psiMinus = 1/sqrt(2)*[0;-1;1;0];

    %rho_sim constraints
    rho_sim = (1-3/2*Q)*phiPlus*phiPlus'+Q/2*(phiMinus*phiMinus'+psiPlus*psiPlus'+psiMinus*psiMinus');

    addExpectations(Q);
    addExpectations(1-Q);
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied channel model end %%%%%%%%%%%%%%%%%%%%%%%%%
    
    channelModel.expectations = expectations;
    %channelModel.errorRate = [ed,ed];
    %channelModel.pSift = [pz^2,(1-pz)^2];
end