%% FUNCTION NAME: generalEC
% This function returns a function handle of error-correction leakage, for general formulation of protocols
% The return function depends on channelModel only, which contains QBER or probability distribution (depending on fine/coarse grain model) and probability of sifting.
%%

function leakageEC = generalEC(names,p)
    %this list varNames should be a subset of the full parameter list declared in the preset file
    varNames=["f","eps","n","Q"];
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    %from here on the parameters specified in varNames can be used like any other MATLAB variables
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %leakageEC = n*f*binaryEntropy(max(eps.EC,0.5))+log2(2/eps.EC);
    leakageEC = n*f*binaryEntropy(Q)+log2(2/eps.EC);
    
end

function [HX_Y, HY_X, IXY] = calculateEC(prob_dist)   	
    %a direct call with single return value will return H(X|Y)

    px_dist = sum(prob_dist,2);
    py_dist = sum(prob_dist)';
    
    HXY = - prob_dist(prob_dist > 0)' * log2(prob_dist(prob_dist > 0));
    HX = - px_dist(px_dist > 0)' * log2(px_dist(px_dist > 0));
    HY = - py_dist(py_dist > 0)' * log2(py_dist(py_dist > 0));
	
    HY_X = HXY - HX;
    HX_Y = HXY - HY;
    IXY  = HY - HY_X;
end