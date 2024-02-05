%% FUNCTION NAME: Step2
% Main entry point function for min-tradeoff function calculation.

function mintradeoff=Step2(gval,zetaEp,names,p)

    %user should supply the list of parameters used in this description/channel file
    %this list varNames should be a subset of the full parameter list declared in the preset file
    %parameters specified in varNames can be used below like any other MATLAB variables
    varNames=["pz"];
    
    %%%%%%%%%%%%%%%%%%%%% interfacing (please do not modify) %%%%%%%%%%%%%%%%%%%%%%%%%
    
    %the functions findVariables and addVariables automatically search the input (names,p) for
    %the parameter values based on varNames, and convert them to MATLAB variables.
    varValues = findVariables(varNames,names,p);
    addVariables(varNames,varValues);
    gamma=(1-pz)^2;
    
    %%%%%%%%%%%%%%%%%%%%% user-supplied description begin %%%%%%%%%%%%%%%%%%%%%%%%%
    %maximum and minimum of g
    gMax=max(gval)-zetaEp;
    gMin=min(gval)-zetaEp;
    
    %coefficient of q(0) and q(1) in g
    gZero=gval(1)-zetaEp;
    gOne=gval(2)-zetaEp;

    %coefficient of q(0),q(1) and q(abort) in f (using eqs(52-53) in
    %george)
    fZero=gMax+(1/gamma)*(gZero-gMax);
    fOne=gMax+(1/gamma)*(gOne-gMax);
    fAbort=gMax;

    %Max(f), Min(f), Min_\Sigma(f), Var(f) (using eqs(54-57) in george)
    fMax=gMax;
    fMin=(1-1/gamma)*gMax+(1/gamma)*gMin;
    fMin_Sigma=gMin; %lower bound of Min_\Sigma(f)
    fVar=(1/gamma)*(gMax-gMin)^2; %upper bound of Var(f)

    %%%%%%%%%%%%%%%%%%%%% user-supplied description end %%%%%%%%%%%%%%%%%%%%%%%%%

    mintradeoff.fZero=fZero;
    mintradeoff.fOne=fOne;
    mintradeoff.fAbort=fAbort;
    mintradeoff.fMax=fMax;
    mintradeoff.fMin=fMin;
    mintradeoff.fMin_Sigma=fMin_Sigma;
    mintradeoff.fVar=fVar;
end





