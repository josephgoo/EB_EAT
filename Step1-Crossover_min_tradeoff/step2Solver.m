%% FUNCTION NAME: step2Solver
% step 2 solver
% The solver solves for the dual problem and finds the lower bound for key rate
%%
%% Idea behind this code
%
%       The basic mathematical framework is based on arXiv:1710.05511
%   (https://arxiv.org/abs/1710.05511).
%   
%    We use the same notation as in the Appendix D, lemma 12. 
%
%   Note we implement a variant of the lemma 12 to take into account inequality constratins   
%
%
%% Syntax
%     [lowerbound, flag] = step2Solver(rho, observables, expectations,
%     krausOp_sp, krausOp_p, options)
%
% Input:
%
% *   rho - a density matrix rho_AB from step 1 calculation
%
% *   observables - a cell of Hermitian operators corresponding to equality
% constraints.
%
% *   expectations - a list of expectation values corresponding to observed statistics
%
% *   krausOp_sp - Joint POVM for secret and public informations
%
% *   krausOp_p- POVM for public information
%
% *   options - a structure that contains options for optimization
%%
% Outputs:
%
% *  lowerbound - the coefficient of the lower bound of "min-tradeoff function" (without error correction term, without rescaling back to log2)
%
% *  zetaEp - error term
%
% *  flag - a flag to indicate whether the optimization problem is solved
% successfully and whether we can trust the numerical result in the variable lowerbound
%
%%

function [lowerbound, zetaEp, status] = step2Solver(rho, observables, expectations, krausOp_sp, krausOp_p, options)
    
%     'REACHED STEP 2'
    warning('off','MATLAB:logm:nonPosRealEig');
    defaultOptions.epsilon = 0; % 0<epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
%    defaultOptions.epsilonprime = 1e-12; % epsilonprime is related to constraint tolerance
    
    if ~isfield(options,'epsilon')
        fprintf("**** solver 2 using default epsilon %f ****\n",defaultOptions.epsilon)
        options.epsilon = defaultOptions.epsilon;
    end
    % if ~isfield(options,'epsilonprime')
    %     fprintf("**** solver 2 using default epsilonprime %f ****\n",defaultOptions.epsilonprime)
    %     options.epsilonprime = defaultOptions.epsilonprime;
    % end
   
    %epsilonprime = options.epsilonprime;
    [fval, epsilon1] = primalfep(rho, krausOp_sp, krausOp_p, options);
    [gradf, epsilon2] = primalDfep(rho, krausOp_sp, krausOp_p, options); % epsilon1 == epsilon2 should hold
    fval = real(fval);
    dprime = size(rho, 1);
    epsilon = max(epsilon1, epsilon2);

    if epsilon> 1/(exp(1)*(dprime-1))
      ME = MException('step2Solver:epsilon too large','Theorem cannot be applied. Please have a better rho to start with');
      throw(ME);
    
    end
   

    % nConstraints = length(observables);
    % expectationsVector = [expectations+epsilonprime;-expectations+epsilonprime;expectations_ineq];
    % minusObservables = cell(nConstraints,1);
    % for i=1:nConstraints
    %    minusObservables{i} = -observables{i}; 
    % end
    % observablesVector = [observables;minusObservables;observables_ineq];
    [dualY, status] = submaxproblem(observables,expectations,gradf);
    
   
    if epsilon == 0
        zetaEp = 0;
    else 
       	%zetaEp = 2 * epsilon * (dprime-1) * log(dprime/(epsilon*(dprime-1)));
        zetaEp = 6 * epsilon * (dprime-1) * log2(dprime/(epsilon*(dprime-1)));
    end
    lowerbound = dualY;
    %pz = 0.9;
    %lowerbound = dualY/(1-pz)^2;
    %lowerbound = sum(expectations.*dualY);
    % note: we use natural log in all the calculation
    % one needs to convert natural log to log2. 
end


function [dualY, status] = submaxproblem(observables, expectations, gradf)
    nConstraints = length(expectations);
    totalDim = size(gradf, 1);
   
    cvx_begin sdp
        variable dualY(nConstraints) 
        maximize  sum(expectations  .* dualY)
        gradf - sdpCondition( dualY, observables) == hermitian_semidefinite(totalDim);
        %-dualY>=0; % dualY should have negative values
    cvx_end
    if strcmp(cvx_status, 'Infeasible') %| strcmp(cvx_status, 'Failed')
        fprintf("**** Warning: step 2 solver exception, submaxproblem status: %s ****\n",cvx_status);
        %checkValue = lambda_min(gradfTranspose - sdpCondition( dualY, GammaVector));
    else 
        %checkValue = lambda_min(gradfTranspose - sdpCondition( dualY, GammaVector));
    end
    
    %record status for debugging
    status = string(cvx_status);
end

function result = sdpCondition(dualY, observables)
   	result =0;
    for iConstraint = 1 : length(dualY)
        result = result + dualY(iConstraint) * observables{iConstraint};
    end  
end