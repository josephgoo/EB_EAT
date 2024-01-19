%%  FUNCTION NAME: primalObjective
% This file contains the primal problem objective function.
%
% %
% $f(\rho)=sum_{s,p}H(\mathcal{K}_{sp}(\rho))-\sum_pH(\mathcal{K}_p(\rho))
%
% Syntax:  fval = primalf(rho,krausOperators_sp,krausOperators_p)
%
% Input: 
%
%  * rho  - density matrix shared by Alice and Bob
% 
%  * krausOperators_sp - Kraus operators for joint secret and public
%  information
%
%  * krausOperators_p - The Kraus operators for public information
%
% Output:
%  
%  * fval - the objective function value. 
%%

function fval = primalf(rho,krausOperators_sp,krausOperators_p)

%check validity of rho and perform perturbation if not valid
[rho,~]=perturbation_channel(rho);

    fval=0;
    for i=1:numel(krausOperators_sp)
        spRho=krausOperators_sp{i}*rho*krausOperators_sp{i}';
        %check validity of spRho and perform perturbation if not valid
        [spRho,~]=perturbation_channel(spRho);
        fval=fval-real(trace(spRho*logm(spRho)));
    end

    for i=1:numel(krausOperators_p)
        pRho=krausOperators_p{i}*rho*krausOperators_p{i}';
        %check validity of pRho and perform perturbation if not valid
        [pRho,~]=perturbation_channel(pRho);
        fval=fval+real(trace(pRho*logm(pRho)));
    end

end