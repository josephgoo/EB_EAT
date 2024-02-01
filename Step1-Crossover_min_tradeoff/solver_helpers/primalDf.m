%% FUNCTION NAME: primalDf
% This function calculates the gradient of primal problem objective
% function.
%%

function dfval = primalDf(rho,krausOperators_sp,krausOperators_p)

    dfval=0;
    for i=1:numel(krausOperators_sp)
        ksp=krausOperators_sp{i};
        spRho=ksp*rho*ksp';
        %check validity of spRho and perform perturbation if not valid
        [spRho,~]=perturbation_channel(spRho);
        dfval=dfval-ksp'*logm(spRho)*ksp;
    end
    
    for i=1:numel(krausOperators_p)
        kp=krausOperators_p{i};
        pRho=kp*rho*kp';
        %check validity of pRho and perform perturbation if not valid
        [pRho,~]=perturbation_channel(pRho);
        dfval=dfval+kp'*logm(pRho)*kp;
    end
end