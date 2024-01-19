%% FUNCTION NAME: primalDf
% This function calculates the gradient of primal problem objective
% function.
%%
% function dfval = primalDf(rho,krausOperators_sp,krausOperators_p)
% 
%     dfval=0;
%     for i=1:numel(krausOperators_sp)
%         spRho=krausOperators_sp{i}*rho*krausOperators_sp{i}';
%         %check validity of spRho and perform perturbation if not valid
%         [spRho,~]=perturbation_channel(spRho);
%         dfval=dfval-krausOperators_sp{i}'*logm(spRho)*krausOperators_sp;
%     end
% 
%     for i=1:numel(krausOperators_p)
%         pRho=krausOperators_p{i}*rho*krausOperators_p{i}';
%         %check validity of pRho and perform perturbation if not valid
%         [pRho,~]=perturbation_channel(pRho);
%         dfval=dfval+krausOperators_p{i}'*logm(pRho)*krausOperators_p;
%     end
% end

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