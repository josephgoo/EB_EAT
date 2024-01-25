%% FUNCTION NAME: primalDfep
%  Calculate $\nabla f_{\epsilon}(\rho)$ function, where $\epsilon$ value is
%  carefully chosen. 
%%
function [Dfval,realEpsilon] = primalDfep(rho, krausOperators_sp, krausOperators_p, options)
  
%     defaultoptions.epsilon = 0; % 0<=epsilon<=1/(e(d'-1)), epsilon is related to perturbed channel
%     defaultoptions.perturbation = 1e-16; % a small value added to the minimum eigenvalue;
%     if nargin == 4
%         if ~isfield(options,'epsilon')
%             options.epsilon = defaultoptions.epsilon;
%         end
%         if ~isfield(options,'perturbation')
%             options.perturbation = defaultoptions.perturbation;
%         end
%     else 
%         options = defaultoptions;
%     end



    [rho,epsilon1] = perturbation_channel(rho);
    Dfval=0;
    epsilon2=[];
    for i=1:numel(krausOperators_sp)
        ksp=krausOperators_sp{i};
        spRho=ksp*rho*ksp';
        %check validity of spRho and perform perturbation if not valid
        [spRho,epsilon]=perturbation_channel(spRho);
        Dfval=Dfval-ksp'*logm(spRho)*ksp;
        epsilon2=[epsilon2,epsilon];
    end
    
    epsilon3=[];
    for i=1:numel(krausOperators_p)
        kp=krausOperators_p{i};
        pRho=kp*rho*kp';
        %check validity of pRho and perform perturbation if not valid
        [pRho,epsilon]=perturbation_channel(pRho);
        Dfval=Dfval+kp'*logm(pRho)*kp;
        epsilon3=[epsilon3,epsilon];
    end
   
    realEpsilon = max([epsilon1,epsilon2,epsilon3]);
    
end