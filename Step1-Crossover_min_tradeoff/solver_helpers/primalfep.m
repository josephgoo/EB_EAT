%% FUNCTION NAME: primalfep
%  Calculate $f_{\epsilon}(\rho)$ function, where $\epsilon$ value is
%  carefully chosen. 
%%

function [fval,realEpsilon] = primalfep(rho,krausOperators_sp,krausOperators_p,options)
%     
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
    fval=0;
    epsilon2=[];
    for i=1:numel(krausOperators_sp)
        spRho=krausOperators_sp{i}*rho*krausOperators_sp{i}';
        %check validity of spRho and perform perturbation if not valid
        [spRho,epsilon]=perturbation_channel(spRho);
        fval=fval-real(trace(spRho*logm(spRho)));
        epsilon2=[epsilon2,epsilon];
    end

    epsilon3=[];
    for i=1:numel(krausOperators_p)
        pRho=krausOperators_p{i}*rho*krausOperators_p{i}';
        %check validity of pRho and perform perturbation if not valid
        [pRho,epsilon]=perturbation_channel(pRho);
        fval=fval+real(trace(pRho*logm(pRho)));
        epsilon3=[epsilon3,epsilon];
    end

    realEpsilon=max([epsilon1,epsilon2,epsilon3]);
end