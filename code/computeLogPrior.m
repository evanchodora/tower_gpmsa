%function model = computeLogPrior(priors,model)
%
% Builds the prior likelihood
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: James R. Gattiker, Los Alamos National Laboratory
%
% This file was distributed as part of the GPM/SA software package
% Los Alamos Computer Code release LA-CC-06-079
%
% © Copyright Los Alamos National Security, LLC.
%
% This Software was produced under a U.S. Government contract 
% (DE-AC52-06NA25396) by Los Alamos National Laboratory, which is operated 
% by the Los Alamos National Security, LLC (LANS) for the U.S. Department 
% of Energy, National Nuclear Security Administration. The U.S. Government
% is licensed to use, reproduce, and distribute this Software. Permission
% is granted to the public to copy and use this Software without charge,
% provided that this Notice and any statement of authorship are reproduced 
% on all copies. Neither the Government nor the LANS makes any warranty, 
% express or implied, or assumes any liability or responsibility for the 
% user of this Software. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function model = computeLogPrior(priors,model)

logprior=0;

  % rhoU - Beta prior
      rhoU= exp(-model.betaU.*(0.5^2));
      rhoU(rhoU>0.999)=0.999;
      logprior=logprior - (1-priors.rhoU.b)*sum(log(1-rhoU(:)));

   % rhoV - Beta prior
      rhoV= exp(-model.betaV.*(0.5^2));
      rhoV(rhoV>0.999)=0.999;
      logprior=logprior - (1-priors.rhoV.b)*sum(log(1-rhoV(:)));

   % LamVz - Gamma prior  
      logprior = logprior + ...
          sum( (priors.lamVz.a-1)*log(model.lamVz) - ...
                priors.lamVz.b*model.lamVz );

   % LamUz - Gamma prior  
      logprior = logprior + ...
          sum( (priors.lamUz.a-1).*log(model.lamUz) - ...
                priors.lamUz.b*model.lamUz );

   % lamWs - Gamma prior  
      logprior = logprior + ...
          sum( (priors.lamWs.a-1).*log(model.lamWs) - ...
                priors.lamWs.b*model.lamWs );
   % lamWOs - Gamma prior  
      logprior = logprior + ...
          sum( (priors.lamWOs.a-1).*log(model.lamWOs) - ...
                priors.lamWOs.b*model.lamWOs );

   % lamOs - Gamma prior  
      logprior = logprior + ...
          (priors.lamOs.a-1).*log(model.lamOs) - ...
           priors.lamOs.b*model.lamOs;

   % theta - Normal prior
      if isfield(priors,'thetaPrior')
        logprior = logprior + feval(priors.thetaPrior,priors,model);
      else
        logprior = logprior - ...
          .5*sum(((model.theta-priors.theta.mean)./ ...
          priors.theta.std).^2);
      end

      model.logPrior=logprior;
