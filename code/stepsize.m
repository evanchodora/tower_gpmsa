% function params = stepsize(params,nburn,nlev)
% compute step sizes from step size data collect run in gpmmcmc
% please see associated documentation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Author: James R. Gattiker, Los Alamos National Laboratory
%         Brian Williams, Los Alamos National Laboratory
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

function params = stepsize(params,nburn,nlev)

pvars=params.mcmc.pvars;
svars=params.mcmc.svars;
wvars=params.mcmc.wvars;

% Grab some stuff to local variables
model=params.model; mcmc=params.mcmc;
p=model.p; q=model.q;
pu=model.pu;
lamVzGnum=model.lamVzGnum;
base=mcmc.base; acc=mcmc.acc;

if ~isfield(mcmc,'pacc');
   mcmc.pacc=1/exp(1);
end

% Setup for logistic regression
logit = log(mcmc.pacc)-log(1.0-mcmc.pacc);
vnburn = ones(nlev,1).*nburn;

jter = (-(nlev-1)/2:(nlev-1)/2)';
for varNum=1:length(svars)
  varName=svars{varNum};
  wvarName=wvars{varNum};
  varwidth=mcmc.(wvarName);
  for k=1:length(varwidth)
    if(varwidth(k) > 0.0)
      logs = log(varwidth(k));
      X = logs+jter.*log(base.(varName)(:,k));
      b = glmfit(X, [acc.(varName)(:,k) vnburn], 'binomial');
	  varwidth(k)=exp((logit-b(1))/b(2));
    end
  end
  %keyboard
  mcmc.(wvarName)=varwidth;
end


params.mcmc = mcmc;
