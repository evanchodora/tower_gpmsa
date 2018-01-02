% function params = gpmmcmc(params,nmcmc,varargin)
%  nmcmc - number of full draws to perform (overridden for stepInit mode)
%  varargs are in string/value pairs
%    'linkedMod' - default 0, 1 ==> process a joint model calibration
%    'noInint'   - default 0, 1 ==> do not init the model
%    'noCounter' - default 0, 1 ==> do not output a counter of interations
%    'initOnly'  - default 0, 1 ==> perform model initialization and return
%    'step'      - default 0, 1 ==> specified step size mode
%    'stepInit'  - default 0, 1 ==> adaptive step size initialization mode
%    'nBurn'     - stepInit mode, number of draws
%    'nLev'      - stepInit mode, number of candidate levels

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

function params = gpmmcmc(params,nmcmc,varargin)

% Process input arguments
linkedMod=0; noInit=0; noCounter=0; initOnly=0;
stepInit=0; step=0; nBurn=0; nLev=0;
parseAssignVarargs({'linkedMod','noInit','noCounter','initOnly', ...
  'step','stepInit','nBurn','nLev'});

% Grab some stuff to local variables
model=params.model; data=params.data; priors=params.priors;
mcmc=params.mcmc;
p=model.p;   q=model.q;
pu=model.pu; pv=model.pv;
lamVzGnum=model.lamVzGnum;
pvars=params.mcmc.pvars; 
svars=params.mcmc.svars; svarSize=params.mcmc.svarSize;
wvars=params.mcmc.wvars; 

% Initialize the models if we are not in a linked model model
if ~noInit
  C.var='all';
  model=computeLogLik(model,data,C); % compute all partial factors
  model=computeLogPrior(params.priors,model);
  model.logPost=model.logPrior+model.logLik;
end
% if we're in initOnly mode, return now
if initOnly
  params.model=model;
  return
end

% set the correct nmcmc value
if stepInit & ~linkedMod; nmcmc=nBurn*nLev; end

% initialize the structure that will record draw info
if numel(params.pvals);
  pvals=params.pvals; poff=length(pvals);
else
  for var=pvars; pvals(1).(var{1})=0; end; poff=0;
end
pvals(poff+nmcmc)=pvals(1);

% Counter will be used and displayed if we are not in linked model mode
if ~noCounter; counter('stime',1,nmcmc,5,10); end

if stepInit
  if ~linkedMod, mcmc.base=[]; end
  mcmcStepInit;
elseif step
  mcmcStep;
else
  mcmcAdaptive;
end

if ~noCounter; counter('end'); end
% And that's it for the main function ....

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function mcmcAdaptive
    % Do a standard mcmc
    for iter=1:nmcmc
      if ~noCounter; counter(iter); end

      C.aCorr=1; % default is no step correction.
      
      for varNum=1:length(svars)
       C.var=svars{varNum};
       switch(C.var)
       case {'theta'}
         for k=1:svarSize(varNum)
           C.index=k;C.val=model.(C.var)(k) + (rand(1)-.5)*mcmc.(wvars{varNum});
           model=mcmcEval(model,data,priors,C);
         end
       case {'betaV','betaU'}
         for k=1:svarSize(varNum)
           cand = exp(-model.(svars{varNum})(k).*(.5^2))+(rand(1)-.5)*mcmc.(wvars{varNum});
           C.index=k;C.val=-log(cand)/(0.5^2);
           model=mcmcEval(model,data,priors,C);
         end
       case {'lamVz','lamUz','lamWs','lamWOs','lamOs'}
         for k=1:svarSize(varNum)
           C.index=k;
		   [C.val C.aCorr]=chooseVal(model.(C.var)(k));
           if C.aCorr; model=mcmcEval(model,data,priors,C); end
         end
         otherwise  
  	       error('Unknown sample variable in gpmmcmc mcmcStep')
  	   end
	  end

      % Save the designated fields
      for var=pvars
        pvals(poff+iter).(var{1})=model.(var{1})(:);
      end

    end

    params.pvals=pvals;
    params.model=model;

    %%%%%%%%%%%%%%%%%%%%%%%%%%
    function [dval,acorr]=chooseVal(cval)
      % select the interval, and draw a new value.
      w=max(1,cval/3);
      dval=cval + (rand*2-1)*w;
      % do a correction, which depends on the old and new interval
      w1=max(1,dval/3);
      if cval > (dval+w1)
        acorr=0;
      else
        acorr=w/w1;
      end
    end
  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function mcmcStep
    
  C.aCorr=1;
      
  for iter=1:nmcmc
    if ~noCounter; counter(iter); end

    for varNum=1:length(svars)
     C.var=svars{varNum};
     switch(C.var)
     case {'theta','lamVz','lamUz','lamWs','lamWOs','lamOs'}
       for k=1:svarSize(varNum)
         C.index=k;C.val=model.(C.var)(k) + (rand(1)-.5)*mcmc.(wvars{varNum})(k);
         model=mcmcEval(model,data,priors,C);
       end
       case {'betaV','betaU'}
        for k=1:svarSize(varNum)
          cand = exp(-model.(svars{varNum})(k).*(.5^2))+(rand(1)-.5)*mcmc.(wvars{varNum})(k);
          C.index=k;C.val=-log(cand)/(0.5^2);
          model=mcmcEval(model,data,priors,C);
        end
       otherwise
	     error('Unknown sample variable in gpmmcmc mcmcStep')
	 end
	end

    % Save the designated fields
    for var=pvars
        pvals(poff+iter).(var{1})=model.(var{1})(:);
    end

  end

  params.pvals=pvals;
  params.model=model;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function mcmcStepInit
    % set up candidate acceptances tracking.
    if ~isempty(mcmc.base)
      ex = mcmc.ex; base = mcmc.base; acc = mcmc.acc;
	  ct = mcmc.ct; ct = mod(ct+1,length(ex)); mcmc.ct = ct;
    else
	  % initialize the step data structures
      ct = 0; mcmc.ct = ct;
      ex = -(nLev-1)/2:(nLev-1)/2; 
	  mcmc.ex = ex;
      for varNum=1:length(svars)
        label=svars{varNum};
        base.(label)=2*ones(nLev,svarSize(varNum));
        acc.(label)=zeros(nLev,svarSize(varNum));
		% if the width parameter is scalar, expand it into a vector
		  wlabel=wvars{varNum};
		  if length(mcmc.(wlabel))==1
		    mcmc.(wlabel)=mcmc.(wlabel)*ones(1,svarSize(varNum));
	      end
      end
   	  % specialize base matrices in some cases
       exGZ = (ex>0);
	   for varNum=1:length(svars)
	    var=svars{varNum};		
  	      switch(var)
		  case 'theta'
            base.(var)(exGZ,:)=20.0^(2.0/(nLev-1));
		  case 'betaV'   
            base.(var)(exGZ,:)=10.0^(2.0/(nLev-1));
          case 'betaU'
            base.(var)(exGZ,:)=10.0^(2.0/(nLev-1));
          case 'lamVz'
          case 'lamUz'
            base.(var)(exGZ,:)=100.0^(2.0/(nLev-1));
          case 'lamWs'
          case 'lamWOs'
		  case 'lamOs'
            base.(var)(exGZ,:)=100.0^(2.0/(nLev-1));
		  otherwise
		    error('invalid var in initialize base, gpmmcmc')
		  end
	   end
       mcmc.base = base;
    end

    nex=length(ex);
    C.aCorr=1;
    
    % Do the mcmc
    for iter=1:nmcmc
      if ~linkedMod, ct = mod(iter-1,nex); end
      jter = ex(ct+1); kter = jter+(1+nLev)/2;
      if ~noCounter; counter(iter); end

      for varNum=1:length(svars)
  	   C.var=svars{varNum};
       switch(C.var)
       case {'theta','lamVz','lamUz','lamWs','lamWOs','lamOs'}
        for k=1:svarSize(varNum)
         C.index=k;C.val=model.(C.var)(k) + (rand(1)-.5)*mcmc.(wvars{varNum})(k)*...
           (base.(C.var)(kter,k)^jter);
         model=mcmcEval(model,data,priors,C);
         acc.(C.var)(kter,k)=acc.(C.var)(kter,k)+model.acc;
        end
       case {'betaV','betaU'}
        for k=1:svarSize(varNum)
          cand = exp(-model.(svars{varNum})(k).*(.5^2))+(rand(1)-.5)*mcmc.(wvars{varNum})(k)*...
            (base.(C.var)(kter,k)^jter);
          C.index=k;C.val=-log(cand)/(0.5^2);
          model=mcmcEval(model,data,priors,C);
          acc.(C.var)(kter,k)=acc.(C.var)(kter,k)+model.acc;
        end
	   otherwise
	     error('Unknown sample variable in gpmmcmc mcmcStepInit')
	   end
	  end

      % Record fields designated
      for var=pvars
        pvals(poff+iter).(var{1})=model.(var{1})(:);
      end
    end

    params.pvals=pvals;
    params.model=model;
    mcmc.acc=acc;
    params.mcmc=mcmc;

  end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function model=mcmcEval(model,data,priors,C)

    pr=0; % print diagnostics
    if pr; fprintf('%s %2d %.2f ',C.var,C.index,C.val); end

    model.acc=0; % we'll record whether this was an accept or not

    %check whether we're in a multi-model mode, we're dealing with
    %thetas, and this model BTW has any links
    linkedVars=linkedMod & strcmp(C.var,'theta');
    if linkedVars
      % if it is a linked model, then check to see whether
      % this theta is noted drawn by another model; if so this mcmc
      % doesn't draw the variable
      if any(params.vLinks.skipList==C.index);
        if pr; fprintf('  Skipping linked var, index %d\n',C.index); end
        return
      end
    end

    % Check hard parameter bounds.
    if (C.val<priors.(C.var).bLower || ...
        priors.(C.var).bUpper<C.val || ...
        ~isreal(C.val));
      if pr; fprintf(' Reject, out of bounds\n'); end
      return
    end

    modelT=model;
    modelT.(C.var)(C.index)=C.val;
    modelT=computeLogPrior(priors,modelT);
    modelT=computeLogLik(modelT,data,C);
    modelT.logPost=modelT.logPrior+modelT.logLik;

    % if theta, consult the constraint function
    if strcmp(C.var,'theta')
      theta=modelT.theta;
      constraintsOK=1;
      for const=priors.theta.constraints
        constraintsOK=constraintsOK & eval(const{1});
      end
      if ~constraintsOK
        if pr; fprintf('  Reject, from theta constraint set\n'); end
        return
      end
    end

    % If we are in a linked model, compute the link correction
    logLikCorr=[];
    if linkedVars
      %check if the current variable is linked
      linkList=find(params.vLinks.sharedIndices==C.index);
      % Get a correction for each link
      if ~isempty(linkList) % may be null
        for kk=1:length(linkList)
          logLikCorr(kk)=...
            params.vLinks.logLikCorrHandle(...
            params.vLinks.vars(linkList(kk)).links, ...
            modelT.theta(C.index) );
        end
        if pr;
          fprintf('  Linked vars; LL change goes from %8.5f by %8.5f\n', ...
            modelT.logLik-model.logLik,sum(logLikCorr));
        end
      end
    end
    logPost=modelT.logPost+sum(logLikCorr);

    if pr; fprintf(':: %.2f %.2f ',modelT.logLik,modelT.logPrior); end

    if ( log(rand)<(logPost-model.logPost + log(C.aCorr)) )
      model=modelT;
      model.acc=1;
      % If we linked to the variable in other models, we need to push
      % out the updated value.
      if linkedVars
        for kk=1:length(linkList)
          params.vLinks.pushValHandle(...
            params.vLinks.vars(linkList(kk)).links, ...
            model.theta(C.index));
        end
      end
      if pr; fprintf(' Accept \n'); end
    else
      if pr; fprintf(' Reject \n'); end
    end

  end % nested function mcmcEval

end %main function gaspmcmc
