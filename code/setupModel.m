% function params=setupModel(obsData,simData,optParms)
% Sets up a gpmsa runnable struct from raw data.
% Please refer to associated documentation

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

function params=setupModel(obsData,simData,optParms)

  % params is optional
    if ~exist('optParms'); optParms=[]; end

  % check for scalar output
    if isfield(optParms,'scalarOutput'); scOut=optParms.scalarOutput;
                                  else scOut=0; end
    model.scOut=scOut;

  % grab some parms to local (short named) vars
    n=length(obsData);
	  if n==0;  % eta-only model ... dummy up some empty fields to satisfy setup
	    obsData(1).x=[]; 
		obsData(1).Dobs=[];
		obsData(1).yStd=[];
	  end
	m=size(simData.yStd,2);
	if n==0  % then all vars are x, otherwise some are x some are theta
	  p=size(simData.x,2);
	  q=0;
    else
      p=length(obsData(1).x);
      q=size(simData.x,2)-p;
    end
    pu=size(simData.Ksim,2);
    pv=size(obsData(1).Dobs,2);

  % check for and process lamVzGroups
    if isfield(optParms,'lamVzGroup'); lamVzGroup=optParms.lamVzGroup;
                                  else lamVzGroup=ones(pv,1); end
    lamVzGnum=length(unique(lamVzGroup));
    if ~isempty(setxor(unique(lamVzGroup),1:lamVzGnum))
      error('invalid lamVzGroup specification in setupModel');
    end

  % put in a Sigy param if not supplied (backward compatability)
    if ~isfield(obsData,'Sigy')
      for k=1:n; obsData(k).Sigy=eye(size(obsData(k).Kobs,1)); end
    end
  % make a local copy of Lamy for use in this routine (do inv() only once)
    for k=1:n; obs(k).Lamy=inv(obsData(k).Sigy); end

  % Construct the transformed obs
    if scOut
       data.x=[]; data.u=[];
       for k=1:n;
         data.x(k,:)=obsData(k).x;
         data.u(k)=obsData(k).yStd;
       end;
    else
      % ridge to be used for stabilization
       DKridge=eye(pu+pv)*1e-4;

       data.x=[]; data.v=[]; data.u=[];
       for k=1:n;
         data.x(k,:)=obsData(k).x;
        % Transform the obs data
         DK=[obsData(k).Dobs obsData(k).Kobs];
         vu=inv( DK'*obs(k).Lamy*DK + DKridge )* ...
            DK'*obs(k).Lamy*obsData(k).yStd;    
         data.v(k,:)=vu(1:pv)';
         data.u(k,:)=vu(pv+1:end)';
       end;
    end
    data.z=simData.x(:,1:p);
    data.t=simData.x(:,p+1:end);
    data.w=(simData.Ksim\simData.yStd)'; % Construct the transformed sim

  % Set initial parameter values
    model.theta=0.5*ones(1,q);           % Estimated calibration variable
    model.betaV=ones(p,lamVzGnum)*0.1;   % Spatial dependence for V discrep
    model.lamVz=ones(lamVzGnum,1)*20;    % Marginal discrepancy precision
    model.betaU=ones(p+q,pu)*0.1;        % Sim PC surface spatial dependence
    model.lamUz=ones(pu,1)*1;            % Marginal precision
    model.lamWs=ones(pu,1)*1000;         % Simulator data precision

% Set up partial results to be stored and passed around;
  % Sizes, for reference:
    model.n=n; model.m=m; model.p=p; model.q=q;
    model.pu=pu; model.pv=pv;
    model.lamVzGnum=lamVzGnum; model.lamVzGroup=lamVzGroup;
  % Precomputable data forms and covariograms.
    if isfield(data,'typeCategorical')
      dataDesc.typeCategorical=data.typeCategorical;
    else
      dataDesc.typeCategorical=zeros(1,p+q);
    end
    model.x0Dist=genDist(data.x,dataDesc);
    model.zDist=genDist([data.z data.t],dataDesc);
    model.w=data.w(:);
    if scOut
       model.uw=[data.u(:);data.w(:)];
       model.u=data.u(:);
    else
       model.vuw=[data.v(:);data.u(:);data.w(:)];
       model.vu=[data.v(:);data.u(:)];
    end

  % compute the PC loadings corrections
    model.LamSim=diag(simData.Ksim'*simData.Ksim);
    
  % initialize the acceptance record field
    model.acc=1;

  % compute LamObs, the u/v spatial correlation
    if scOut
      LO = zeros(n*pu);
      for kk=1:n
         ivals = (1:pu)+(kk-1)*pu;
         LO(ivals,ivals) = obs(kk).Lamy;
      end
      rankLO = rank(LO);
    else
      LO = zeros(n*(pv+pu));
      for kk=1:n
        DK = [obsData(kk).Dobs obsData(kk).Kobs];
        ivals = (1:pv+pu)+(kk-1)*(pv+pu);
        LO(ivals,ivals) = DK'*obs(kk).Lamy*DK;
      end
      rankLO = rank(LO);
      for kk=1:n
        ivals = (1:pv+pu)+(kk-1)*(pv+pu);
        LO(ivals,ivals) = LO(ivals,ivals) + DKridge;
      end
      % now reindex LamObs so that it has the v's first and the
      % u's 2nd.  LamObs is n*(pu+pv) in size and indexed in
      % the order v1 u1 v2 u2 ... vn un.  We want to arrange the
      % order to be v1 v2 ... vn u1 u2 ... un.  
      inew = [];
      for kk=1:pv
        inew = [inew; (kk:(pu+pv):n*(pu+pv))'];
      end
      for kk=1:pu
        inew = [inew; ((pv+kk):(pu+pv):n*(pu+pv))'];
      end
      LO = LO(inew,inew);
    end
    % compute the Penrose inverse of LO
    model.SigObs=pinv(LO)+1e-6*eye(size(LO,1));

  % Set prior distribution values
    priors.lamVz.a=1;  priors.lamVz.b=0.0010;
    priors.lamUz.a=5;  priors.lamUz.b=5;
    priors.lamWOs.a=5; priors.lamWOs.b=0.005;
    priors.lamWs.a=3;  priors.lamWs.b=0.003;
    priors.lamOs.a=1;  priors.lamOs.b=.001;
    priors.rhoU.b=0.1;
    priors.rhoV.b=0.1;
    priors.theta.mean=0.5; priors.theta.std=10;

    if isfield(optParms,'priors')
       if isfield(optParms.priors,'lamWOs')
          priors.lamWOs.a=optParms.priors.lamWOs.a;
          priors.lamWOs.b=optParms.priors.lamWOs.b;
       end
       if isfield(optParms.priors,'lamOs')
          priors.lamOs.a=optParms.priors.lamOs.a;
          priors.lamOs.b=optParms.priors.lamOs.b;
       end
    end

  % Prior correction for lamOs and lamWOs prior values (due to D,K basis xform)
    %for lamOs, need DK basis correction
    totElements=0; 
    for ii=1:length(obsData); 
      totElements=totElements+length(obsData(ii).yStd);
    end
    aCorr=0.5*(totElements-rankLO);

    bCorr=0;
    if ~scOut
       for ii=1:n
         DKii = [obsData(ii).Dobs   obsData(ii).Kobs];
         vuii = [data.v(ii,:)'; data.u(ii,:)'];
         resid=obsData(ii).yStd(:) - DKii*vuii;  
         bCorr=bCorr+0.5*sum(resid'*obs(ii).Lamy*resid);
       end
    end
    priors.lamOs.a=priors.lamOs.a+aCorr;
    priors.lamOs.b=priors.lamOs.b+bCorr;
    
    %for lamWOs, need K basis correction
    aCorr=0.5*(size(simData.yStd,1)-pu)*m;
    ysimStdHat = simData.Ksim*data.w';
    bCorr=0.5*sum(sum((simData.yStd-ysimStdHat).^2));

    priors.lamWOs.a=priors.lamWOs.a+aCorr;
    priors.lamWOs.b=priors.lamWOs.b+bCorr;

  % Set the initial values of lamOs and lamWOs based on the priors.
    model.lamWOs=max(100,priors.lamWOs.a/priors.lamWOs.b);
    model.lamOs=max(20, priors.lamOs.a/priors.lamOs.b);

  % Set prior bounds 
    priors.lamVz.bLower=0;      priors.lamVz.bUpper=Inf;
    priors.lamUz.bLower=0.3;    priors.lamUz.bUpper=Inf;
    priors.lamWs.bLower=60;     priors.lamWs.bUpper=1e5;
    priors.lamWOs.bLower=60;    priors.lamWOs.bUpper=1e5;
    priors.lamOs.bLower=0;      priors.lamOs.bUpper=Inf;
    priors.betaU.bLower=0;      priors.betaU.bUpper=Inf;
    priors.betaV.bLower=0;      priors.betaV.bUpper=Inf;
    priors.theta.bLower=0;      priors.theta.bUpper=1;
    % if thetaConstraintFunction supplied, use that, otherwise
    % use a dummy constraint function
      if isfield(optParms,'thetaConstraints')
        priors.theta.constraints=optParms.thetaConstraints;
        % select a valid starting theta
        ii=0;
        while (ii<1e6) & ~tryConstraints(priors.theta.constraints,model.theta)
          model.theta=rand(size(model.theta));
        end
        if ii==1e6; error('unable to draw theta within constraints'); end
      else
        priors.theta.constraints={};
      end

      function constraintsOK=tryConstraints(constraints,theta)
        constraintsOK=1;
        for const=priors.theta.constraints
          constraintsOK=constraintsOK & eval(const{1}); 
        end
      end
      
  % Set mcmc step interval values
    mcmc.thetawidth=0.2;
    mcmc.rhoUwidth=0.1;
    mcmc.rhoVwidth=0.1;
    mcmc.lamVzwidth=10;
    mcmc.lamUzwidth=5;
    mcmc.lamWswidth=100;
    mcmc.lamWOswidth=100;
    mcmc.lamOswidth=model.lamOs/2;
  % set up control var lists for sampling and logging
    % pvars is the list variables to log
    % svars is the list variables to log
    % svarSize is the length of each svar variable
    % wvars is the list of corresponding mcmc width names
    if n>0 % if there's obsData, do the full deal.
	  mcmc.pvars={'theta','betaV','betaU','lamVz','lamUz','lamWs', ...
                  'lamWOs','lamOs','logLik','logPrior','logPost'};
	  mcmc.svars={'theta','betaV','betaU','lamVz', ...
	              'lamUz','lamWs','lamWOs','lamOs'};
   	  mcmc.svarSize=[q              % theta
	                 p*lamVzGnum    % betaV
				     pu*(p+q)       % betaU
				  	 lamVzGnum      % lamVz
					 pu             % lamUz
					 pu             % lamWs
					 1              % lamWOs
					 1]';           % lamOs
	  mcmc.wvars={'thetawidth','rhoVwidth','rhoUwidth','lamVzwidth', ...
	              'lamUzwidth','lamWswidth','lamWOswidth','lamOswidth'};
    else % we're doing just an eta model, so a subset of the params. 
	  mcmc.pvars={'theta','betaU','lamUz','lamWs', ...
                      'lamWOs','logLik','logPrior','logPost'};
	  mcmc.svars={'theta','betaU','lamUz','lamWs','lamWOs'};
   	  mcmc.svarSize=[q              % theta
				     pu*(p+q)       % betaU
					 pu             % lamUz
					 pu             % lamWs
					 1]';           % lamWOs
	  mcmc.wvars={'thetawidth','rhoUwidth', ...
	              'lamUzwidth','lamWswidth','lamWOswidth'};
    end				  

% Over and out
  params.data=data;
  params.model=model;
  params.priors=priors;
  params.mcmc=mcmc;
  params.obsData=obsData;
  params.simData =simData;
  params.pvals=[];  % initialize the struct

end
