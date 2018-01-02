%function out=optimizeTheta(pout,varargin)
% varargs are
%   mode - 'optimize' or 'characterize'
%   theta0 - starting point for optimization; marginal median if not given
%   nPvals - number of points to average over, default is 10
%   errType- type of error to optimize, 'RMSE' (default), 'Profile' or
%            'Posterior'
%   theta0Mode - if this is set to 'random', the choose starting thetas from
%               tRange
%   tRange - random theta starting point range, default [0 1];
%   rmseOption - use posterior mean, median or user specified parameter values
%                for RMSE optimization, 'mean', 'median' or 'user'
%   rmsePvals - parameter values associated with rmseOption = 'user'
%   clist - indicates the index of each variable in each model (0 if
%           variable is not in model)
  
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

function out=optimizeTheta(pout,varargin)

  % parse the arguments
    theta0=[];
    nPvals=10;
    errType='RMSE';
    verbose=false;
    mode='optimize';
    optType='simplex';
    optRuns=1;
    tRange=[];
    tNoise=[];
    theta0Mode='specified';
    rmseOption='';
    rmsePvals=[];
    clist=[];
    parseAssignVarargs({'mode','theta0','nPvals','errType','verbose', ...
                        'tRange','tNoise','optRuns','theta0Mode', ...
                        'optType', 'rmseOption', 'rmsePvals', 'clist'});

  if clist, n=size(clist,1); else n=pout(1).model.q; end

  thetaMed=zeros(n,1);
  for ii=1:length(pout)
     P(ii)=pout(ii).model.p; Q(ii)=pout(ii).model.q;
     PU(ii)=pout(ii).model.pu; PV(ii)=pout(ii).model.lamVzGnum;
     nBurn(ii)=pout(ii).nburn; nLev(ii)=pout(ii).nlev;
     p=P(ii);q=Q(ii);pu=PU(ii);pv=PV(ii);nburn=nBurn(ii);nlev=nLev(ii);
     pout(ii).model.betaV=median([pout(ii).pvals(nburn*nlev+1:end).betaV]');
     pout(ii).model.betaV=reshape(pout(ii).model.betaV,[],pv);
     pout(ii).model.betaU=median([pout(ii).pvals(nburn*nlev+1:end).betaU]');
     pout(ii).model.betaU=reshape(pout(ii).model.betaU,[],pu);
     pout(ii).model.lamVz=median([pout(ii).pvals(nburn*nlev+1:end).lamVz]');
     pout(ii).model.lamUz=median([pout(ii).pvals(nburn*nlev+1:end).lamUz]');
     pout(ii).model.lamWs=median([pout(ii).pvals(nburn*nlev+1:end).lamWs]');
     pout(ii).model.lamWOs=median([pout(ii).pvals(nburn*nlev+1:end).lamWOs]');
     pout(ii).model.lamOs=median([pout(ii).pvals(nburn*nlev+1:end).lamOs]');
     if clist, jj=find(clist(:,ii) ~= 0); else jj=1:n; end
     thetaMed(jj)=median([pout(ii).pvals(nburn*nlev+1:end).theta]');
  end

  iss=[];
  switch rmseOption
    case 'mean'
      nPvals=1;
      for ii=1:length(pout)
         pu=PU(ii);pv=PV(ii);nburn=nBurn(ii);nlev=nLev(ii);
         pvals(ii).theta=mean([pout(ii).pvals(nburn*nlev+1:end).theta]');
         pvals(ii).betaV=mean([pout(ii).pvals(nburn*nlev+1:end).betaV]');
         pvals(ii).betaV=reshape(pvals(ii).betaV,[],pv);
         pvals(ii).betaU=mean([pout(ii).pvals(nburn*nlev+1:end).betaU]');
         pvals(ii).betaU=reshape(pvals(ii).betaU,[],pu);
         pvals(ii).lamVz=mean([pout(ii).pvals(nburn*nlev+1:end).lamVz]');
         pvals(ii).lamUz=mean([pout(ii).pvals(nburn*nlev+1:end).lamUz]');
         pvals(ii).lamWs=mean([pout(ii).pvals(nburn*nlev+1:end).lamWs]')';
         pvals(ii).lamWOs=mean([pout(ii).pvals(nburn*nlev+1:end).lamWOs]');
         pvals(ii).lamOs=mean([pout(ii).pvals(nburn*nlev+1:end).lamOs]');
      end
    case 'median'
      nPvals=1;
      for ii=1:length(pout)
         pvals(ii).theta=thetaMed; pvals(ii).betaV=pout(ii).model.betaV;
         pvals(ii).betaU=pout(ii).model.betaU; pvals(ii).lamVz=pout(ii).model.lamVz;
         pvals(ii).lamUz=pout(ii).model.lamUz; pvals(ii).lamWs=pout(ii).model.lamWs';
         pvals(ii).lamWOs=pout(ii).model.lamWOs; pvals(ii).lamOs=pout(ii).model.lamOs;
      end
    case 'user'
      nPvals=1;
      pvals=rmsePvals;
    otherwise
      iss=zeros(length(pout),nPvals);
      for ii=1:length(pout)
         nburn=nBurn(ii);nlev=nLev(ii);
         iss(ii,:)=int32(linspace(nburn*nlev+1,length(pout(ii).pvals),nPvals));
      end
  end

  %fprintf('%5.4f ',theta0); fprintf('\n');
  
  switch mode
    case 'optimize'
      if optRuns>1; counter('stime',1,optRuns,5,10); end
      for ii=1:optRuns;
        if optRuns>1; counter(ii); end
        if ~isempty(tRange)
           theta0=rand(1,n) * (tRange(2)-tRange(1)) + tRange(1);
        elseif ~isempty(tNoise)
           theta0=thetaMed+(rand(size(thetaMed))-0.5)*tNoise;
           theta0(theta0<0)=0;
           theta0(theta0>1)=1;
        elseif isempty(theta0)
           theta0=thetaMed;
        end
        out(ii)=optMode;
         %out(ii).theta0=theta0;
         %out(ii).thetaMed=thetaMed;
      end
      if optRuns>1; counter('end'); end
    case 'characterize'
      out=charMode;
    otherwise error('invalid mode specified in optimizeTheta');
  end

  function out=optMode

    options=[];
    %options.Display='iter';
    switch optType
      case 'simplex'
        theta=fminsearch(@calcErr,theta0,options);
      case 'hillClimb'
        theta=theta0;
        err=calcErr(theta0); errLast=Inf;
        nsla=0;
        while abs(err-errLast) > (1e-4)
          thetaNew=theta+(rand(size(theta))-0.5)* (0.1/log((nsla+2)));
          thetaNew(thetaNew<0)=0; thetaNew(thetaNew>1)=1;
          errNew=calcErr(thetaNew);
             %fprintf('err=%f try=%f\n',err,errNew);
             %if errNew>0; keyboard; end
          if errNew<err
             theta=thetaNew;
             errLast=err;
             err=errNew;
             nsla
             nsla=sqrt(nsla);
             fprintf('Accept err=%f last=%f\n',err,errLast);
          else
            nsla=nsla+1;
          end
        end
    end
    out.theta=theta;
    out.err=calcErr(theta);
    out.theta0=theta0;
    out.t0err=calcErr(theta0);
    
    if verbose
      R0=  calcErr(theta0,'k');
      Rout=calcErr(out.theta,'r');
      subplot(2,1,1); 
      title(['t0=[' sprintf('%4.3f ',theta0) ']; tout=[' sprintf('%4.3f ',out.theta) ']']); 
      subplot(2,1,2); 
      title(sprintf('R0=%f, Rout=%f',R0,Rout));
    end
  end

  function infs=charMode
    % pace over various line lengths and find the number of inflections
    % along each, in the error function
    nsamps=10; ntrials=50;
    len=0.05:0.05:0.5;
    counter('stime',1,length(len)*ntrials,10,10);
    for ilen=1:length(len)
      for ii=1:ntrials;
      counter(sub2ind([ntrials length(len)],ii,ilen));
      % choose endpoints
        p1=rand(1,n);
        p2=-1;
        while any(p2<0) | any(p2>1) % must choose a second point in the range
          offset=rand(1,n) - 0.5; % get an offset vector
          offset=len(ilen)*offset/sqrt(sum(offset.^2)); % and normalize it
          p2=p1+offset;
        end
      % make the line points
        for jj=1:n; 
          spoint(:,jj)=linspace(p1(jj),p2(jj),nsamps); 
        end
      % get the response for all the points
        for jj=1:nsamps
          err(jj)=calcErr(spoint(jj,:));
        end
      % calculate and record the number of inflections
      % (first diff is up/down move direction, second is infection present)
        udmove=sign(diff(err));
        infs(ilen,ii)=sum(diff(udmove)~=0);
      end
    end
    counter('end');
  end

  function err=calcErr(theta,plotColor)
    switch(errType)
      case 'RMSE'
        if any(theta<0) | any(theta>1); err=realmax; return; end
        yRMSE=zeros(nPvals,1); err=0;
        for ii=1:length(pout)
           if clist, jj=find(clist(:,ii) ~= 0); else jj=1:n; end
           if iss,
              pred=gPred(0.5,pout(ii).pvals(iss(ii,:)),pout(ii).model,...
                         pout(ii).data,'wpred',theta(jj));
           else
              pred=gPred(0.5,pvals(ii),pout(ii).model,pout(ii).data,...
                         'wpred',theta(jj));
           end
      
           yobsStd=pout(ii).obsData.yStd;
           yobs=pout(ii).obsData.orig.y;
           ypredStd=pred.Myhat*pout(ii).obsData.Kobs';
           ypred=ypredStd*pout(ii).simData.orig.ysd + ...
                 repmat(pout(ii).obsData.orig.ymean',size(ypredStd,1),1);
  
           yRMSE=yRMSE+sum( (ypredStd-repmat(yobsStd',size(ypred,1),1)).^2 ,2);
           err=err+size(pout(ii).obsData.Kobs,1);
    
           if exist('plotColor')
             figure(ii);
             subplot(2,1,1);
             plot(mean(ypredStd),plotColor); hold on; plot(yobsStd,'*');
             subplot(2,1,2);
             plot(mean(ypred),plotColor); hold on; plot(yobs,'*');
           end    
        end
        err=mean(sqrt(yRMSE./err));
      case 'Profile'
        if any(theta<0) | any(theta>1); err=realmax; return; end
        C.var='all';
        err=0;
        for ii=1:length(pout)
           model=pout(ii).model;
           if clist, jj=find(clist(:,ii) ~= 0); else jj=1:n; end
           model.theta=theta(jj);
           model.logLik=0; model.logPrior=0; model.logPost=0;
           modll=computeLogLik(model,pout(ii).data,C);
           modpr=computeLogPrior(pout(ii).priors,model);
           err=err-(modll.logLik+modpr.logPrior);
        end
        %fprintf('%20.10f\n',err);
        %keyboard
      case 'Posterior'
        if any(theta<0) | any(theta(1:n)>1); err=realmax; return; end
	C.var='all';
        err=0; skip=0;
        for ii=1:length(pout)
	   model=pout(ii).model;
           p=P(ii);q=Q(ii);pu=PU(ii);pv=PV(ii);
           if clist, jj=find(clist(:,ii) ~= 0); else jj=1:n; end
	   model.theta=theta(jj);
           model.betaV=reshape(theta(skip+n+1:skip+n+p*pv),[],pv);
           model.betaU=reshape(theta(skip+n+p*pv+1:skip+n+q*pu+p*(pu+pv)),[],pu);
           model.lamVz=theta(skip+n+q*pu+p*(pu+pv)+1:skip+n+q*pu+p*pu+pv*(p+1));
           model.lamUz=theta(skip+n+q*pu+p*pu+pv*(p+1)+1:skip+n+q*pu+(pu+pv)*(p+1));
           model.lamWs=theta(skip+n+q*pu+(pu+pv)*(p+1)+1:skip+n+pu*(q+1)+(pu+pv)*(p+1));
           model.lamWOs=theta(skip+n+pu*(q+1)+(pu+pv)*(p+1)+1);
           model.lamOs=theta(skip+n+pu*(q+1)+(pu+pv)*(p+1)+2);
	   model.logLik=0; model.logPrior=0; model.logPost=0;
	   modll=computeLogLik(model,pout(ii).data,C);
	   modpr=computeLogPrior(pout(ii).priors,model);
	   err=err-(modll.logLik+modpr.logPrior);
           skip=skip+pu*(q+1)+(pu+pv)*(p+1)+2;
        end
	%fprintf('%20.10f\n',err);
	%keyboard
    end
  end
end
