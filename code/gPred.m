%function pred=gPred(xpred,pvals,model,data,mode,theta)
%  Predict using a gpmsa constructed model. 
%  please see associated documentation

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

function pred=gPred(xpred,pvals,model,data,mode,theta)

if isfield(data,'typeCategorical')
  dataDesc.typeCategorical=data.typeCategorical;
else
  dataDesc.typeCategorical=zeros(1,model.p+model.q);
end

if exist('theta','var'); ncal=1; else ncal=0; end

switch mode
  case 'uvpred'
    if ncal
      pred=uvPred(xpred,pvals,model,data,dataDesc,theta);
    else
      pred=uvPred(xpred,pvals,model,data,dataDesc);
    end
  case 'wpred'
    if ncal
      pred=wPred(xpred,pvals,model,data,dataDesc,theta);
    else
      pred=wPred(xpred,pvals,model,data,dataDesc);
    end
  case 'etamod'
  	pred=wEtaModelPred(xpred,pvals,model,data,dataDesc);
  otherwise
    error('invalid model in gPred');
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pred=wEtaModelPred(xpred,pvals,model,data,dataDesc)

m=model.m;p=model.p;q=model.q;pu=model.pu;

npred=size(xpred,1);

diags1=diagInds(m*pu);
diags2=diagInds(npred*pu);

nreal=length(pvals);
tpred=zeros([nreal,npred*pu]);

for ii=1:length(pvals)
  theta=pvals(ii).theta'; 
  betaU=reshape(pvals(ii).betaU,p+q,pu);
  lamUz=pvals(ii).lamUz;
  lamWs=pvals(ii).lamWs; lamWOs=pvals(ii).lamWOs;

  xpredDist=genDist(xpred,dataDesc);  
  zxpredDist=genDist2([data.z data.t],xpred,dataDesc);

    SigW=zeros(m*pu);
      for jj=1:pu
        bStart=(jj-1)*m+1; bEnd=bStart+m-1; 
        SigW(bStart:bEnd,bStart:bEnd)=...
            covarioGram(model.zDist,betaU(:,jj),lamUz(jj));
      end
      SigW(diags1)=SigW(diags1)+ ...
          kron(1./(model.LamSim*lamWOs)',ones(1,m)) + ...
          kron(1./(lamWs)',ones(1,m)) ;

    SigWp=zeros(npred*pu);
      for jj=1:pu
        bStart=(jj-1)*npred+1; bEnd=bStart+npred-1;
        SigWp(bStart:bEnd,bStart:bEnd)= ...
            covarioGram(xpredDist,betaU(:,jj),lamUz(jj));
      end
      SigWp(diags2)=SigWp(diags2)+ ...
           kron(1./(lamWs)',ones(1,npred)) ;
           % kron(1./(model.LamSim*lamWOs)',ones(1,npred)) + ...
  
    SigWWp=zeros(m*pu,npred*pu);
      for jj=1:pu
        bStartI=(jj-1)*m+1; bEndI=bStartI+m-1;
        bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
        SigWWp(bStartI:bEndI,bStartJ:bEndJ)=...
            covarioGram(zxpredDist,betaU(:,jj),lamUz(jj));
      end
      
  SigData=SigW;
  SigPred=SigWp;
  SigCross=SigWWp;
  % Get the stats for the prediction stuff. 
    W=(SigCross')/SigData;
    Myhat=W*(data.w(:));
    Syhat=SigPred-W*SigCross;
    % And do a realization
    tpred(ii,:)=rmultnorm(1,Myhat,Syhat)';

   % add the distribution params
    pred.Myhat(ii,:)=Myhat;
    pred.Syhat{ii}=Syhat;
end

% Reshape the pred matrix to 3D:
%  first dim  - (number of realizations [pvals])
%  second dim - (number of principal components)
%  third dim  - (number of points [x,theta]s) 
  pred.w=zeros(length(pvals),pu,npred);
  for ii=1:pu
    pred.w(:,ii,:)=tpred(:,(ii-1)*npred+1:ii*npred);
  end

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pred=wPred(xpred,pvals,model,data,dataDesc,thetapred)

if exist('thetapred','var')
   ncal=1;
   if any(thetapred<0)
     randcal=1;
   else
     randcal=0;
   end
 else 
   ncal=0;
   randcal=0;
end

m=model.m;p=model.p;q=model.q;pu=model.pu;

npred=size(xpred,1);

diags1=diagInds(m*pu);
diags2=diagInds(npred*pu);

nreal=length(pvals);
tpred=zeros([nreal,npred*pu]);

for ii=1:length(pvals)
  theta=pvals(ii).theta'; 
  betaU=reshape(pvals(ii).betaU,p+q,pu);
  lamUz=pvals(ii).lamUz;
  lamWs=pvals(ii).lamWs; lamWOs=pvals(ii).lamWOs;

  if ncal
    if randcal
      xpredt=[xpred rand(size(thetapred))];
    else
      xpredt=[xpred thetapred];
    end
  else
    xpredt=[xpred repmat(theta,npred,1)];
  end
  
  xpredDist=genDist(xpredt,dataDesc);  
  zxpredDist=genDist2([data.z data.t],xpredt,dataDesc);

    SigW=zeros(m*pu);
      for jj=1:pu
        bStart=(jj-1)*m+1; bEnd=bStart+m-1; 
        SigW(bStart:bEnd,bStart:bEnd)=...
            covarioGram(model.zDist,betaU(:,jj),lamUz(jj));
      end
      SigW(diags1)=SigW(diags1)+ ...
          kron(1./(model.LamSim*lamWOs)',ones(1,m)) + ...
          kron(1./(lamWs)',ones(1,m)) ;

    SigWp=zeros(npred*pu);
      for jj=1:pu
        bStart=(jj-1)*npred+1; bEnd=bStart+npred-1;
        SigWp(bStart:bEnd,bStart:bEnd)= ...
            covarioGram(xpredDist,betaU(:,jj),lamUz(jj));
      end
      SigWp(diags2)=SigWp(diags2)+ ...
           kron(1./(lamWs)',ones(1,npred)) ;
           % kron(1./(model.LamSim*lamWOs)',ones(1,npred)) + ...
  
    SigWWp=zeros(m*pu,npred*pu);
      for jj=1:pu
        bStartI=(jj-1)*m+1; bEndI=bStartI+m-1;
        bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
        SigWWp(bStartI:bEndI,bStartJ:bEndJ)=...
            covarioGram(zxpredDist,betaU(:,jj),lamUz(jj));
      end
      
  SigData=SigW;
  SigPred=SigWp;
  SigCross=SigWWp;
  % Get the stats for the prediction stuff. 
    W=(SigCross')/SigData;
    Myhat=W*(data.w(:));
    Syhat=SigPred-W*SigCross;
    % And do a realization
    tpred(ii,:)=rmultnorm(1,Myhat,Syhat)';

   % add the distribution params
    pred.Myhat(ii,:)=Myhat;
    pred.Syhat{ii}=Syhat;

end

% Reshape the pred matrix to 3D:
%  first dim  - (number of realizations [pvals])
%  second dim - (number of principal components)
%  third dim  - (number of points [x,theta]s) 
  pred.w=zeros(length(pvals),pu,npred);
  for ii=1:pu
    pred.w(:,ii,:)=tpred(:,(ii-1)*npred+1:ii*npred);
  end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function pred=uvPred(xpred,pvals,model,data,dataDesc,thetapred)
n=model.n;m=model.m;p=model.p;q=model.q;pu=model.pu;pv=model.pv;
lamVzGnum=model.lamVzGnum; lamVzGroup=model.lamVzGroup;

if exist('thetapred','var')
   ncal=1;
   if any(thetapred<0)
     randcal=1;
   else
     randcal=0;
   end
 else 
   ncal=0;
   randcal=0;
end

npred=size(xpred,1);

diags0=diagInds(n*pu);
diags1=diagInds(m*pu);
diags2=diagInds(npred*pu);

x0Dist=genDist(data.x,dataDesc);
xpred0Dist=genDist(xpred,dataDesc);
xxpred0Dist=genDist2(data.x,xpred,dataDesc);

nreal=length(pvals);
tpred=zeros([nreal,npred*(pv+pu)]);

for ii=1:length(pvals)
  theta=pvals(ii).theta';
  betaV=reshape(pvals(ii).betaV,p,lamVzGnum); 
  betaU=reshape(pvals(ii).betaU,p+q,pu);
  lamVz=pvals(ii).lamVz; lamUz=pvals(ii).lamUz; lamWOs=pvals(ii).lamWOs;
  lamWs=pvals(ii).lamWs; lamOs=pvals(ii).lamOs;

  if ncal
    if randcal
      xpredt=[xpred rand(size(thetapred))];
    else
      xpredt=[xpred thetapred];
    end
  else
    xpredt=[xpred repmat(theta,npred,1)];
  end
  
  xDist=genDist([data.x repmat(theta,n,1)],dataDesc);
  zDist=genDist([data.z data.t],dataDesc);
  xzDist=genDist2([data.x repmat(theta,n,1)],[data.z data.t],dataDesc);
  xpredDist=genDist(xpredt,dataDesc);  
  xxpredDist=genDist2([data.x repmat(theta,n,1)],xpredt,dataDesc);
  zxpredDist=genDist2([data.z data.t],xpredt,dataDesc);
  
  % Generate the part of the matrix related to the data
  % Four parts to compute: Sig_v, Sig_u, Sig_w, and the Sig_uw crossterm
    SigV=zeros(n*pv);
      for jj=1:lamVzGnum;
        vCov(jj).mat=covarioGram(x0Dist, betaV(:,jj), lamVz(jj));
      end
      for jj=1:pv
        bStart=(jj-1)*n+1; bEnd=bStart+n-1;
        SigV(bStart:bEnd,bStart:bEnd)=vCov(lamVzGroup(jj)).mat;
      end
    SigU=zeros(n*pu);
      for jj=1:pu
        bStart=(jj-1)*n+1; bEnd=bStart+n-1;
        SigU(bStart:bEnd,bStart:bEnd)= ...
            covarioGram(xDist,betaU(:,jj),lamUz(jj));
      end
      SigU(diags0)=SigU(diags0)+...
           kron(1./(lamWs)',ones(1,n)) ;
    SigW=zeros(m*pu);
      for jj=1:pu
        bStart=(jj-1)*m+1; bEnd=bStart+m-1; 
        SigW(bStart:bEnd,bStart:bEnd)=...
            covarioGram(zDist,betaU(:,jj),lamUz(jj));
      end
      SigW(diags1)=SigW(diags1)+ ...
          kron(1./(model.LamSim*lamWOs)',ones(1,m)) + ...
          kron(1./(lamWs)',ones(1,m)) ;
    SigUW=zeros(n*pu,m*pu);
      for jj=1:pu
        bStartI=(jj-1)*n+1; bEndI=bStartI+n-1;
        bStartJ=(jj-1)*m+1; bEndJ=bStartJ+m-1;
        SigUW(bStartI:bEndI,bStartJ:bEndJ)=...
            covarioGram(xzDist,betaU(:,jj),lamUz(jj));
      end
      if model.scOut
        SigData=[ SigU+SigV    SigUW; ...
                  SigUW'       SigW ];
        SigData(1:n*pu,1:n*pu) = ...
          SigData(1:n*pu,1:n*pu) + model.SigObs*1/lamOs;
      else
        SigData=[SigV                 zeros(n*pv,(n+m)*pu);  ...
                 zeros((n+m)*pu,n*pv) [ SigU    SigUW; ...
                                        SigUW'  SigW  ] ];
        SigData(1:n*(pv+pu),1:n*(pv+pu)) = ...
          SigData(1:n*(pv+pu),1:n*(pv+pu)) + model.SigObs*1/lamOs;
      end

  % Generate the part of the matrix related to the predictors
  % Parts to compute: Sig_vpred, Sig_upred
    SigVp=zeros(npred*pv);
      for jj=1:lamVzGnum;
        vpCov(jj).mat=covarioGram(xpred0Dist, betaV(:,jj), lamVz(jj));
      end
      for jj=1:pv
        bStart=(jj-1)*npred+1; bEnd=bStart+npred-1;
        SigVp(bStart:bEnd,bStart:bEnd)=vpCov(lamVzGroup(jj)).mat;
      end
      %SigVp(diagInds(npred*pv))=SigVp(diagInds(npred*pv))+1;
    SigUp=zeros(npred*pu);
      for jj=1:pu
        bStart=(jj-1)*npred+1; bEnd=bStart+npred-1;
        SigUp(bStart:bEnd,bStart:bEnd)= ...
            covarioGram(xpredDist,betaU(:,jj),lamUz(jj));
      end
      SigUp(diags2)=SigUp(diags2)+...
           kron(1./(lamWs)',ones(1,npred)) ;
    SigPred=[SigVp                     zeros(npred*pv,npred*pu);  ...
             zeros(npred*pu,npred*pv)  SigUp  ];

  % Now the cross-terms. 
    SigVVx=zeros(n*pv,npred*pv);    
      for jj=1:lamVzGnum;
        vvCov(jj).mat=covarioGram(xxpred0Dist, betaV(:,jj), lamVz(jj));
      end
      for jj=1:pv
        bStartI=(jj-1)*n+1; bEndI=bStartI+n-1;
        bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
        SigVVx(bStartI:bEndI,bStartJ:bEndJ)=vvCov(lamVzGroup(jj)).mat; 
      end
    SigUUx=zeros(n*pu,npred*pu);
      for jj=1:pu
        bStartI=(jj-1)*n+1; bEndI=bStartI+n-1;
        bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
        SigUUx(bStartI:bEndI,bStartJ:bEndJ)=...
            covarioGram(xxpredDist,betaU(:,jj),lamUz(jj));
      end
    SigWUx=zeros(m*pu,npred*pu);
      for jj=1:pu
        bStartI=(jj-1)*m+1; bEndI=bStartI+m-1;
        bStartJ=(jj-1)*npred+1; bEndJ=bStartJ+npred-1;
        SigWUx(bStartI:bEndI,bStartJ:bEndJ)=...
            covarioGram(zxpredDist,betaU(:,jj),lamUz(jj));
      end
    if model.scOut
      SigCross=[SigVVx                 SigUUx; ...
                zeros(m*pu,npred*pv)   SigWUx];
    else
      SigCross=[SigVVx                 zeros(n*pv,npred*pu); ...
                zeros(n*pu,npred*pv)   SigUUx; ...
                zeros(m*pu,npred*pv)   SigWUx];
    end

  if 0
    figure(3)
    subplot(2,2,1); imagesc(gScale(SigData,'sqrt'))
    subplot(2,2,2); imagesc(gScale(SigCross,'sqrt'))
    subplot(2,2,3); imagesc(gScale(SigCross','sqrt'))
    subplot(2,2,4); imagesc(gScale(SigPred,'sqrt'))
    keyboard
  end
  
  % Get the stats for the prediction stuff. 
    W=(SigCross')/SigData;
    if model.scOut, Myhat=W*model.uw; else Myhat=W*model.vuw; end
    Syhat=SigPred-W*SigCross;
    % And do a realization
    tpred(ii,:)=rmultnorm(1,Myhat,Syhat)';
  % log the distribution params
    pred.Myhat(ii,:)=Myhat;
    pred.Syhat{ii}=Syhat;
end

% Reshape the pred matrix to 3D, for each component:
%  first dim  - (number of realizations [pvals])
%  second dim - (number of principal components)
%  third dim  - (number of points [x,theta]s) 

  pred.v=zeros(length(pvals),pv,npred);
  pred.u=zeros(length(pvals),pu,npred);
  for ii=1:pv
    pred.v(:,ii,:)=tpred(:,(ii-1)*npred+1:ii*npred);
  end
  for ii=1:pu
    pred.u(:,ii,:)=tpred(:,pv*npred+((ii-1)*npred+1:ii*npred) );
  end

end


%
% Helper function rmultnorm computes multivariate normal realizations
function rnorm = rmultnorm(n,mu,cov)
  R=chol(cov);
  rnorm = repmat(mu,1,n) + R' * randn(size(mu,1),n);
end 
