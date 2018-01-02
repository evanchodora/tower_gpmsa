% function model = computeLogLik(model,data,C)
% 
% Builds the log likelihood of the data given the parameters.
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

function model = computeLogLik(model,data,C)

 n=model.n;   m=model.m;
pu=model.pu; pv=model.pv;
 p=model.p;   q=model.q;
lamVzGnum=model.lamVzGnum; lamVzGroup=model.lamVzGroup;

% validate and process the changed field
do.theta=0;do.betaV=0;do.lamVz=0;do.betaU=0;do.lamUz=0;
do.lamWs=0;do.lamWOs=0;

if strcmp(C.var,'theta') || strcmp(C.var,'all') % update distances
       if isfield(data,'typeCategorical')
         dataDesc.typeCategorical=data.typeCategorical;
       else
         dataDesc.typeCategorical=zeros(1,p+q);
       end
       model.xDist=genDist([data.x repmat(model.theta,n,1)],dataDesc);
       model.xzDist=genDist2([data.x repmat(model.theta,n,1)], ...
                              [data.z data.t],dataDesc);
end

switch C.var
    case 'all'
       do.theta=1; do.betaV=1; do.lamVz=1; do.betaU=1; 
       do.lamUz=1; do.lamWs=1; do.lamWOs=1; do.lamOs=1; 
       model.SigWl=zeros(pu,1);       
    case 'theta'; do.theta=1;
    case 'betaV'; do.betaV=1;
    case 'lamVz'; do.lamVz=1;
    case 'betaU'; do.betaU=1;
    case 'lamUz'; do.lamUz=1;
    case 'lamWs'; do.lamWs=1;
    case 'lamWOs'; do.lamWOs=1;
    case 'lamOs'; do.lamOs=1;
    otherwise
        error('Invalid Subtype in computeLogLik');
end

betaV=model.betaV;   lamVz=model.lamVz;
betaU=model.betaU;   lamUz=model.lamUz;
lamWs=model.lamWs;   lamWOs=model.lamWOs; lamOs=model.lamOs;

% Four parts to compute: Sig_v, Sig_u, Sig_w, and the Sig_uw crossterm

  if (do.theta || do.betaV || do.lamVz)
	SigV=[];
    for jj=1:lamVzGnum
      SigV(jj).mat=covarioGram(model.x0Dist, betaV(:,jj), lamVz(jj));
    end
    model.SigV=SigV;
  else
    SigV=model.SigV;
  end

  if (do.theta || do.betaU || do.lamUz || do.lamWs)
    SigU(pu).mat=[];
    diags1=diagInds(n);
    for jj=1:pu
        SigU(jj).mat=covarioGram(model.xDist,betaU(:,jj),lamUz(jj));
        SigU(jj).mat(diags1)=SigU(jj).mat(diags1)+1/lamWs(jj);
    end
    model.SigU=SigU;
  else
    SigU=model.SigU;
  end

  if (do.betaU || do.lamUz || do.lamWs || do.lamWOs)
    diags2=diagInds(m);
    switch C.var
      case 'all'
        jinds=1:pu;
      case 'betaU'
        jinds=ceil( C.index/(p+q) );
      case {'lamUz','lamWs'}
        jinds=C.index;
      case 'lamWOs'
        jinds=1:pu;
    end
    for jj=jinds
      cg=covarioGram(model.zDist,betaU(:,jj),lamUz(jj));
      cg(diags2)=cg(diags2)+1/(model.LamSim(jj)*lamWOs) + 1/lamWs(jj);

      % calculate the SigW likelihood for each block
        model.SigWl(jj)=doLogLik(cg,model.w((jj-1)*m+1:jj*m));
      % calculate the SigW inverse for each block 
        model.SigWi(jj).mat=inv(cg);
    end
  end

  if (do.theta || do.betaU || do.lamUz)
    SigUW(pu).mat=[];
    for jj=1:pu
        SigUW(jj).mat=covarioGram(model.xzDist,betaU(:,jj),lamUz(jj));
    end
    model.SigUW=SigUW;
  else
    SigUW=model.SigUW;
  end
  
  % The computation is decomposed into the likelihood of W, 
  %  and the likelihood of VU|W. 

  % Compute the likelihood of the W part (already done the blocks)
    LogLikW=sum(model.SigWl);

  % Compute the likelihood of the VU|W
  % This requires using the gaussian model estimation stuff.
    % It is complicated because of shortcuts allowed by lack of correlation
    % between W and V

    % do these ops, on the block diagonal blocks: 
    %    W=SigUW*model.SigWi;
    %    SigUgW=SigU-W*SigUW';
      W(pu).mat=[];
      SigUgW(pu).mat=[];
      for ii=1:pu
        W(ii).mat=SigUW(ii).mat*model.SigWi(ii).mat;
        SigUgW(ii).mat=SigU(ii).mat-W(ii).mat*SigUW(ii).mat';
      end

    %for scalar output: SigVUgW=[SigV+SigUgW] ...
    %                          + model.SigObs/lamOs;
    %otherwise:         SigVUgW=[SigV             zeros(n*pv,n*pu); ...
    %                            zeros(n*pu,n*pv) SigUgW         ] ...
    %                          + model.SigObs/lamOs;
      SigVUgW=model.SigObs/lamOs;
      for ii=1:pv
        blkRange=(ii-1)*n+1:ii*n;
        SigVUgW(blkRange,blkRange)=SigVUgW(blkRange,blkRange)+ ...
                                                SigV(lamVzGroup(ii)).mat;
      end
      if model.scOut
        for ii=1:pu
          blkRange=(ii-1)*n+1:ii*n;
          SigVUgW(blkRange,blkRange)=SigVUgW(blkRange,blkRange)+SigUgW(ii).mat;
        end
      else
        for ii=1:pu
          blkRange=n*pv+[(ii-1)*n+1:ii*n];
          SigVUgW(blkRange,blkRange)=SigVUgW(blkRange,blkRange)+SigUgW(ii).mat;
        end
      end

    % do this op: MuVUgW =W*model.w;
      MuVUgW=zeros(n*pu,1);
      for ii=1:pu
        blkRange1=(ii-1)*n+1:ii*n;
        blkRange2=(ii-1)*m+1:ii*m;
        MuVUgW(blkRange1)=W(ii).mat*model.w(blkRange2);
      end
    
    % for scalar output:  MuDiff=   [u] - [MuVUgW]
    % otherwise:          MuDiff= [v;u] - [0;MuVUgW] 
      if model.scOut
         MuDiff=model.u; MuDiff=MuDiff-MuVUgW;
      else
         MuDiff=model.vu; MuDiff(pv*n+1:end)=MuDiff(pv*n+1:end)-MuVUgW;
      end

    % Now we can get the LL of VU|W
      LogLikVUgW=doLogLik(SigVUgW,MuDiff);
  
  % Final Answer, LL(VU) = LL(VU|W) + LL(W)
    model.logLik=LogLikVUgW+LogLikW;
  

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function L=doLogLik(Sigma,data)
    chCov=chol(Sigma);
    logDet=sum(log(diag(chCov)));
    p1=(chCov')\data;
    L=-logDet-0.5*(p1'*p1);
  end

end
