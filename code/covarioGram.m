% function Scov = covarioGram(dist,beta,lamz,lams)
% given n x p matrix x of spatial coords, and dependence parameters
% beta p x 1, this function returns a matrix that gives the GPM
% correlation
%     Scov_ij = exp{- sum_k=1:p beta(k)*(x(i,k)-x(j,k))^2 } ./lamz

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

function Scov = covarioGram(dist,beta,lamz,lams)

% check for case of a null dataset
%if isempty(dist.d); Scov=[]; return; end
  
switch dist.type
 case 1 
  n=dist.n;
  Scov=zeros(n);
  if n>0  % if it's not a null dataset, do the distances
    % Scov(dist.indm)=exp(-(dist.d*beta))./lamz; %%% this is faster: 
       t=exp(-(dist.d*beta))./lamz;
       Scov(dist.indm)=t;
    Scov=Scov+Scov';
    diagInds = (0:n:n*(n-1)) + (1:n);
    Scov(diagInds)=1/lamz;
    if exist('lams','var')
      Scov(diagInds)=Scov(diagInds) + 1/lams;
    end
  end
 case 2
  n=dist.n; m=dist.m;
  Scov=zeros(n,m);
  if n*m >0 % if it's not a null dataset, do the distances
    %Scov(dist.indm)=exp(-(dist.d*beta))./lamz; %%% this is faster: 
      t=exp(-(dist.d*beta))./lamz;
      Scov(dist.indm)=t;
  end
 otherwise  
  error('invalid distance matrix type in gaspcov');
end
