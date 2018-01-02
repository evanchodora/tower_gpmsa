% function d = gendist2(data1,data2,dataDesc);
%   generates the nxmxp distance array values and supporting
%   information, given the nxp matrix data1 and mxp data2

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

function d = genDist2(data1,data2,dataDesc)

d.type=2;
[n p1] = size(data1);
[m p2] = size(data2);
p=max(p1,p2);

%generate & store the list of n*m distance indices
  inds=n*m;
  indi=repmat(1:n,1,m);
  indj=repmat(1:m,n,1); indj=indj(:)';
 
  d.n=n; d.m=m; d.p=p;
  d.indi=indi; d.indj=indj;
  d.indm=indi + n*(indj-1);

if any([p1 p2]==0); d.d=[]; return; end % if either dataset is empty

if ~exist('dataDesc','var'); cat=zeros(1,p);
                  else cat=[dataDesc.typeCategorical(1:p)];
end
cat0=find(~cat); cat1=find(cat);


if isempty(cat1);
  d.d=(data1(indi,:)-data2(indj,:)).^2;
else
  d.d=zeros(inds,p);
  d.d(:,cat0)=(data1(indi,cat0)-data2(indj,cat0)).^2;
  d.d(:,cat1)=(data1(indi,cat1)~=data2(indj,cat1))*0.5;
end
