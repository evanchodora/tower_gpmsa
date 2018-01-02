% function d = gendist(data,dataDesc);
%   generates the nxnxp distance array values and supporting
%   information, given the nxp location matrix x
%   or if a d is passed in, just update the distances

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

function d = genDist(data,dataDesc)

d.type=1;

[n p] = size(data);

%generate the list of (n-1)*(n-1) distance indices
  inds=n*(n-1)/2;
  indi=zeros(inds,1);indj=zeros(inds,1);
  ind=1;for ii=1:n-1; indi(ind:ind+n-ii-1)=ii; indj(ind:ind+n-ii-1)=ii+1:n; ind=ind+n-ii; end;
 
  d.n=n; d.p=p;
  d.indi=indi; d.indj=indj;
  d.indm=indi + n*(indj-1);

if p==0; d.d=[]; return; end

if ~exist('dataDesc','var')    
   cat=zeros(1,p);
else    
   cat=[dataDesc.typeCategorical(1:p)];
end
cat0=find(~cat);
cat1=find(cat);


  if isempty(cat1)
    d.d=(data(indj,:)-data(indi,:)).^2;
  else
    d.d=zeros(inds,p);
    d.d(:,cat0)=(data(indj,cat0)-data(indi,cat0)).^2;
    d.d(:,cat1)=(data(indj,cat1)~=data(indi,cat1))*0.5;
  end
