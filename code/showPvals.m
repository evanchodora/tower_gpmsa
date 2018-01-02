% function showPvals(pvals, skip)
%   skip = the beginning index to display; optional

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

function showPvals(pvals, skip)

if ~exist('skip'); skip = 1; end;

fprintf('Plotting pval struct from index %d to %d\n',skip,length(pvals));

f=fieldnames(pvals);
flen=length(f);

x=skip:length(pvals);
pvals=pvals(skip:end);

clf
for ii=1:flen
     y=[pvals.(f{ii})];
     h(ii)=pvalSubplot(flen,ii);
     if length(x)==length(y)
       plot(x,y);
     end
     ylabel(f{ii});
end
set(h(1:end-1),'XTick',[]);

end


% An internal function is needed to get the subplots to use more of the
% available figure space

function h=pvalSubplot(n,i) 

sep=0.25;
left  =0.1;
sizelr=0.8;
bottom=(1-(i/n))*0.8+0.1;
sizetb=(1/n)*0.8*(1-sep/2);

h=subplot('position',[left bottom sizelr sizetb]);
end
