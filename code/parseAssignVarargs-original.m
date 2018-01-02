% function parseAssignVarargs(validVars)
% assigns specified caller varargs to the corresponding variable name
%   in the calling workspace. vars not specified are not assigned.
% validVars is a cell array of strings that represents possible
%    arg names, and the variable name in the workspace (identical)
% varargs is the varargin cell array to parse

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

function parseAssignVarargs(validVars)

pr=0;

bargs=evalin('caller','varargin');

for ii=1:2:length(bargs)
  varID=find(strcmp(validVars,bargs{ii}));
  if length(varID)>1; error('bad parse list to parseAssignVarargs'); end
  if pr; fprintf('Assigning: %s\n', validVars{varID}); bargs{ii+1}, end 
  assignin('caller',validVars{varID},bargs{ii+1});
  
end
 