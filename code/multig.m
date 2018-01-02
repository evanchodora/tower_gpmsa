%function params=multig(params,clist,nmcmc,varargin)
% Run simultaneous GPM models with simultaneous calibration of shared
% variables

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

function params=multig(params,clist,nmcmc,varargin)

% Process input arguments
stepInit=0; step=0; nBurn=0; nLev=0;
parseAssignVarargs({'step','stepInit','nBurn','nLev'});
if stepInit; nmcmc=nBurn*nLev; end

pr=0; % print diagnostics

% get count of params structures
  models=length(params);

% The clist structure contains lists of thetas that are in common.
% It is a matrix, each row is for one shared variable (theta). A row is a
% list of indicators the same length as the number of models. A zero in a
% corresponding position indicates it is not in the corresponding model, a
% nonzero entry indicates the index of the variable. 

%initialize new structure as necessary
  for ii=1:models;
    params(ii).vLinks.sharedIndices=[];
    params(ii).vLinks.skipList=[];
    params(ii).vLinks.vars=[];
  end

% Set up params' reference structures
for ii=1:size(clist,1); %number of cross-model variables

  mLinks=find(clist(ii,:)); % these are the models linked by that var
  
  %first, push out common value for this shared variable to all models
    for jj=2:length(mLinks)
      initialVal=params(mLinks(1)).model.theta(clist(ii,mLinks(1)));
      params(mLinks(jj)).model.theta(clist(ii,mLinks(jj)))=initialVal;
      params(mLinks(jj)).vLinks.skipList=...
          [params(mLinks(jj)).vLinks.skipList clist(ii,mLinks(jj))];
    end

  % Now, build the link structure for the model managing the variable
    jj=mLinks(1); %the manager is the first one
    otherMLinks=mLinks(2:end); % the list of other models

    % add a struct in each models params structure detailing the index of
    % the linked models and the variable ID within that model
    % vIndices are the linked models' indices
    % vars.links are, for each model link, the theta indices
    params(jj).vLinks.sharedIndices= ...
         [params(jj).vLinks.sharedIndices clist(ii,jj)];
    vCount=length(params(jj).vLinks.sharedIndices);
    params(jj).vLinks.vars(vCount).links=...
         [otherMLinks' clist(ii,otherMLinks)'];
    params(jj).vLinks.logLikCorrHandle=@getExtLikVal;
    params(jj).vLinks.pushValHandle=@updateModels;
  
end

if pr; interpretVLinks(params); end

% Initialize all model partial computations so that linking is valid
for ii=1:models
  params(ii) = gpmmcmc(params(ii),0,'initOnly',1);
  % modelRef is a temporary place to hold a computed linked model,
  %  until the draw is accepted by the master model (or not)
  modelRef(ii)=params(ii).model;
  if stepInit, params(ii).mcmc.base=[]; end
end

% Run the models. They will callback the getExtLikVal nested function when
% they need to compute the external links
 counter('stime',1,nmcmc,5,10);
 for ii=1:nmcmc
  counter(ii);
  for jj=1:models
    if pr; fprintf('Iteration %d, model %d\n',ii,jj); end;
    params(jj)=gpmmcmc(params(jj),1,'linkedMod',1,'noInit',1,...
      'noCounter',1,'stepInit',stepInit,'step',step,'nLev',nLev);
  end
 end
 counter('end');

% Push out mcmc acceptance tracking vectors for shared variables
 if stepInit
   for ii=1:models
     sInd=params(ii).vLinks.sharedIndices;
     for jj=1:length(sInd)
       vLinks=params(ii).vLinks.vars(jj).links;
       for kk=1:size(vLinks,1)
         modInd=vLinks(kk,1);
         modVar=vLinks(kk,2);
         params(modInd).mcmc.acc.theta(:,modVar)= ...
         params(ii).mcmc.acc.theta(:,sInd(jj));
       end
     end
   end
 end

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function LLval=getExtLikVal(vLinks,newVal);
  % Nested function callback from function handle
  % Compute a loglikelihood correction for linked variables

    LLval=0; LL1=0;
    for kk=1:size(vLinks,1)
      modInd=vLinks(kk,1); 
      modVar=vLinks(kk,2);

      if pr; fprintf('\nComputing link to model %d var %d\n',modInd,modVar); end

      model=params(modInd).model;
      LL1=model.logLik;

      if pr; fprintf('theta from %f to %f\n',model.theta(modVar),newVal); end;

      C.var='all';
      model.theta(modVar)=newVal;
      modelRef(modInd)=computeLogLik(model,params(modInd).data,C);
      LL2=modelRef(modInd).logLik;

      LLval=LLval+(LL2-LL1);

      if pr; fprintf('Result is %f-%f=%f\n',LL2,LL1,LLval); end;
    end
  end %nested function

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function LLval=updateModels(vLinks,newVal);
  % Nested function callback from function handle
  % Push an accepted shared variable update to the linked models

    for kk=1:size(vLinks,1)
      modInd=vLinks(kk,1);
      modVar=vLinks(kk,2);
      if pr
        fprintf('Updating from model %d to model %d, var %d, link index %d\n', ...
               jj,modInd,modVar,kk);
      end
      params(modInd).model=modelRef(modInd);
    end
  end %nested function

end %m-file function


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function interpretVLinks(params)
% Interprets and displays the variable link structure in a multi-model
% gasp setup, for diagnostic purposes. 

  fprintf('\n')

  for ii=1:length(params)
    fprintf('Model %d linked variables:', ii);
    fprintf(' %d',params(ii).vLinks.sharedIndices);

    fprintf('\n');

    for jj=1:length(params(ii).vLinks.vars)
      fprintf(' Index %d \n',jj);

      for kk=1:size(params(ii).vLinks.vars(jj).links,1);
        fprintf('  Model %d variable %d\n',...
            params(ii).vLinks.vars(jj).links(kk,1), ...
            params(ii).vLinks.vars(jj).links(kk,2) );
      end
    end
  end
  fprintf('\n')
end
