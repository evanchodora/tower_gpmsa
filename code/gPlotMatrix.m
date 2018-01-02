function [h bigAx]=gPlotMatrix(data,varargin)
% function [h bigAx]=gPlotMatrix(data,varargin)
%  data - contains vectors for scatterplots
%    each row is an vector, as expected for plotmatrix
% varargs include
%  'Pcontours' are the percentile levels for the contour plot
%  'ngrid' is axis grid size (symmetric) (a good guess is 25)
%  'labels', a cell array of variable names [optional]
%  'ttl', an overall plot title [optional]
%  'axRange', a 2-vector of axis scalings, default [0,1] or data range
%  'ksd', the sd of the contour smoothing kernel
%  'Pcontours', probability contours, default [0.5 0.9]
%  'ustyle', 'lstyle' is the type of the off-diagonal plots
%     'scatter' is xy scatterplots [default]
%     'imcont' is a smoothed image (2d est. pdf) with contours
%  'shade' causes the scatterplots to to shade from blue to red over
%     the input sequence of points

  [n p]=size(data);

  % defaults for varargs
    labels=[];
    if any((data(:)<0)|any(data(:)>1));
      axRange=[min(data)' max(data)'];
    else
      axRange=repmat([0 1],p,1);
    end
    ksd=0.05;
    ngrid=10;
    ustyle='scatter'; lstyle='scatter';  
    Pcontours=[0.5;0.9];
    shade=0;
    ttl=[];
   parseAssignVarargs({'labels','axRange','ngrid','ksd', ...
                       'Pcontours','ustyle','lstyle','shade','ttl'});  

  Pcontours=Pcontours(:);                      
  ncont=length(Pcontours);
  
  % if shading is enabled, set up the shade structs. 
   if shade
    shgroups=min(n,100); % need to set up groups if n is large
    sls=linspace(1,n,shgroups+1)'; slc=linspace(0,1,shgroups);
    for shi=1:shgroups; % define a range and a color for each group
      sh(shi).ix=ceil(sls(shi)):floor(sls(shi+1)); 
      sh(shi).color=(1-slc(shi))*[0 0 1] + slc(shi)*[1 0 0];
    end
   else
    shgroups=1; sh.ix=1:n; sh.color=[0 0 1];
   end
  
  % Put the data into a standard [0,1] range
    data=(data-repmat(axRange(:,1)',n,1)) ./ ...
              repmat((axRange(:,2)-axRange(:,1))',n,1);  
   
  % Generate a grid and supporting data structures
    gridvals = linspace(0,1,ngrid);
    [g1 g2] = meshgrid(gridvals,gridvals);
    g1v = g1(:); g2v = g2(:);
    gvlen = length(g1v);
    dens = zeros(gvlen,1);

  % begin
    clf;
    
  % establish the subplots
    for ii=1:p; for jj=1:p; 
      h(ii,jj)=gPackSubplot(p,p,ii,jj);
    end; end
  
  % Put in the histograms
    for ii=1:p
      axes(h(ii,ii));
      hist(data(:,ii)); 
      axisNorm(h(ii,ii),'x',[0 1]);
    end
  
   % Go through the 2D plots
   for ii=1:p-1; for jj=ii+1:p
     % compute the smooth and contours, if it's called for either triangle
     if any(strcmp({ustyle,lstyle},'imcont'))
       % compute the smooth response
         for i=1:gvlen
           f = normpdf(data(:,jj),g1v(i),ksd).*normpdf(data(:,ii),g2v(i),ksd);
           dens(i) = sum(f);
         end
       % normalize dens
         dens = dens/sum(dens);
       % get the contours
         for j=1:ncont
          hlevels(j) = fzero(@(x) getp(x,dens)-Pcontours(j),[0 max(dens)]);
         end
     end

     %% Do the upper triangle plots
       axes(h(ii,jj));
       switch ustyle
       case 'scatter'
         for shi=1:shgroups
           plot(data(sh(shi).ix,jj),data(sh(shi).ix,ii),'.', ...
                'MarkerSize',6,'Color',sh(shi).color);
           hold on;
         end
       case 'imcont'
         imagesc(g1v,g2v,reshape(dens,ngrid,ngrid)); axis xy; hold on;
         colormap(repmat([.9:-.02:.3]',[1 3]));
         contour(g1,g2,reshape(dens,ngrid,ngrid),hlevels,'LineWidth',1.0,'Color','b'); 
       otherwise
         error('bad specification for lstyle');
       end
       axis([0 1 0 1]);

     %% Do the lower triangle plots
       axes(h(jj,ii)); 
       switch lstyle
       case 'scatter' 
         for shi=1:shgroups
           plot(data(sh(shi).ix,ii),data(sh(shi).ix,jj),'.', ...
                'MarkerSize',6,'Color',sh(shi).color);
           hold on;
         end
         hold on;
       case 'imcont'
         imagesc(g1v,g2v,reshape(dens,ngrid,ngrid)'); axis xy; hold on;
         colormap(repmat([.9:-.02:.3]',[1 3]));
         contour(g1,g2,reshape(dens,ngrid,ngrid)',hlevels,'LineWidth',1.0,'Color','b'); 
       otherwise
         error('bad specification for lstyle');
       end
       axis([0 1 0 1]);
     
   end; end

   set(h,'XTick',[],'YTick',[]);
   %for ii=1:p;
   %   set(h(end,ii),'XTick',axRange(ii,:));
   %end

   % labels
   if ~isempty('labels')
     for ii=1:p
       title(h(1,ii),labels{ii});
       ylabel(h(ii,1),labels{ii});
     end
   end
   
   % if a title was supplied, put it up relative to an invisible axes
    if ~isempty(ttl)
      bigAx=axes('position',[0.1 0.1 0.8 0.8],'visible','off'); hold on;
      text(0.5,1.05,ttl,'horizontalalignment','center','fontsize',14);
    end 
   
end


% function to get probability of a given level h
function pout = getp(h,d);
    iabove = (d >= h);
    pout = sum(d(iabove));
end

function parseAssignVarargs(validVars)
% assigns specified caller varargs to the corresponding variable name
%   in the calling workspace. vars not specified are not assigned.
  bargs=evalin('caller','varargin');
  for ii=1:2:length(bargs)
    varID=find(strcmp(validVars,bargs{ii}));
    if length(varID)>1; error('bad parse list to parseAssignVarargs'); end
    assignin('caller',validVars{varID},bargs{ii+1});
  end
end
   
