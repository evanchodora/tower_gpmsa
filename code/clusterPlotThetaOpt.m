function [clust trec]=clusterThetaOpt(trec)

% if there is no theta record, then create some
if ~exist('trec')
  for ii=1:100; 
    trec(ii)=optimizeTheta(poutH(1), ...
                           'mode','optimize', ...
                           'errType','Likelihood', ...
                           'theta0',rand(1,3)/2+0.25); 
  end
end
  
t=zeros(length(trec),length(trec(1).theta));
t0=zeros(length(trec),length(trec(1).theta0));
err=zeros(length(trec),1);
for ii=1:length(trec);
  t(ii,:)=trec(ii).theta;
  t0(ii,:)=trec(ii).theta0;
  err(ii)=trec(ii).err;
end

D = pdist(t); 
L = linkage(D); 
c = cluster(L,'cutoff',0.1,'criterion','distance');

fprintf('%d clusters\n',max(c));

colors={'b*','g*','r*','k*','c*','m*'};
figure(1); clf;
for ii=1:size(t,1)
   plot3([t(ii,1) t0(ii,1)],[t(ii,2) t0(ii,2)],[t(ii,3) t0(ii,3)]); hold on;
end
for ii=1:max(c)
  tt=t(c==ii,:);
  clust(ii).theta= mean(tt,1);
  clust(ii).err  = mean(err(c==ii));
  plot3(clust(ii).theta(1),clust(ii).theta(2),clust(ii).theta(3), ...
        colors{mod((ii-1),length(colors))+1},'markersize',10); 
end
axis([0 1 0 1 0 1]);
grid on;

figure(2); clf; 
dendrogram(L);
