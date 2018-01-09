% commands to run the tower example...
addpath(fullfile('code'))

% read data
dat=sc(1);

% specify scalar output
optParms.scalarOutput=1;

% initial set-up
optParms.priors.lamOs.a=1000; 
optParms.priors.lamOs.b=1000;
params=setupModel(dat.obsData,dat.simData,optParms);

% modifications to defaults
params.priors.lamVz.a=1;
params.priors.lamVz.b=0.0001;
params.priors.lamWs.a=1;
params.priors.lamWs.b=0.0001;
params.priors.lamWOs.a=2;
params.priors.lamWOs.b=1e-6;

% modifications to initial values
params.model.lamWOs=1e6;
params.model.lamOs=1;

% modifications to initial step widths
params.mcmc.lamWOswidth=0;

% step size
nburn=500; nlev=21;
params=gpmmcmc(params,0,'stepInit',1,'nBurn',nburn,'nLev',nlev);
params=stepsize(params,nburn,nlev);

% mcmc
nmcmc=10000;
pout=gpmmcmc(params,nmcmc,'step',1);
save pout pout;

nmcmc=nmcmc+nburn*nlev;
pvec=floor(linspace(nburn*nlev+1,nmcmc,500));
pout.pvec=pvec;

% plots
% load pout;
scPlots(pout,pvec,1:4);

% calibration parameters
tmp = [pout.pvals.theta]';
save 'theta' tmp '-ascii';

% % sensitivity analysis
% sens=gSens(pout,'pvec',pvec);
% 
% pout.sens=sens;
% save pout pout;
% 
% tmp=[sens.smePm;sens.stePm]';
% save 'sa_pm' tmp '-ascii';
% 
% scPlots(pout,pvec,5);
