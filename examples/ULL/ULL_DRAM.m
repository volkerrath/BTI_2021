
clear all
close all
clc


% SET PATHS
pltpath='./';
datpath='./';
srcpath='../';
utlpath='../';

addpath([srcpath,'src']);
addpath([srcpath,'src/mcmc/']);
addpath([srcpath,'tools'])
addpath([datpath])

yeartosec=31557600;sectoyear=1/yeartosec;


parjobs=4;

mypool=gcp;
if isempty(mypool)
    mypool=parpool(parjobs);
else
    delete(mypool);    
    mypool=parpool(parjobs);
end 


dfmt=1;ffmt='.zip';
archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);
yeartosec=31557600;sectoyear=1/yeartosec;


% BOREHOLES8
site        = 'ULL';
props       = lower(site);    %

addpath([srcpath,strcat(['src/props/',lower(props)])])


meth        = 'DRAM';
prepstr       = '';
CovFun='G'; L=3;


%GRAPHICS
set_graphpars
%plotfmt='epsc2';
plotfmt='png';

name=strcat([site prepstr,'_GSTQ_',meth,'_',CovFun,num2str(L),'_',datestr(now,1)]);

save('common',...
    'srcpath','utlpath','datpath','parjobs','name');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  NUMERICAL PARAMETERS FOR FWD SOLUTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta           =   1.;         % time steping weight 1/FI .5/CN
% NONLINEAR ITERATION
maxitnl         =   4;                          % number of nl iterations
tolnl           =   0.00001;
relaxnl           =  1.;
freeze          =  1;                           % include freezin/thawing


numpar=mstruct(theta,maxitnl,tolnl,relaxnl,freeze);
F=strcat([name,'_NumPar.mat']);
save (F,'numpar')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER FOR INVERSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% DEFINE LOGARITHMIC INVERSION GRID

disp(['   ']);
disp([' ...set up parametrization for paleoclimate inversion  ']);

nsteps          =  21;                 % number of steps
base            =   0.;                 % base
tstart          =  105000*yeartosec;         % from
tend            =  30*yeartosec;             % to



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE MESHES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES SET HERE OUTSIDE OKU_MESH OVERWRITE DEFAULTS INSIDE!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% F=strcat([name,'_Mesh_in.mat']);
% nt = 401; nz = 301;% zend = 5000;
% mesh_in=mstruct(nz,zend,nt,tstart,tend);
% save(F,'Mesh_in');
disp(strcat([' generate meshes for ' name]));
C=strcat([site,'_Mesh(name);']);
eval(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE PHYSICAL MODEL AND DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES SET HERE OUTSIDE OKU_PREP OVERWRITE DEFAULTS INSIDE!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plotit=0;
prep_in=mstruct(plotit,props);
F=[name,'_Prep_in'];
save(F,'prep_in');
disp(strcat([' generate model for ' name]));
C=strcat([site,'_Prep(name);']);
eval(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE INITIAL VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES SET HERE OUTSIDE OKU_INIT OVERWRITE DEFAULTS INSIDE!
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

disp(strcat([' generate initial values for ' name]));
% plotit=0;
% init_type='e';
% GSTH0=7.;
% POM=-4.;
% init_in=mstruct(plotit,init_type,pom,sval,props);
% F=[name,'_Init_in'];
% save(F,'init_in');
C=strcat([site,'_Init(name);']);eval(C);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MCMC CTRL PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F=strcat([name '_SitePar.mat']);load(F);
DATA.sitepar    = sitepar;mstruct(sitepar)
F=strcat([name,'_NumPar.mat']);load(F);
DATA.numpar     = numpar;
F=strcat([name,'_Init.mat']);load(F)
DATA.initial    = Tinit;

load([name,'_TimeGrid.mat']);
[gsth,pt]=set_mgsth(t,base,tstart,tend,nsteps);
gsthpar=mstruct(gsth,pt,nsteps);
F=strcat([name,'_GSTHPar.mat']);save (F,'gsthpar');
save (F,'gsthpar');
DATA.gsthpar = gsthpar;


% SELECT THE METHOD USED
method = lower(meth); % 'dram' 'mh','am','dr', or 'dram', see below

% PARAMETERS FOR MCMC
switch method
    case 'mh'
        nsimu    = 100000; drscale  = 0; adaptint = 0; 
    case 'dr'
        nsimu    = 200000; drscale  = 2; adaptint = 0; 
    case 'am'
        nsimu    = 500000; drscale  = 0; adaptint = 25000; 
    case 'dram'
        nsimu    = 200000; drscale  = 2; adaptint = 10000; 
end


% [RESULTS,CHAIN,S2CHAIN,SSCHAIN] = MCMCRUN(MODEL,DATA,PARAMS,OPTIONS)
MODEL.ssfun         = @GSTH_MCObjFunGauss;
MODEL.N             = length(id);
MODEL.sigma2        = 0.1^2;

OPTIONS.updatesigma = 1;        % update error variance. s2 sampled with updatesigma=1
OPTIONS.nsimu       = nsimu;                       % size of the chain


OPTIONS.burnintime  = 1;        % burnin is dealt with from full chain
OPTIONS.verbosity   = 1;

OPTIONS.printint    =  100;                        % output interval
OPTIONS.saveint     =  300;                        % save interval
OPTIONS.residout    =  1;




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE MCMC GSTH / PARS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% SET RANDOM GENERATOR
%rng('simdTwister')
rng('shuffle');
%randn('state',sum(100*clock));
%
GSTprior = 4.; GSTerr=2.5;
Gprior=[GSTprior*ones(nsteps,1)]; ns=length(Gprior);
GErr=GSTerr*ones(size(Gprior(:)'));
Qprior = 42.;QErr = 3;
Hprior  = 1.5;HErr = 0.3;
Pprior   = [Gprior'  Qprior  Hprior ];
Perr=[ GErr QErr HErr];

Pcutoff = 3.;
Pmin    = Pprior-Pcutoff*Perr;Pmax=Pprior+Pcutoff*Perr;
Pavg    = Pprior+0.5*Perr.*randn(size(Pprior));

% active parameters
Pact    = 1*ones(size(Pavg)); Pact(ns+1)=1; Pact(ns+2)=1;
DATA.mactive = Pact;
np = length(Pact);

switch lower(CovFun)
    case{'e' 'markov' 'exponential'}
        C0G=CovarExpnl(GErr,L);
    case{'g' 'gauss'}
        C0G=CovarGauss(GErr,L);
    case{'m' 'matern'}
        C0G=CovarMatern(GErr,L);
    otherwise
        error([lower(CovFun), ' not implemented. STOP']);
end
COQ=[QErr^2];
C0H=[HErr^2];
C0=blkdiag(C0G,COQ,C0H);
OPTIONS.qcov     = C0;     % initial proposal covariance


for iparam =1:np
    PARAMS{iparam} =  ...
        {['par' num2str(iparam)],...
        Pavg(iparam),Pmin(iparam),Pmax(iparam),Pprior(iparam),Perr(iparam)};
end


mcmcpar=mstruct(MODEL,DATA,PARAMS,OPTIONS);
F=strcat([name,'_MCMCPar.mat']);
disp([' ']);disp([' >>>>> mcmc configuration written to: ' F]);
save (F,'mcmcpar');


parfor job=1:parjobs
        RunMCMC(job,name,DATA,OPTIONS,PARAMS,MODEL);
end

delete(mypool)
