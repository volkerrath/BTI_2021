
clear all
close all
clc

% SET PATHS
pltpath='./';
datpath='./';
srcpath='../';
utlpath='../';

addpath([srcpath,filesep,'src']);
addpath([srcpath,filesep,strcat(['tools'])]);

% ONLY FOR PARRALLEL  EXECUTION
run_parallel    =  1;
parcors         =  6;


save('common','srcpath','utlpath','datpath','pltpath','run_parallel','parcors'),



dfmt=1;ffmt='.zip';
archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);
yeartosec=31557600;sectoyear=1/yeartosec;

%GRAPHICS

set_graphpars
%plotfmt='epsc2';
plotfmt='png';

% ONLY FOR PARRALLEL  EXECUTION
run_parallel=1;
parcors=   4;

% SET RANDOM GENERATOR
rng('shuffle');
%randn('state',sum(100*clock));

%
site       = 'OKU';
props       = 'oku';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER FOR FORWARD MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
theta           =   1.;                   % time steping weight 1/FI .5/CN
maxitnl         =   4;                    % number of nl iterations
tolnl           =   0.00001;
relaxnl         =  1.;
freeze          =  1;                     % include freezing/thawing

num_par=mstruct(theta,maxitnl,tolnl,relaxnl,freeze);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% INVERSION PARAMETER
% DEFINE LOGARITHMIC GSTH INVERSION GRID
nsteps          =   24;                     % number of steps
base            =   0.;                     % base
tstart          =  105000*yeartosec;        % from
tend            =  10*yeartosec;            % to
m_apr_set       = 5;                        % prior value
m_ini_set       = 5;                        % initial value

% SENSITVITY CALCULATION PARAMETERS
diffmeth        =   'FD';                   % how to calculate Sensitivity
dp              =   0.001;                  % FD perturbation
% diffmeth        =   'AD';                 % AutoDiff (disabled)

% ITERATION PARAMETER
tol_solve       =   0.00001;                % tolerance for solver
maxiter_solve   =   32;                     % number of solver iterations
tol_inv         =   [0.0001,0.00001];       % tolerance for inversion rms
maxiter_inv     =   100;                    % maximal number of iterations

% REGULARIZATION PARAMETER FOR TIKHONOV
% reg_opt           =  'UPRE';          % Unbiased Predictive Risk Estimator
% reg_opt           =  'GCV';           % Generalized Cross-Validation
% reg_opt           =  'MLE';           % Maximum Likelihood  (disabled)
% reg_opt           =  'LC';            % L-Curve (disabled)
% reg_opt         =   'FIX';            % fixed regularization parameter
reg_opt         =   'GCV';              % fixed regularization parameter
start_regpar    =    5;                  % when to start search
modul_regpar    =    1;                 % how often
mregpar_adaptint =   1;
%
msteps_regpar    =   48;                % number of test values
regpar0=[1 1 1];
reg0par=[0.1];                         % logspace(-1,1,10);31
reg1par=logspace(-3.,3,msteps_regpar);
reg2par=[0];

% RELAXATION
relax           =   1;                  % 0.98;
start_relax     =   2;
modul_relax     =   1;
min_relax       =   0.1;
reg_shift=1;

outsteps=0;

for Qb=-[40 41 42 43 44]*1e-3
    prepstr       = ['_GCV_Q' num2str(abs(Qb*1000))];
    
    name=[site prepstr]
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERATE MESHES
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % VARIABLES SET HERE OUTSIDE OKU_MESH OVERWRITE DEFAULTS INSIDE!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % F=strcat([name,'_mesh_in.mat']);
    % mesh_in=mstruct();
    % save(F,'mesh_in');
    disp(strcat([' generate meshes for ' name]));
    C=strcat([site,'_Mesh(name);']);
    eval(C);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % GENERATE PHYSICAL MODEL
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % VARIABLES SET HERE OUTSIDE OKU_PREP OVERWRITE DEFAULTS INSIDE!
    %!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    plotit=0;
    prep_in=mstruct(plotit,Qb);
    F=[name,'_prep_in'];
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
    plotit=0;
    init_in=mstruct(plotit,Qb);
    F=[name,'_init_in'];
    save(F,'init_in');
    disp(strcat([' generate initial values for ' name]));
    C=strcat([site,'_Init(name);']);eval(C);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % START inversion
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    fwdpar=mstruct(theta,maxitnl,tolnl,relaxnl,freeze);
    F=strcat([name,'_FwdPar.mat']);
    save (F,'fwdpar')
    
    
    F=strcat([name,'_TimeGrid.mat']);load(F);
    [gsth,pt]=set_mgsth(t,base,tstart,tend,nsteps);
    
    invpar=mstruct(...
        gsth,pt,nsteps,m_apr_set,m_ini_set,...
        diffmeth,dp,...
        tol_solve,maxiter_solve,...
        tol_inv,maxiter_inv,...
        reg_opt,start_regpar,modul_regpar,mregpar_adaptint,msteps_regpar,...
        regpar0,reg0par,reg1par,reg2par,reg_shift,...
        relax,start_relax,modul_relax,min_relax,outsteps);
    
    F=strcat([name,'_InvPar.mat']);
    save (F,'invpar');
    
    
    Tikh_gsth(name);
    
end