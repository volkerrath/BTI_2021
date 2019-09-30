function [ierr]=SITE_InvPar(name)

% Prepare Data and Model for inversion

yeartosec=31557600;sectoyear=1/yeartosec;

ierr = 0;

load('common.mat')


% SET PATHS

addpath([srcpath,'/src']);
addpath([srcpath,'/tools']);
addpath([datpath]);
%
% dfmt=1;ffmt='.zip';
% archive(mfilename,strcat([mfilename '_' datestr(now,dfmt)]),ffmt);

% GENERAL SETTINGS

%
site       = 'SYNB';
props       = 'syn';

% DEFINE LOGARITHMIC GSTH INVERSION GRID
nsteps          =   21;                     % number of steps
base            =   0.;                     % base
tstart          =  110000*yeartosec;        % from
tend            =  30*yeartosec;            % to
m_apr_set       = 1;                        % prior value
m_ini_set       = 1;                        % initial value

% SENSITVITY CALCULATION PARAMETERS
diffmeth        =   'FD';                   % how to calculate Sensitivity
dp              =   0.001;                  % FD perturbation
% diffmeth        =   'AD';                 % AutoDiff (disabled)

% ITERATION PARAMETER
tol_solve       =   0.00001;                % tolerance for solver
maxiter_solve   =   32;                     % number of solver iterations
tol_inv         =   [0.0001,0.00001];       % tolerance for inversion rms
maxiter_inv     =   150;                    % maximal number of iterations

% REGULARIZATION PARAMETER FOR TIKHONOV
% reg_opt           =  'UPRE';          % Unbiased Predictive Risk Estimator
% reg_opt           =  'GCV';           % Generalized Cross-Validation
% reg_opt           =  'MLE';           % Maximum Likelihood  (disabled)
% reg_opt           =  'LC';            % L-Curve (disabled)
% reg_opt         =   'FIX';            % fixed regularization parameter
reg_opt         =   'GCV';              % fixed regularization parameter
start_regpar    =   25;                  % when to start search
modul_regpar    =    1;                 % how often
mregpar_adaptint =   1;
%
msteps_regpar    =   48;                % number of test values
regpar0=[0.001 1 0];
reg0par=[0.001];                         % logspace(-1,1,10);31
reg1par=logspace(-3.,3,msteps_regpar);
reg2par=[0];

% RELAXATION
relax           =   1;                  % 0.98;
start_relax     =   2;
modul_relax     =   1;
min_relax       =   0.1;
reg_shift=1;

outsteps=0;

%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% VARIABLES OUTSIDE TMP_PREP OVERWRITE DEFAULTS ABOVE!
F=strcat([name,'_InvPar_in.mat']);
if exist(F)
    disp([mfilename ' defaults overwritten!'])
    load(F); mstruct(invpar_in);
end
%!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



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


