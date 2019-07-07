% PARAMETERS FOR INVERSION 
clear all;close all;clc
debug=1;
iterout = 0;

disp([' ELK1 - Permafrost '])
disp([' run at: ',date])
disp([' ']);

%==========================================================================
% set general and default parameter 
%==========================================================================

% setup model
zstart=10;zend=5000;nz=251;type='logarithmic';dir=1;
disp([' ']);disp([' ...set up ' type ' spatial mesh ']);
[z,dz]=set_mesh(zstart,zend,nz,type,dir,0);nc=nz-1;nu=nc;
% setup timesteps
year2sec=31557600;
tstart=100000*year2sec;tend=10*year2sec;nt=257;type='logarithmic';dir=-1;
disp([' ']);disp([' ...set up ' type ' temporal mesh ']);
[t,dt]=set_mesh(tstart,tend,nt,type,dir,0);ns= nt- 1;
% set inversion grid
disp(['   ']);disp([' ...set up parametrization for paleoclimate inversion  ']);
nmod=24; base=0.;
[mod,it]=set_paleo_grid(t,base,tstart,tend,nmod);np= length(mod);
% set prior 
step_t=[100000 15000 200];step_T=[base base base];L=5;
[prior]=set_prior(step_T,step_t,L,tstart,tend,nmod,0);


% store into structures & set default values
% model
model.k(1: nu)          = 3.;
model.kA(1: nu)         = 0.;
model.kB(1: nu)         = 0.;
model.p(1: nu)          = 0.1;
model.h(1: nu)          = 1.e-6;
model.r(1: nu)          = 2300;
model.c(1: nu)          = 1000;
model.z                 = z;
model.dz                = dz;
model.nz                = length(z);
model.ip(1: nc)         = [1:nc];
model.qb                = .060;
model.Ts                = mod(it);
model.it                = it;

% time step
timestep.dt             =dt;
timestep.nt             =nt;  
timestep.theta(1:ns)    = 1.;
timestep.out        	= 0;

% linear solver
linsolve.solver         = 'direkt';

% fixed point iteration
nonlinear.mean          = 'g';
nonlinear.maxiter       = 2;
nonlinear.tol           = 1.e-4;
nonlinear.relax         = 1.;

%freezing
freeze.flag             = 'y';
freeze.Lh               = 333600;
freeze.Tf               = 0.;
freeze.w                = 1.;

% inversion 
inverse.mod             = mod; 
inverse.np              = length(mod); 
inverse.prior           = prior;
inverse.tolrms          = 0.01;      % tolerance for rms
inverse.maxiter         = 16;        % maximal number of iterations
% cgls control
inverse.maxcgls         = 12;        % maximal number of CGLS iterations
inverse.stoprule        = 2;         % 1= use tol 2=use ACB-criterium
inverse.tolcgls         = .00001;    % CGLS tolerance  
inverse.reorth          = 1;         % reorthogonalization
% regularization parameter
inverse.tau0            = 0.5;       % TiKhonov L_0
inverse.tau1            = 3;         % TiKhonov L_1 (flatness penalty)
inverse.tau2            = 0;         % TiKhonov L_2 (smoothness penalty)
% FD jacobian control
inverse.dp               = 0.01;
inverse.maxitdp          = 1;
inverse.toldp            = 1e-5;
%


%==========================================================================
%  JOB
%==========================================================================


tau0=0.5;
tau1=3;
tau2=0;

inverse.tau0=tau0;inverse.tau1=tau1;inverse.tau2=tau2;
freeze.flg='yes';
borehole={'Elk1_prep'}; bore=borehole{1};  NAME=strcat(['ELK1-y-' num2str(tau0) '-'  num2str(tau1) '-'  num2str(tau2)]); NEW_TIPAL;
freeze.flg='no';
borehole={'Elk1_prep'}; bore=borehole{1};  NAME=strcat(['ELK1-n-' num2str(tau0) '-'  num2str(tau1) '-'  num2str(tau2)]); NEW_TIPAL;



