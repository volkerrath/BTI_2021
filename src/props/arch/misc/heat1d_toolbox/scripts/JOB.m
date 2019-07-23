% SCRIPT FOR 
clear;close all;clc


borehole={'Kola-2400','Kola-3200','Kola-3209','Kola-3356'};
% borehole={'SG3_97','2400','2915','3200','3209','3356','3359'};



% GENERATE SPATIAL MESH 
zstart=10;zend=6000;nz=301;type='logarithmic';dir=1;
disp([' ']);disp([' ...set up ' type ' spatial mesh ']); 
[z,dz]=set_mesh(zstart,zend,nz,type,dir,1);

% GENERATE TEMPORAL MESH
year2sec=31557600;
tstart=150000*year2sec;tend=10*year2sec;nt=257;type='logarithmic';dir=-1;
disp([' ']);disp([' ...set up ' type ' temporal mesh ']); 
[t,dt]=set_mesh(tstart,tend,nt,type,dir,1);

% PARAMETER FOR TRANSIENT FWD CALCULATIONS
timestep.theta          = 0.5*ones(1,nt);    % time stepping control parameter,
timestep.dt             = dt;    % time stepping control parameter,
timestep.theta          = 0.5*ones(1,nt);    % time stepping control parameter,
timestep.theta(1:10)    = 1.;                % (.5 = Crank-Nicholson, 1=Backward Euler)
timestep.T0(1:nz)       = 0.                 % initial temperatures


% PARAMETER FOR NONLINEAR ITERATION 
nonlinear.mean      = 'g';              % n/a/g/s calculation of effective properties
                                        % (fix, arithmetic, geometric, square-root)
nonlinear.maxiter   = 2;                % maximal number of nonlinear Picard iterations
nonlinear.tol       = 0.001;            % stopping tolerance for Picard iterations
                                        % (based on max deviation from last iteration)
nonlinear.relax     = 1.                % relaxation parameter for Picard iterations
                                        % (usually 0 < relax <1)  

% DEFINE LOGARITHMIC INVERSION GRID

% PARAMETERS FOR INVERSION AND JACOBIAN CALCULATIONS
nsteps=32; base=0.;
[Ts,it]=set_paleo_grid(t,base,tstart,tend,nsteps);
np=length(pt);
model.Ts          = Ts;             % values of GSTH steps
model.it          = it;             % pointer associating pt values to temporal grid points   	
jacpal.dp          = 1;             % perturbation for FD jacobian calculations
jacpal.maxitdp     = 1;             % max.number of FD refinements for jacobian calculations
jacpal.toldp       = 0.01;          % stoppping tolerance for jacobian calculations
inverse.solver      = 'cgls';        % solver
inverse.tolrms      = 0.5;           % tolerance for rms
inverse.maxiter     =16;             % maximal number of iterations
inverse.tau0 = 0.001;                % intial/center regularization parameters
inverse.tau1 = 5;                    % Gradient
inverse.tau2 = .0;                   % Laplacian
inverse.cgls_iter = 12;              % number of CGLS iterations
inverse.cgls_reorth = 1;             % number of CGLS iterations
inverse.cgls_tol    = 0.0001;        % stoppping tolerance for CGLS iterations
inverse.cgls_stoprule = 1;           % 1/std  2/ACB  3/MULTI
inverse.m_apr       = zeros(1:np);   % priors for GSTH  
inverse.m           = zeros(1:np);   % actual model   
% PARAMETERS FOR PERMAFROST MODELING
freeze.flag         = 'yes';        % yes/no
freeze.Tf           = 0.;           % freezing temperature (default 0C)
freeze.w            = 1.;           % width of freezing interval (in K, ususlly 1)
freeze.Lh           = 333600;       % latent heat of freezing
% END OF PARAMETER DEFINITIONS




% SECIFIC SETTINGS FOR THIS RUN
inverse.tau0 = 0.001;                % intial/center regularization parameters
inverse.tau1 = 5;                    % Gradient
inverse.tau2 = .0;                   % Laplacian
JOBNAME=strcat(['Kola-1']);
NTIPAL_Kola;



