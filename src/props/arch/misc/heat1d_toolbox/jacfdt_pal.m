function[J,Tc,pointer]=...
    jacfdt_pal(inverse,model,timestep,linsolve,nonlinear,freeze,out)
% JACFDT_PAL calculates transient Jacobians with respect to paleoclimate.
%
% [J]=jacfdt_pal(jacpal,mod,timestep,nonlinear,freeze,out) caclulates
% Jacobians of recent (index nt+1) borehole temperatures with respect to
% paleoclimate values given in pt(it) by forward differencing.
%
%  structure inverse:
%       inverse.mod             = model
%       inverse.dp              = model perturbation
%       inverse.maxitdp         = maximal number of FD refinements for
%                                 jacobiancalculations
%       inverse.toldp           = stopping tolerance for jacobian calculations
%
%  structure mod:
%       model.k (1:nu)          = lambda
%       model.kA, kB(1:nu)      = thermal conductivity coefficient A and B
%       model.h(1:nu)		    = volumetric heat production
%       model.qb			    = basal heat flow
%       model.Ts            	= time-dependent boundary temperatures
%                     		      at the top (e.g., values of GSTH steps
%       model.it                = pointer associating pt values to temporal
%                                 grid points
%       model.index (1:nc)      = pointer to assign parameters to cells
%                                   (nc gridsize in cells)
%       model.dz(1:nc)          = cell size (m)
%
%  structure timestep:
%       timestep.dt (1:nt)      = time step (s)
%       timestep.theta(1:nt-1) 	= time stepping control parameter,
%                      	     	  (.5 = Crank-Nicholson, 1=Backward Euler)
%       timestep.T0(1:nc+1)    	= initial temperatures
%
%  structure linsolve:
%       linsolve.solver         = 'direct','bicgstab','gmres'
%       if not 'direct':
%           linsolv.liniter         = number of iterations if not 'direct'
%           linsolv.lintol          = stopping tollerance  if not 'direct'
%           linsolv.precon          = preconditioner (ILU only)
%           linsolv.pretol          = tollerance for ILU
%           linsolv.prefreq         = frequency of preconditioner updates
%           linsolv.restart         = restart (for 'gmres')
%  structure nonlinear:
%       nonlinear.mean          = n/a/g/s calculation of effective properties
%                                   (fix, arithmetic, geometric, square-root)
%       nonlinear.maxiter       = maximal number of nonlinear Picard iterations
%       nonlinear.tol           = stopping tolerance for Picard iterations
%                                   (based on max deviation from last iteration)
%
%  structure freeze:
%       freeze.flag             = yes/no
%       freeze.Tf               = freezing temperature (default 0C)
%       freeze.w                = width of freezing interval (in K, ususlly 1)
%
%  out                          = debug output
%
% On output:
% J                             = Jacobian with repect to paeoclimate steps
% Tc                            = central value for temperature.
%
% vr, April 3, 2005

if nargin<7, out= 'no';end
if nargin<6,
    disp([' STOP. all parameters are obligatory. ']);
    return;
end

np          = inverse.np;
dp          = inverse.dp;
m           = inverse.mod;
maxit       = inverse.maxitdp;
tol         = inverse.toldp;

it          = model.it;
Ts          = model.Ts;
gt          = model.gt;
nt          = timestep.nt;ni=nt-1;
dt          = timestep.dt;
T0          = timestep.T0;

% calculate central value
timestep.out=1;
Tc=heat1dt(model,timestep,linsolve,nonlinear,freeze,out);

%      tic
ipoint=0;
timestep.out=0;
for i=1:np
    tinter=(find(it==i));
    if isempty(tinter),
        % set to zero if interval not present
        J (:,i)= zeros(nz,1);
    else
        % modify initial values (causality)
        istart      = min(tinter);
        timestep.T0 = Tc(:,istart);
        timestep.dt = dt(istart:ni);
        timestep.nt = length((istart:ni))+1;
        % perturb model
        m(i)        = m(i)+dp;
        model.Ts    = m(it(istart:nt))+gt;
        T=heat1dt(model,timestep,linsolve,nonlinear,freeze,out);
        % store to Jacobian J
        ipoint=ipoint+1;pointer(ipoint)=i;
        J (:,i)= (T(:)-Tc(:,nt))/dp;
        % reset model
        m(i)        = m(i)-dp;
    end
end
%      toc
% timestep.nt = nt;
% timestep.dt = dt;
% timestep.T0 = T0;
% model.Ts    = Ts;
% 
% % surface temperature
% model.gt = model.gt+dp;
% T=heat1dt(model,timestep,linsolve,nonlinear,freeze,out);
% J (:,np+1)= (T(:)-Tc(:,nt))/dp;
% model.gt = model.gt-dp;
% 
% % basal heat flow
% dp =1.e-3;
% model.qb = model.qb*(1.d0+dp);
% T=heat1dt(model,timestep,linsolve,nonlinear,freeze,out);
% J (:,np+2)= (T(:)-Tc(:,nt))/dp;
% model.gt = model.gt/(1.d0+dp);
