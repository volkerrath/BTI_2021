function[J,Tc]=sensfds_pet(jacpet,mod,timestep,nonlinear,freeze,out)
% JACFDT_PAL calculates transient Jacobians with respect to paleoclimate. 
%     [J]=jacfdt_pal(jacpet,mod,timestep,nonlinear,freeze,out) 
% caclulates Jacobians of recent (index nt+1) borehole temperatures with 
% respect to paleoclimate values given in pt(it) by forward differencing.
%   
%  structure jacpet:
%       jacpet.pt          = values of gsth steps
%       jacpet.it	   = pointer assiciating pt values to temporal
%                            grid points   	
%       jacpet.maxitdp     = maximal number of FD refinements for jacobian
%			     calculations
%       jacpet.toldp       = stopping tolerance for jacobian calculations
%       jacpet.flag        = controls the treatment of parameters: 
%                            (0/no differentiation, 1/linear, 2/logarithimic) 
%%  
%  structure mod:
%       mod.k (1:nu)        = lambda
%       mod.kA, kB(1:nu)	= thermal conductivity coefficient A and B
%       mod.h(1:nu)		    = volumetric heat production
%       mod.qb			    = basal heat flow
%       mod.Ts(1:nt+1)    	= time-dependent boundary temperatures 
%                     		   at the top (e.g., paleoclimate)
%       mod.index (1:nc)    = pointer to assign parameters to cells
%                   	   	  (nc gridsize in cells)
%       mod.dz(1:nc)      	= cell size (m)
%
%  structure timestep:
%       timestep.dt (1:nt)      = time step (s)
%       timestep.theta(1:nt-1) 	= time stepping control parameter,
%                      		      (.5 = Crank-Nicholson, 1=Backward Euler)
%       timestep.T0(1:nc+1)    	= initial temperatures

%  structure nonlinear:
%       nonlinear.mean      = n/a/g/s calculation of effective properties
%                             (fix, arithmetic, geometric, square-root)
%       nonlinear.maxiter   = maximal number of nonlinear Picard iterations
%       nonlinear.tol       = stopping tolerance for Picard iterations
%                             (based on max deviation from last iteration)
%  
%  structure freeze:
%       freeze.flag         = yes/no 
%       freeze.Tf           = freezing temperature (default 0C)
%       freeze.w            = width of freezing interval (in K, ususlly 1)
%
%  out                      = debug output
%
% On output:
% J           		    = Jacobian with repect to 
% Tc                        = central value for temperature.
%
% vr, April 3, 2005
 
if nargin<5, out= 'no';end
if nargin<4, disp([' STOP. all parameters are obligatory. ']); 
break; end

dz      = timestep.dz;nt=length(dz)+1;
dt	= timestep.dt;nt=length(dt)+1;
pt	= jacpet.pt;np=length(pt);
dp	= jacpet.dp;
flag    = jacpet.flag;
nunits  = length(mod.kl);


Tc = heat1ds(mod,nonlinear,freeze,out); 
for i=1:nunits 
      dpi=dp;
      for k=1:maxitdp
         mp=m;
	 case()
	 % forward perturbation
         mp(ipar)=m(ipar)+dpi;
	     pp=m2p(mp,flags(2,:));
         Tp=heat1dns(pp,nu,ic,dz,maxitnl,tolnl);
         dTp=(Tp-Tc);
	 % backward perturbation
         mp(ipar)=m(ipar)-dpi;
	     pp=m2p(mp,flags(2,:));
         Tm=heat1dns(pp,nu,ic,dz,maxitnl,tolnl);
         dTm=(Tc-Tm);
	 % check convergence
	     test=max(abs(dTp-dTm)/Tc);
         if  test < jacpet.toldp break; end
         dpi=dpi/10;
       end
       J(:,    i)=0.5*(Tp-Tm)/dpi;
%       fprintf(' parameter no. %3i :  %g  at delta= %g   \n',...
%               i,test,dpi);
end

