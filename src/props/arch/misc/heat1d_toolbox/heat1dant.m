function T=heat1dat(par,nu,index,dz,dt,T0,Ts)
% HEAT1DNT solves analytically time-dependent heat equation (N&B,89)
%
% T = heat1dnt(par,nu,index,dz,dt,T0,Ts,methods) calculates 
% temperatures for e givern set of timesteps, given a model for thermal 
% conductivity, heat production and rho*c. 
% transient heat conduction with time-dependent source term is assumed. 
% Input :
% par (     1:  nu) = lambda 
% par (  nu+1:2*nu) = volumetric heat production
% par (2*nu+1:3*nu) = rho c
% par (3*nu+1)      = basal heat flow 
% nu                = number of specified units
% index (1:nc)      = pointer to assign parameters to cells 
%                    (nc gridsize in cells)
% dz(1:nc)          = cell size (m)
% dt(1:nt-1)        = time step (s)
% T0(1:nc+1)        = initial temperatures
% Ts(1:nt+1)        = time-dependent boundary temperatures 
%                     at the top (e.g., paleoclimate)
% methods(1:nt-1)      = methods(1:nt) is time stepping control parameter,
%                      .5 =Crank-Nicholson, 1=Bacward Euler
%
% Output:
% T(1:nc+1,1:nt )  = temperatures at given time steps dt
% 
% V. R., Oct. 27, 2001
display('heat1dat is not yet implemented')
