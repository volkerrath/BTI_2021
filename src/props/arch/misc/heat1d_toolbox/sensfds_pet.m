function[J,Tc]=sensfds_pet(m,nu,ic,dz,...
                   maxitnl,tolnl,flags,dp,maxitdp,toldp)
% SENSFDS calculates stationary nonlinear Jacobian 
% with respect to layer parameters and boundary values
% by brute-force divided differences 
%
% [J]=sensfds_nl(par,nu,ic,dz,maxitnl,tolnl,flags,dp,maxitdp,toldp)
%
% On input: 
% m               set of parameter, possibly logarithmic 
% nu              number of specified units,
% dz ,ic          describe model geometry.
% maxitnl,tolnl   for nonlinear iteration
% dp              relative parameter perturbation
% maxitdp, toldp  number of refinements, required tolerance. 
% flag            controls the treatment of parameters: 
%                 flag(1,:) determines wether the parameter is to be
%                 differntiated 
%                 flag(2,:) logarithimic or linear differntiation. 
%
% On output: 
% J               Jacobian. 
% Tc              central value for temperature.
%   
%
% V. R.  Nov. 5, 2002 


if nargin <11, out=0;end
if nargin< 8, error('SENSFDS_PET: no differentiation control'); end

p=m2p(m,flags(2,:));
Tc = heat1dns(p,nu,ic,dz,maxitnl,tolnl); 

list=find(flags(1,:)==1);nl=length(list);
for i=1:nl 
      dpi=dp;ipar=list(i);
      for k=1:maxitdp
         mp=m;
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
	     test=max(abs(dTp-dTm));
         if  test < toldp break; end
         dpi=dpi/10;
       end
       J(:,    i)=0.5*(Tp-Tm)/dpi;
%       fprintf(' parameter no. %3i :  %g  at delta= %g   \n',...
%               i,test,dpi);
end

