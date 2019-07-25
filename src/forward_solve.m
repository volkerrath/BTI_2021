function [Tcalc,dT,Q,K] = forward_solven(m,k,kA,kB,h,r,c,rc,p,ip,dz,qb,gts,...
    it,dt,T0,theta,maxitnl,tolnl,freeze,out)
% FORWARD_SOLVE solves forward problem
%   Detailed explanation goes here
GST=m+gts;POM=GST(1);

if isempty(T0)
    T0=heat1dns(k, kA, kB,h,r,p,qb,POM,dz,ip,maxitnl,tolnl,freeze,out);
end

[Tcalc,dT,Q,K]=heat1dnt(k,kA,kB,h,r,c,rc,p,qb,...
    dz,ip,dt,it,GST,T0,theta,maxitnl,tolnl,freeze,1);

end

