function [Tcalc,ipor] = forward_solve1(m,k,kA,kB,h,p,c,r,rc,ip,dz,qb,gt,...
                                 it,dt,theta,maxitnl,tolnl,freeze,out)
% FORWARD_SOLVE solves forward problem
%   Detailed explanation goes here
GST=m+gt;POM=GST(1);
T0=heat1dns(k, kA, kB,h,p,qb,POM,dz,ip,maxitnl,tolnl,freeze,out);
[Tcalc,ipor]=heat1dnt(k,kA,kB,h,rc,p,qb,...
    dz,ip,dt,it,GST,T0,theta,maxitnl,tolnl,freeze,1);
end

