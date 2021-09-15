function [m_loc,T_loc,r_loc,theta_d,theta_m,r_norm,m_norm,gcv,upr]= ...
               solvereg(Jw,L0,L1,L2,m,m_apr,regpar,tol,maxiter,locpar)

mstruct(locpar);

A=[Jw;...
    sqrt(regpar(3))*L2; ...
    sqrt(regpar(2))*L1; ...
    sqrt(regpar(1))*L0 ];
rhs=[r_itr(:);...
    -sqrt(regpar(3))*L2*(m-m_apr);...
    -sqrt(regpar(2))*L1*(m-m_apr);...
    -sqrt(regpar(1))*L0*(m-m_apr)];
[delta_m,flag,relres] = lsqr(A,rhs,tol,maxiter);

m_loc= m + delta_m ; 
      
[T_loc]=forward_solve(m_loc,...
    k,kA,kB,h,r,c,p,ip,dz,qb,gts,...
    it,dt,T0,theta,maxitnl,tolnl,freeze,0);

r_loc=Tobs(:)-T_loc(id,length(dt)+1);
nd=length(r_loc);

% rms=sqrt(sum(res.*res)/nd));
%CALCULATE OBJECTIVE FUNCTION THETA
Wm          = sqrt(regpar(1))*L0 + sqrt(regpar(2))*L1 + sqrt(regpar(3))*L2;
theta_m     = norm(Wm*(m_loc-m_apr))^2;
Wd          = spdiags(1./Terr(:),0,nd,nd);
theta_d     = norm(Wd*r_loc)^2;


rms         = sqrt(norm(Wd*r_loc)^2/nd);
r_norm      = norm(Wd*r_loc);
m_norm1     = norm(m_loc);
m_norm2     = norm(Wm*(m_loc-m_apr));
m_norm      = m_norm2;

Jdag        = inv(Jw'*Jw+regpar(1)*(L0'*L0)+regpar(2)*(L1'*L1)+regpar(3)*(L2'*L2))*Jw';
M           = Jw*Jdag;
gcv         = nd*norm(Wd*r_loc)^2/trace(eye(size(M))-M)^2;
upr         = norm(Wd*r_loc).^2/nd -2*norm(Terr)^2 + 2*norm(Terr)^2*trace(Jw*Jdag)/nd;
                






end