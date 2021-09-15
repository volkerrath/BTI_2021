function[J,Tc,Ti]=sensfdt_palp(k,kA,kB,h,r,c,p,qb,...
    dz,ip,dt,it,gst,T0,theta,maxitnl,tolnl,dp,freeze,out)
% SENSFDT_PAL calculates transient Jacobians with respect to paleoclimate
%
% [J]=sensfdt(kl,kAl,kBl,hl,rhoml,cpml,porl,rhocpf,rhoi,cpi,qb,...
%     dz,ip,dt,it,GST,T0,theta,...
%     maxitnl,tolnl,dp,maxitdp,toldp,freeze,out)
% caclulates Jacobians of recent borehole temperatures wrt paleoclimate
% values given in gst(it), and heat flow density qb.
%
% v. r., jan 1, 2013

tic;
nz=length(dz)+1;nt=length(dt)+1;ngst=length(gst);
J = 0*ones(nz,ngst+1);
dTref=1;del=dp*dTref;
% CALCULATE CENTRAL VALUE

if isempty(T0)
    T0=heat1dns(k,kA,kB,h,r,p,qb,...
        gst(1),dz,ip,maxitnl,tolnl,freeze,1);
end

Tc=heat1dnt(k,kA,kB,h,r,c,p,qb,...
    dz,ip,dt,it,gst,T0,theta,maxitnl,tolnl,freeze,1);

parfor icol=1:ngst
    istart=find(it==icol, 1 );
    %disp([num2str(icol),' ',num2str(istart)])
    dgst=gst;deli=del;
    dgst(icol)=dgst(icol)+deli;
    its=it(istart:nt);
    dts=dt(istart:nt-1);
    Ts=Tc(:,istart);
    %     %     if icol==1,
    %     %          T0=heat1dns(k,kA,kB,h,p,qb,...
    %     %              dgst(1),dz,ip,maxitnl,tolnl,freeze,1);
    %     %     end
    T=heat1dnt(k,kA,kB,h,r,c,p,qb,dz,ip,...
        dts,its,dgst,Ts,...
        theta,maxitnl,tolnl,freeze,0);
    %
    %     %find(T(:)-Tc(:,nt))==0)
    J(:,icol)= (T(:)-Tc(:,nt))/deli;
end

disp([ ' cpu time for parallel FD transient jacobian    :',num2str(toc),' s ']);
if out>0,
    Jfdp=J;
    filename=strcat(['pal_fdp_', num2str(dp),'_',num2str(iter),'.mat']);
    if mod(iter,out)==0, save(filename,'Jfdp'); end
end

end

