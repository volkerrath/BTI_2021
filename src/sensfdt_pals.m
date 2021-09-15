function[J,pointer,Tc]=sensfdt_pal(k,kA,kB,h,r,c,p,qb,...
    dz,ip,dt,it,gst,T0,theta,maxitnl,tolnl,dp,freeze,out);
% SENSFDT_PAL calculates transient Jacobians with respect to paleoclimate
%
% [J]=sensfdt(kl,kAl,kBl,hl,rhoml,cpml,porl,rhocpf,rhoi,cpi,qb,...
%     dz,ip,dt,it,GST,T0,theta,...
%     maxitnl,tolnl,dp,maxitdp,toldp,freeze,out)
% caclulates Jacobians of recent borehole temperatures wrt paleoclimate 
% values given in gst(it), and heat flow density qb.
%
% V. R., March 12, 2003
tic;

nz=length(dz)+1;
nt=length(dt)+1;
ngst=length(gst);
J = 0*ones(nz,ngst+1);
dTref=1;
% CALCULATE CENTRAL VALUE
if isempty(T0)
    T0=heat1dns(k,kA,kB,h,r,p,qb,...
        gst(1),dz,ip,maxitnl,tolnl,freeze,1);
end

Tc=heat1dnt(k,kA,kB,h,r,c,p,qb,...
    dz,ip,dt,it,gst,T0,theta,maxitnl,tolnl,freeze,1);
ipoint=0;
for i=1:ngst
    if isempty(find(it==i)),
        %disp([' GTH value no is not in given time interval']);
        J (:,i)= zeros(nz,1);
    else
        %disp([' GTH value no ',num2str(i),' = ',num2str(gst(i))]);
        %disp([' cells: ',num2str(find(it==i))]);
        istart=min(find(it==i));
        T0i =Tc(:,istart);
        del=dp*dTref;
        gst(i)=gst(i)+del;
        iti=it(istart:nt);
        Tsi  = gst(it(istart:nt));
        ni=length(Tsi);
        dti =  dt(istart:nt-1);
        T=heat1dnt(k,kA,kB,h,r,c,p,qb,...
            dz,ip,dti,iti,gst,T0i,theta,maxitnl,tolnl,freeze,0);
        ipoint=ipoint+1;
        J (:,ipoint)= (T(:)-Tc(:,nt))/del;
        gst(i)=gst(i)-del;
    end
end
% CALCULATE SENSITVITY WRT QB
del=dp*qb;
qbi=qb+del;
T=heat1dnt(k,kA,kB,h,r,c,p,qbi,...
    dz,ip,dti,iti,gst,T0,theta,maxitnl,tolnl,freeze,0);
ipoint=ipoint+1;
J (:,ipoint)= (T(:)-Tc(:,nt))/del;


disp([ ' cpu time for FD transient jacobian    :',num2str(toc),' s '])
%Jfd=J;
%save fd1 fdtime Jfd

end
