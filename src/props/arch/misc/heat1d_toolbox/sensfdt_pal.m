function[J,pointer,tc]=sensfdt_pal_pf(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
    index,dz,dt,T0,pt,it,theta,maxitnl,tolnl,dp,maxitdp,toldp,freeze,out)
% SENSFDT_ calculates transient Jacobians with respect to paleoclimate 
%
% [J]=sensfdt(kl,kAl,kBl,hl,rhoml,cpml,porl,rhof,cpf,rhoi,cpi,qb,...
%     index,dz,dt,T0,Ts,,pt,it,theta,...
%     maxitnl,tolnl,dp,maxitdp,toldp,freeze,out)
% caclulates Jacobians of recent (index nt+1) borehole temperatures with respect to 
% paleoclimate values given in pt(it).
%
%
% 
% V. R., March 12, 2003 


nz=length(dz)+1;
nt=length(dt)+1;
npt=length(pt);

% calculate central value 
Ts=pt(it);
Tc=heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
    index,dz,dt,T0,Ts,theta,maxitnl,tolnl,freeze,-1);
ipoint=0;
for i=1:npt
    if isempty(find(it==i)), 
        %disp([' GTH value no is not in given time interval']);
        J (:,i)= zeros(nz,1);
    else
        %disp([' GTH value no ',num2str(i),' = ',num2str(pt(i))]);
        %disp([' cells: ',num2str(find(it==i))]);
        istart=min(find(it==i));  
        T0i =Tc(:,istart);
        del=dp;
        pt(i)=pt(i)+del;
        Tsi  = pt(it(istart:nt));ni=length(Tsi);
        dti =  dt(istart:nt-1); 
        T=heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
            index,dz,dti,T0i,Tsi,theta,maxitnl,tolnl,freeze,0);
        
        ipoint=ipoint+1;
        pointer(ipoint)=i;
        J (:,i)= (T(:)-Tc(:,nt))/del;
        pt(i)=pt(i)-del;
    end
end
%      toc

