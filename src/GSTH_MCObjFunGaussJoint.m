function [s,r,c]=GSTH_ObjFunJoint(m,info)


% LOAD MODEL PARAMETERS & OBSERVATIONS
sitepar=info.sitepar;
mstruct(sitepar);
dz=diff(z);nz=length(z);
dt=diff(t);nt=length(t);
numpar=info.numpar;
mstruct(numpar);
mactive=info.mactive;


% MODIFY MODEL
nm   =length(m); 
ncnd =length(k);
ngst = nm-(ncnd+2); 
GST=m(1:ngst)+gts; 

if mactive(ngst+1:ngst+ncnd),   KB=m(ngst+1:ngst+ncnd); else KB=k;end
if mactive(ngst+ncnd+1),        QB=-m(nm-1)*1e-3; else QB=qb;end 
if mactive(ngst+ncnd+2),        HB=m(nm-0)*1.e-6;H=HB*ones(size(h)); else H=h;end 
    
if isempty(Tinit)
    POM=GST(1);
    T0=heat1dns(KB, kA, kB,H,p,QB,POM,dz,ip,maxitnl,tolnl,freeze,0);
else
    T0=Tinit;
end


[Tcalc]=heat1dnt(KB,kA,kB,H,rc,p,QB,...
    dz,ip,dt,it,GST,T0,theta,maxitnl,tolnl,freeze,1);

% CALCULATE RESIDUAL
nd=length(Tobs);
%Wd=spdiags(1./err,0,nd,nd);
%res=Wd*(Tobs-Tcalc(id,nt));
c=Tcalc(id,nt);
r=(Tobs-Tcalc(id,nt));
s=sum(r.^2);

