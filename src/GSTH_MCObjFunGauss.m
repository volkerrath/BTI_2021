function [s,r,c]=GSTH_ObjFunInit1(m,data)


% LOAD MODEL PARAMETERS, OBSERVATIONS, INITIAL VLUES
sitepar=data.sitepar;
mstruct(sitepar);dz=diff(z);nz=length(z);dt=diff(t);nt=length(t);
numpar=data.numpar;
mstruct(numpar);
T0=data.initial;
gsthpar=data.gsthpar;
it=gsthpar.pt;
% MODIFY MODEL IF NECESSARY
mactive=data.mactive;
nm=length(m);
GST=m(1:nm-2)+gts; 
%mactive

if mactive(nm-1),QB=-m(nm-1)*1e-3; else QB=qb;end 
if mactive(nm-0),H=m(nm-0)*1.e-6;H=H*ones(size(h)); else H=h;end 

out=0;POM=-4;
T0=heat1dns(k, kA, kB,H,p,QB,GST(1)+POM,dz,ip,maxitnl,tolnl,freeze,out);

% RUN MODEL
[Tcalc,G]=heat1dnt(k,kA,kB,H,rc,p,QB,...
    dz,ip,dt,it,GST,T0,theta,maxitnl,tolnl,freeze,1);
% CALCULATE RESIDUAL
nd=length(Tobs);
%Wd=spdiags(1./err,0,nd,nd);
%res=Wd*(Tobs-Tcalc(id,nt));
c=Tcalc(id,nt);
r=(Tobs-Tcalc(id,nt));
s=sum(r.^2);
% disp([num2str(min(r)) ' - ' num2str(max(r))]);
% save('MCMC','k','kA','kB','H','rc','p','QB','dz','ip','dt','it','GST','T0','theta','c','G')
% 
% linwdt=2;
% fontwg='normal';
% fontsz=16;
% plot(c,z(id),'.b','LineWidth',1); hold on
% plot(Tobs,zobs,'-r','LineWidth',1); hold on
% xlim([0 50]);
% ylim([0 2750]);
% set(gca,'ydir','rev','FontSize',fontsz,'FontWeight',fontwg)
% xlabel('T (C)','FontSize',fontsz,'FontWeight',fontwg)
% ylabel('depth (m)','FontSize',fontsz,'FontWeight',fontwg)
% % title(strcat([site, ' temperatures']),'FontSize',fontsz,'FontWeight',fontwg);
% S1=strcat(['orig']);
% S2=strcat(['mod']);
% legend(S1,S2,'location', 'southwest')
% grid on
% saveas(gcf,'test','png')