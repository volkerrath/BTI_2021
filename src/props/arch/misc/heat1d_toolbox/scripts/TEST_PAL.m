% GENERATE TEMPORAL MESH
clear all;close all;clc

opts = struct('Format','psc2','Color','rgb','Resolution',600);
disp(['Nonlinear Paleoclimate Modeling'])
disp(['- Movies -'])
disp(['V. R., ',date] )
disp([' ']);

% GENERATE SPATIAL MESH 
zstart=10;zend=3000;nz=601;
methz='linear';dirz=1;
%disp([' ']);disp([' ...set up ' methz ' spatial mesh ']); 
[z,dz]=set_mesh(zstart,zend,nz,methz,dirz,0);
year2sec=31557600;
tstart=200000*year2sec;tend=100*year2sec;nt=801;
metht='logarithmic';dirt=-1;
%disp([' ']);disp([' ...set up ' metht ' temporal mesh ']); 
[t,dt]=set_mesh(tstart,tend,nt,metht,dirt,0);

% SETUP PALEOCLIMATE
GTemp=[-5 -10 0];GTime=[-100000 -10000]*year2sec;
[Ts,T]=paleo_boxcar_smooth(t,GTemp,GTime,8);

% figure; 
% filename=strcat('Paleoclimate.ps');
% plot(-t/year2sec,Ts,'-b', 'LineWidth',3); grid on;hold on;
% plot(-t/year2sec,T,'--r','LineWidth',3);
% xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
% ylim([-15 5]); xlim([100 2e5]);
% set(gca,'XScale','log','XDir','reverse')
% title(['Paleoclimate'],'FontSize',14)
% grid on;
% exportfig(gcf,filename,opts)

% PARAMETER FOR FWD CALCULATIONS
theta=1.*ones(nt,1);
% theta(1:10)=1.;
% NONLINEAR ITERATION PARAMETERS
maxitnl=3;
tolnl=0.0001;



% READ UNITS AND PROPERTIES FOR BOREHOLE
filename=strcat('PALMovie.prp');
[unit] = get_props(filename);units=length(unit);
disp([' ...properties from    : ' filename]);

% READ MODEL GEOMETRY FOR BOREHOLE
filename=strcat('PALMovie1.mod');
[model]=get_model(filename,z);
disp([' ...model geometry from: ' filename]);

% SETUP INITIAL VALUES FOR TEMPERATURE AT TSTART
disp([' ']); disp([' ...setup initial values ' ]); 
kl=[unit.k];hl=[unit.h];kAl=[unit.kA];kBl=[unit.kB];porl=[unit.p];
cpml=[unit.c];rhoml=[unit.r];
qb=model.qb;gt=model.gt;ip=model.ip;
gt1=gt+Ts(1);
time0=cputime;
%Ti=heat1dns(kl, kAl, kBl,hl,porl,qb,gt1,ip,dz,maxitnl,tolnl,'no');
%init.T1=Ti;
Ti=heat1dns(kl, kAl, kBl,hl,porl,qb,gt1,ip,dz,maxitnl,tolnl,'yes');
init.T2=Ti;
disp([ ' cpu time for initial model <yes> :',num2str(cputime-time0),' s '])  


disp([' ']); disp([' ...calculate model ' ]); 
Ts=Ts+gt;
% FORWARD MODEL
time0=cputime;
% Tini=init.T1;
% [Tcalc1,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]= ...
%     heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
%     ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,'no');
% disp([ ' cpu time for model <no> :',num2str(cputime-time0),' s '])  

time0=cputime;
Tini=init.T2;
[Tcalc2,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]= ...
    heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
    ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,'yes');
disp([ ' cpu time for model <yes> :',num2str(cputime-time0),' s '])

save PALMovie1
% 
% % READ MODEL GEOMETRY FOR BOREHOLE
% filename=strcat('PALMovie2.mod');
% [model]=get_model(filename,z);
% disp([' ...model geometry from: ' filename]);
% 
% % SETUP INITIAL VALUES FOR TEMPERATURE AT TSTART
% disp([' ']); disp([' ...setup initial values ' ]); 
% kl=[unit.k];hl=[unit.h];kAl=[unit.kA];kBl=[unit.kB];porl=[unit.p];
% cpml=[unit.c];rhoml=[unit.r];
% qb=model.qb;gt=model.gt;ip=model.ip;
% gt1=gt+Ts(1);
% Ti=heat1dns(kl, kAl, kBl,hl,porl,qb,gt1,ip,dz,maxitnl,tolnl,'no');
% init.T1=Ti;
% Ti=heat1dns(kl, kAl, kBl,hl,porl,qb,gt1,ip,dz,maxitnl,tolnl,'yes');
% init.T2=Ti;
% 
% disp([' ']); disp([' ...calculate model ' ]); 
% Ts=Ts+gt;
% % FORWARD MODEL
% time0=cputime;
% Tini=init.T1;
% [Tcalc1,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]= ...
%     heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
%     ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,'no');
% disp([ ' cpu time for model <no> :',num2str(cputime-time0),' s '])  
% 
% time0=cputime;
% Tini=init.T2;
% [Tcalc2,zout,tout,N,k_eff,rc_eff,ipor,lheat,rci]= ...
%     heat1dnt(kl,kAl,kBl,hl,rhoml,cpml,porl,qb,...
%     ip,dz,dt,Tini,Ts,theta,maxitnl,tolnl,'yes');
% disp([ ' cpu time for model <yes> :',num2str(cputime-time0),' s '])
% 
% save PALMovie2
