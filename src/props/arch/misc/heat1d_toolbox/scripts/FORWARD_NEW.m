
clear all; close all; clc
% addpath('/home/volker/Matlab//heat1d_toolbox/');
out= 1;
debug=1;

% GENERATE SPATIAL MESH
zstart= 10; zend= 5000; nz= 251; type= 'logarithmic'; dir= 1;
[z, dz]= set_mesh(zstart, zend, nz, type, dir, debug);
% GENERATE TEMPORAL MESH
year2sec= 31557600;
tstart= 100000* year2sec; tend= 10* year2sec; nt= 257; type= 'logarithmic'; dir= -1;
[t, dt]= set_mesh(tstart, tend, nt, type, dir, debug);
% GENERATE INVERSE MESH
GTemp=[-10       -20      -2  2 ];
GTime=[-140000 -12000 -400 -10 ]*year2sec;
L=5;
[Ts]= paleo_boxcar_smooth(t,GTemp,GTime,L,debug);



% SETUP STRUCTURES
nu= 1; nc= nz- 1;
model.k(1: nu)= 3.;
model.kA(1: nu)= 0.;
model.kB(1: nu)= 0.;
model.p(1: nu)= 0.1;
model.h(1: nu)= 1.e-6;
model.r(1: nu)= 2300;
model.c(1: nu)= 1000;
model.dz= dz;
model.ip(1: nc)= 1;
model.qb= .040;
model.Ts= Ts(nt);
model.it= [1:nt];

ns= nt- 1;
timestep.dt(1: nt- 1)= dt;
timestep.theta(1: nt- 1)= 1.;
timestep.out    	= 0;

linsolv.solver='direkt';

nonlinear.mean= 'g';
nonlinear.maxiter= 2;
nonlinear.tol= 1.e-4;
nonlinear.relax= 1.;

freeze.flag= 'y';
freeze.Lh= 333600;
freeze.Tf= 0.;
freeze.w= 1.;


time0=cputime;
Ti= heat1ds(model,linsolv,nonlinear,freeze,debug);
disp([ ' cpu time for initial temperatures :',num2str(cputime-time0),' s '])

timestep.T0= Ti;
model.Ts= Ts;

time0=cputime;
%profile on
T=heat1dt(model,timestep,linsolv,nonlinear,freeze,debug);
%profile viewer

disp([ ' cpu time for transient model      :',...
    num2str(cputime-time0),' s (',num2str(nt),' time steps)'])

opts = struct('Format','psc2','Color','rgb','Resolution',600);
figure;
plot(-t/year2sec',Ts, 'LineWidth',2,'Color','r');
set(gca,'XScale','log','XDir','reverse')
grid on;ylim([-20 5]);xlim([10,100000]);
xlabel('time b. p.(a)','FontSize',14);ylabel('\Delta T','FontSize',14);
grid on;
filename=strcat(['head1d_test1.ps']);
exportfig(gcf,filename,opts)
% close(gcf)

figure;
plot([Ti,T],z, 'LineWidth',2);
set(gca,'YDir','reverse');grid on;ylim([0 5000]);xlim([-25,100]);
xlabel('T (^\circ C)');ylabel ('z (m)');
grid on;
legend('steady-state','transient','Location','northeast');
filename=strcat(['head1d_test2.ps']);
exportfig(gcf,filename,opts)
% close(gcf)

