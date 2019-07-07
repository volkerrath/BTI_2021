
clear all; close all; clc
addpath('/home/volker/Matlab//heat1d_toolbox/');
out= 1; 

% GENERATE SPATIAL MESH
zstart= 10; zend= 5000; nz= 251; type= 'logarithmic'; dir= 1; 
[z, dz]= set_mesh(zstart, zend, nz, type, dir, 0); 
% GENERATE TEMPORAL MESH
year2sec= 31557600; 
tstart= 100000* year2sec; tend= 10* year2sec; nt= 257; type= 'logarithmic'; dir= -1; 
[t, dt]= set_mesh(tstart, tend, nt, type, dir, 0); 
% GENERATE INVERSE MESH
gts= 0; 
[Ts, it]= paleo_haenel(t, -gts); 
%old: [Ts,it]=paleo_transylvania(t,-gts);


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
model.Ts= 10; 

ns= nt- 1; 
timestep.dt(1: nt- 1)= dt; 
timestep.theta(1: nt- 1)= 1.; 

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
Ti= heat1ds(model,linsolv,nonlinear, freeze, out); 
disp([ ' cpu time for initial temperatures:',num2str(cputime-time0),' s '])  

timestep.T0= Ti;
model.Ts= Ts; 
model.it= it; 

time0=cputime;
T=heat1dt(model,timestep,linsolv,nonlinear,freeze,out);
disp([ ' cpu time for transient model     :',num2str(cputime-time0),' s (',num2str(nt),' time steps)'])  
whos
opts = struct('Format','psc2','Color','rgb','Resolution',600);
figure;
plot([Ti,T],z, 'LineWidth',2);
set(gca,'YDir','reverse');grid on;ylim([0 5000]);xlim([0,100]);
xlabel('T (^\circ C)'):ylabel ('z (m)');
grid on;
legend('steady-state','transient','Location','northeast');
filename=strcat(['head1d_test.ps']);
exportfig(gcf,filename,opts)
% close(gcf)

