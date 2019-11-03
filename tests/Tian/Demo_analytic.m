clear all;
close all;
clc
srcpath='../../src/';
addpath([srcpath]);

y2s=3600*24*365.25;s2y=1./y2s;
truecolor=[0.9167    0.9167    0.9167]*.9;
% set up meshes and GSTH

disp([' ']);disp([' ...setup spatial grid'])
zstart=10; zend=6000; nz=601; type= 'logarithmic'; dir= 1;
[z, dz]= set_mesh(zstart, zend, nz, type, dir,0);

disp([ ' ']);disp([ ' ...setup temporal grid'])
tstart= 80000*y2s; tend= 10*y2s; nt=2001; type= 'logarithmic'; dir= -1;
[t, dt]= set_mesh(tstart, tend, nt, type, dir,0);

disp([ ' ']);disp([ ' ...setup true GSTH'])
L=3;
GTemp=[  -2       -6         1        0      1     1.5]-1.;
GTime=[    -60000     -12000     -440     -150    -50   -1 ]*y2s;
[Ts]= set_gst(t,GTemp,GTime,L,0);
Ts_true=Ts;

% call 1D  forward model
out= 1;
debug=0;
% time axis
tlog=  2000.;
refyr= 2009;
% physical parameters
k=2.5;
rho   = 2500;
cp    = 1000;
kappa = k/(rho*cp);
Qb   = 0.05;
T0   = 0.;
gtime=flipud(-t');
gsth=flipud(Ts);


parvec=[gsth;T0;Qb];
[T,Q]=heat1dana(parvec,gtime,kappa,k,z',tlog,refyr,out);

T0=T.val;z=T.z;
save('LGM_all_AN.mat'); 
