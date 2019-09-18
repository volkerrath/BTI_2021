% TIHUNVINV_D performs joint inversion in parameter space
% V.R. : Oct.16, 2001 
close all;clear all

% GENERATE WEAKLY GRADED MESH 
dz=logspace(1,1.5,161);z=[0,cumsum(dz)];nz=length(z);

% GENERATE MODEL PARAMETER 
[ip,np]   =setmodel('./testmodel_synthetic',z);  
ic=ip;
properties_molasse;

% 
par=[pc ph qb gt];
T0=heat1dns(par,ic,dz);

% 
maxnl=25;tolnl=0.01;
par=[pc ph qb gt,pcA,pcB];
T=heat1dns_nl(np,par,ic,dz,maxnl,tolnl); 

subplot(1,2,1)
plot(T0 ,-z,'-b');grid on;hold on;
plot(T  ,-z,'-r');grid on;hold on;
title(strcat(['Nonlinear \lambda (t)']));
xlabel('Temperature [\circ C]');
ylabel('Depth [m]');legend('constant','nonlinear');

subplot(1,2,2)
plot(T-T0  ,-z,'-r');grid on;hold on;
xlabel('\Delta T [K]'); ylabel('Depth [m]');
title('Residuals');


