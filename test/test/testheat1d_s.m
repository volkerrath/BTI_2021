% TESTHEAT1D
%
% Script um HEAT1DN zu testen
% v.r. 8/01                    last change: Aug.15,2001 
clear all;close all;

disp('Test stationary heat conduction')

z=[0:10:5000];      %[0:200:2000]';      			
dz=diff(z);
nz=length(z);
gt=  10.  ; % 10; %
qb= 60e-3;               				
pc =  [ 2 ];
ph =  2.5e-6;
np=length(pc);
ic=[ ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
    1 1 1  1 1  1 1 1  1 1  1 1 1  1 1   1 1 1  1 1  1 1 1  1 1 ...
]'; 


figure

par=[pc ph qb gt];
Tn=heat1dns(par,ic,dz); 
Ta=heat1dai(par,ic,dz) ;

subplot(2,2,1)
plot([Ta' Tn] ,-z);grid on;
title(strcat(['Heat production ',num2str(ph*1e6),' \muW/m^3']));
%xlabel('Temperature [\circ C]');
ylabel('Depth [m]');
legend(' analytic','numeric');
subplot(2,2,2)
plot(Tn-Ta',-z);grid on;
%xlabel('T_{num} - T_{ana} [K]');
%ylabel('Depth [m]');
%title('Heat production');

par=[pc 0. qb gt];
Tn=heat1dns(par,ic,dz); 
Ta=heat1dai(par,ic,dz) ;
subplot(2,2,3)
plot([Ta' Tn] ,-z);grid on;
title('No  Heat production');
xlabel('Temperature [\circ C]');ylabel('Depth [m]');
legend(' analytic','numeric');
subplot(2,2,4)
plot(Tn-Ta',-z);grid on;
xlabel('T_{num} - T_{ana} [K]');
%ylabel('Depth [m]');
%title('No Heat production');
