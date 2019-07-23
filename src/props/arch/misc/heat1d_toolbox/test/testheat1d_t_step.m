% TESTHEAT1D
%
% Script um HEAT1DN zu testen
% v.r. 8/01                    last change: Aug.16,2001 
clear all;close all;

disp('Test transient heat conduction')

z=[0:20:2000];      %[0:200:2000]';      			
dz=diff(z);
nz=length(z);
gt=  10.  ; % 10; %
qb= 60e-3;               				
pc =  [1 2 3 4 5];  %[ 2 2 2 2 2];
ph =  0. *ones(1,5)*1.e-6;
pr =  2.3*ones(1,5)*1.e6;
np=length(pc);
ic=[2 2 2 2 2  2 2 2 2 2  2 2 2 2 2   2 2 2 2 2  2 2 2 2 2 ...
    4 4 4 4 4  4 4 4 4 4  4 4 4 4 4   4 4 4 4 4  4 4 4 4 4 ...
    2 2 2 2 2  2 2 2 2 2  2 2 2 2 2   2 2 2 2 2  2 2 2 2 2 ...
    3 3 3 3 3  3 3 3 3 3  3 3 3 3 3   3 3 3 3 3  3 3 3 3 3 ...
]'; 

% calculate inital model 
par=[pc ph qb gt];

T0=heat1dns(par,ic,dz); 

%figure
%plot(T0,-z)
%title('Temperature (FD)');grid on;
%xlabel('Temperature [\circ C]');ylabel('Depth [m]');



clear par;
par=[pc ph pr qb];

year=31557600;
t=[0:100:30000];nt=length(t);
t=t*year;dt=diff(t);
%  set step function
tt=10000*year;
[Ts]=gt+paleo_step(t,5,tt);


beta=0.5;

%figure
%plot(t/year,Ts);
%title('Surface Temperature History');grid on;
%xlabel('Time [a]');ylabel('Temperature [\circ C]');


tic
T=heat1dnt(par,ic,dz,dt,T0,Ts,beta);
toc


for i=1:nt
 Tr(:,i)=T(:,i)-T0;
end

figure;plot(Tr(:,[100:10:200]),-z)
title(' Reduced Temperature (FD)');grid on;
xlabel('Temperature [K]');ylabel('Depth [m]');
legend('t =  0 ka','  =  1 ka','  =  2 ka','  =  3 ka','  =  4 ka', ...
       '  =  5 ka','  =  6 ka','  =  7 ka','  =  8 ka','  =  9 ka', ...
       '  = 10 ka')  

figure;
imagesc(T); colorbar;
title('Temperature (FD)');grid on;
xlabel('Time [steps]');ylabel('Depth [cells]');
figure;
imagesc(Tr); colorbar;
title('Residual Temperature (FD)[\circ C]');grid on;
xlabel('Time [steps]');ylabel('Depth [cells]');
