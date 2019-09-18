% TESTHEAT1D
%
% Script um HEAT1DN zu testen
% v.r. 8/01                    last change: sept. 4,2001 
clear all;close all;

disp('Test transient heat conduction with heat production')

% physical parameters
z=[0:10:5000];       			
dz=diff(z);
nz=length(z);
gt= 10.  ; 
qb= 0e-3;               				
pc =  [2];  
ph =  2.e-6 ;
pr =  2.e6;
np=length(pc);
ic=ones(1,length(dz)); 
par=[pc ph pr qb];



% initial values
T0= gt*ones(nz,1);

% time steps &time dependent surface temperatures
year=31557600;
t=[0:100:30000];nt=length(t);
t=t*year;dt=diff(t);
[Ts]=zeros(1,nt);Ts(1)=gt;


beta=0.5*ones(1,nt);
beta(1:10)=1.;

T=heat1dnt(par,ic,dz,dt,T0,Ts,beta);
par=[pc 0. pr qb];
T2=heat1dnt(par,ic,dz,dt,T0,Ts,beta);





% numerical solution 
figure
plot(T(2:10,1:50)');
title(' Temperature node 2-10');grid on;
ylabel('Temperature [K]');xlabel('time step');
legend(' 2',' 3',' 4',' 5',' 6',' 7',' 8',' 9','10',4)  ;
text(30,1, strcat('Theta = ',num2str(beta))) ;  

% analytical solution 
% initial values

a=gt;
b=0.;%03;
T0=a+b*z';

kappa=pc/pr;
AK=ph/pc;

Tc(:,1)=T0;
for i=2:nt
   term1=(a+kappa*AK*t(i)+0.5*AK*z'.^2).*erf(z'.*0.5./sqrt(kappa*t(i)));
   term2= AK*z'*sqrt(kappa*t(i)/pi).*exp(-z'.^2/(4*kappa*t(i)));
   term3= - 0.5*AK*z'.^2;
   term4= b*z';

   Tc(:,i)= term1+term2+term3+term4;

end


% plot results and analytical solution 
figure;
plot(T(:,[1:50:300]),-z);hold on;
title(' Temperature (FD)');grid on;
xlabel('Temperature [K]');ylabel('Depth [m]');
xlim([-1 12]);
plot(Tc(1:10:nz,[1:50:300]),-z(1:10:nz),'o')
legend('t =  0 ka','  =  5 ka','  =  10 ka','  =  15 ka','  =  20 ka', ...
       '  =  25 ka',4)  

text(3,-2500, ...
 'Analytical Solution: Carslaw & Jaeger (1947)')
text(3,-3000, strcat('Theta = ',num2str(beta))) ;  

% plot residuals  
Tres= T(:,[1:50:300])-Tc(:,[1:50:300]);figure;plot(Tres,-z); grid on;    
legend('t =  0 ka','  =  5 ka','  =  10 ka','  =  15 ka','  =  20 ka', ...
       '  =  25 ka',3)  
title(' Temperature Residuals')
xlabel('Temperature [K]');ylabel('Depth [m]');% xlim([-0.03 0.03]);
% ylim([-1000 0]);
text(0.005,-600, strcat('Theta = ',num2str(beta))) ; 
 
figure;
subplot(1,2,1);
TT=T-Tc;imagesc(TT);colorbar   
title(' T-T_a');xlabel('timestep');ylabel('Depth [nodes]');
subplot(1,2,2);
TT=log10(abs(T-Tc)+eps);imagesc(TT);colorbar   
title(' log (|T-T_a|)');xlabel('timestep');
suptitle(strcat('Theta = ',num2str(beta)));


figure;
plot(T(:,[1:50:300]),-z,'o');hold on;
plot(Tc(1:10:nz,[1:50:300]),-z(1:10:nz),'+')
plot(T2(1:10:nz,[1:50:300]),-z(1:10:nz),'-')
title(' Temperatures');grid on;
xlabel('Temperature [K]');ylabel('Depth [m]');
xlim([-1 12]);
%legend('t =  0 ka','  =  5 ka','  =  10 ka','  =  15 ka','  =  20 ka', ...
 %      '  =  25 ka',3)  
